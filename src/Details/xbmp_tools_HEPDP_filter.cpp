namespace xbmp::tools::filters
{
    void gaussian_blur(xcore::bitmap& dst, const xcore::bitmap& src ) noexcept
    {
        constexpr auto kernel = std::array
        {
            std::array{0.003f, 0.013f, 0.022f, 0.013f, 0.003f},
            std::array{0.013f, 0.060f, 0.098f, 0.060f, 0.013f},
            std::array{0.022f, 0.098f, 0.162f, 0.098f, 0.022f},
            std::array{0.013f, 0.060f, 0.098f, 0.060f, 0.013f},
            std::array{0.003f, 0.013f, 0.022f, 0.013f, 0.003f}
        };

        const int width  = static_cast<int>(src.getWidth());
        const int height = static_cast<int>(src.getHeight());
        const int k_size = 5;
        const int k_half = k_size / 2;

        dst.CreateBitmap(width, height);

        auto SrcData = src.getMip<xcore::icolor>(0);
        auto DstData = dst.getMip<xcore::icolor>(0);

        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                xcore::fcolor AccPixel = { 0, 0, 0, 0 };
                for (int ky = -k_half; ky <= k_half; ++ky) 
                {
                    for (int kx = -k_half; kx <= k_half; ++kx) 
                    {
                        const int   sy      = std::max(0, std::min(height - 1, y + ky));
                        const int   sx      = std::max(0, std::min(width  - 1, x + kx));
                        const auto  pixel   = SrcData[sy * width + sx];
                        const float w       = kernel[ky + k_half][kx + k_half];

                        AccPixel += pixel * w;
                    }
                }

                auto& pixel = DstData[y * width + x];
                pixel = AccPixel;
                pixel.m_A = 255;
            }
        }
    }

    // Hybrid Edge Propagation with Detail Preservation(HEPDP)

    void MakeBitmapTilableHEPDP(xcore::bitmap& seamless, const xcore::bitmap& src_image, float WidthOverlapPercentage, float HeightOverlapPercentage) noexcept
    {
        const int width       = static_cast<int>(src_image.getWidth());
        const int height      = static_cast<int>(src_image.getHeight());
        const int blend_cols  = std::max(1, static_cast<int>(width * WidthOverlapPercentage));
        const int blend_rows  = std::max(1, static_cast<int>(height * HeightOverlapPercentage));

        //
        // Step 1: Decompose into base and detail layers
        //
        xcore::bitmap base;
        xcore::bitmap detail;

        detail.CreateBitmap(width, height);
        gaussian_blur(base, src_image);

        {
            auto src_image_pixel = src_image.getMip<xcore::icolor>(0);
            auto base_pixel      = base.getMip<xcore::icolor>(0);
            auto detail_pixel    = detail.getMip<xcore::icolor>(0);

            for (int y = 0; y < height; ++y) 
            {
                for (int x = 0; x < width; ++x) 
                {
                    auto& orig      = src_image_pixel[y * width + x];
                    auto& blurred   = base_pixel[y * width + x];
                    auto& detail    = detail_pixel[y * width + x];

                    detail.m_R = (orig.m_R - (int)blurred.m_R + 128);
                    detail.m_G = (orig.m_G - (int)blurred.m_G + 128);
                    detail.m_B = (orig.m_B - (int)blurred.m_B + 128);
                    detail.m_A = (orig.m_A - (int)blurred.m_A + 128);
                }
            }
        }

        //
        // Step 2 & 3: Horizontal and Vertical propagation on base layer
        //
        {
            const float     k = 1.14f;  // Decay rate

            // copy image
            seamless.CreateBitmap(width, height);

            auto seamless_pixels = seamless.getMip<xcore::icolor>(0);
            auto base_pixels     = base.getMip<xcore::icolor>(0);

            memcpy(seamless_pixels.data(), base_pixels.data(), base.getDataSize());

            for (int y = 0; y < height; ++y) 
            {
                for (int x = 0; x < blend_cols; ++x) 
                {
                    float alpha = pow( float(x) / blend_cols, 1.0f/k);

                    auto& L = base_pixels[y * width + x];
                    auto& R = base_pixels[y * width + (width - x - 1)];

                    auto& l_seamless_pixel = seamless_pixels[y * width + x];
                    auto& r_seamless_pixel = seamless_pixels[y * width + (width - x - 1)];

                    l_seamless_pixel.m_R = static_cast<std::uint8_t>((1 - alpha) * R.m_R + alpha * L.m_R);
                    l_seamless_pixel.m_G = static_cast<std::uint8_t>((1 - alpha) * R.m_G + alpha * L.m_G);
                    l_seamless_pixel.m_B = static_cast<std::uint8_t>((1 - alpha) * R.m_B + alpha * L.m_B);
                    l_seamless_pixel.m_A = static_cast<std::uint8_t>((1 - alpha) * R.m_A + alpha * L.m_A);

                    r_seamless_pixel.m_R = static_cast<std::uint8_t>((1 - alpha) * l_seamless_pixel.m_R + alpha * R.m_R);
                    r_seamless_pixel.m_G = static_cast<std::uint8_t>((1 - alpha) * l_seamless_pixel.m_G + alpha * R.m_G);
                    r_seamless_pixel.m_B = static_cast<std::uint8_t>((1 - alpha) * l_seamless_pixel.m_B + alpha * R.m_B);
                    r_seamless_pixel.m_A = static_cast<std::uint8_t>((1 - alpha) * l_seamless_pixel.m_A + alpha * R.m_A);
                }
            }

            for (int x = 0; x < width; ++x) 
            {
                for (int y = 0; y < blend_rows; ++y) 
                {
                    const float     alpha = pow(float(y) / blend_rows, 1.0f/k);

                    auto& T = seamless_pixels[y * width + x];
                    auto& B = seamless_pixels[x + (height - y - 1) * width];

                    auto& t_seamless_pixel = seamless_pixels[y * width + x];
                    auto& b_seamless_pixel = seamless_pixels[x + (height - y - 1) * width];

                    t_seamless_pixel.m_R = static_cast<std::uint8_t>((1 - alpha) * B.m_R + alpha * T.m_R);
                    t_seamless_pixel.m_G = static_cast<std::uint8_t>((1 - alpha) * B.m_G + alpha * T.m_G);
                    t_seamless_pixel.m_B = static_cast<std::uint8_t>((1 - alpha) * B.m_B + alpha * T.m_B);
                    t_seamless_pixel.m_A = static_cast<std::uint8_t>((1 - alpha) * B.m_A + alpha * T.m_A);

                    b_seamless_pixel.m_R = static_cast<std::uint8_t>((1 - alpha) * t_seamless_pixel.m_R + alpha * B.m_R);
                    b_seamless_pixel.m_G = static_cast<std::uint8_t>((1 - alpha) * t_seamless_pixel.m_G + alpha * B.m_G);
                    b_seamless_pixel.m_B = static_cast<std::uint8_t>((1 - alpha) * t_seamless_pixel.m_B + alpha * B.m_B);
                    b_seamless_pixel.m_A = static_cast<std::uint8_t>((1 - alpha) * t_seamless_pixel.m_A + alpha * B.m_A);
                }
            }
        }

        //
        // Step 4: Reintroduce details with tapering
        //
        {
            auto seamless_pixels = seamless.getMip<xcore::icolor>(0);
            auto detail_pixel    = detail.getMip<xcore::icolor>(0);

            for (int y = 0; y < height; ++y)
            {
                for (int x = 0; x < width; ++x) 
                {
                    const float dist_x = std::min(x, width - x - 1) / (float)blend_cols;
                    const float dist_y = std::min(y, height - y - 1) / (float)blend_rows;
                    const float taper  = std::min(1.0f, std::max(dist_x, dist_y));
                    auto& base_px      = seamless_pixels[y * width + x];
                    auto& detail_px    = detail_pixel[y * width + x];

                    xcore::icolor c; 
                    c.m_R = (uint8_t)std::min(255, std::max(0, (int)(base_px.m_R + (detail_px.m_R - 128) * taper)));
                    c.m_G = (uint8_t)std::min(255, std::max(0, (int)(base_px.m_G + (detail_px.m_G - 128) * taper)));
                    c.m_B = (uint8_t)std::min(255, std::max(0, (int)(base_px.m_B + (detail_px.m_B - 128) * taper)));
                    c.m_A = (uint8_t)std::min(255, std::max(0, (int)(base_px.m_A + (detail_px.m_A - 128) * taper)));

                    base_px = c;
                }
            }
        }
    }
} // namespace xbmp::tools::filters