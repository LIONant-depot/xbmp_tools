#include <omp.h>
#include <intrin.h>
#include <cstdint>

namespace xbmp::tools::mips
{

    // Modified Bessel functions of the first kind I_nu(x)
    // - bessel_i0(x): I_0(x)
    // - bessel_i1(x): I_1(x)
    // - cyl_bessel_i(nu, x): general order (nu >= 0)
    //
    // Notes:
    //   - For small and medium |x|, these use ascending power series.
    //   - For large |x|, I0 and I1 use simple asymptotic expansions.
    //   - cyl_bessel_i handles nu >= 0. For negative non-integer nu, a full
    //     implementation requires K_nu(x); this header does not provide K_nu.
    inline double eps() noexcept { return std::numeric_limits<double>::epsilon(); }
    inline double inf() noexcept { return std::numeric_limits<double>::infinity(); }
    inline double nan() noexcept { return std::numeric_limits<double>::quiet_NaN(); }

    // Use our own PI to avoid platform macros.
    static constexpr double PI = 3.141592653589793238462643383279502884;

    // I0(x): Modified Bessel of the first kind, order 0.
    // Series: I0(x) = sum_{k=0..inf} (x/2)^(2k) / (k!)^2
    inline double cyl_bessel_i0(double x)
    {
        const double ax = std::fabs(x);

        if (ax <= 8.0) {
            const double y = (x * x) * 0.25; // (x/2)^2
            double term = 1.0;               // k = 0
            double sum = 1.0;
            for (int k = 1; k < 200; ++k) {
                term *= y / (static_cast<double>(k) * static_cast<double>(k));
                sum += term;
                if (std::fabs(term) < std::fabs(sum) * eps() * 4.0) break;
            }
            return sum;
        }

        // Asymptotic: I0(x) ~ exp(x) / sqrt(2*pi*x) * (1 + 1/(8x) + 9/(128x^2) + 225/(3072x^3))
        const double invx = 1.0 / ax;
        const double pre = std::exp(ax) / std::sqrt(2.0 * PI * ax);
        const double p = 1.0 + invx * (1.0 / 8.0) + invx * invx * (9.0 / 128.0) + invx * invx * invx * (225.0 / 3072.0);
        return pre * p;
    }

    //---------------------------------------------------------------------------------------------
    // Kaiser windowed sinc filter implementation based on research
    // Beta ~4.0 for good balance between sharpness and ringing (from Godot proposal and general practices)
    // Window width: 3 lobes (support radius 3), for high quality without excessive cost
    //---------------------------------------------------------------------------------------------
    inline float kaiser(float x, float beta, float half_width)
    {
        if (std::abs(x) > half_width) return 0.0f;
        float t = x / half_width;
        float b = beta * std::sqrt(1.0f - t * t);
        return static_cast<float>(cyl_bessel_i0(b) / cyl_bessel_i0(beta)); // Modified Bessel function of order 0
    }

    constexpr int   kKaiserRadius = 3; // Support: -3 to 3, 7 taps
    constexpr float kKaiserBeta = 4.0f;
    constexpr float kKaiserHalfWidth = static_cast<float>(kKaiserRadius);

    // Precompute kernel weights for downsampling (separable, normalized)
    std::array<float, 2 * kKaiserRadius + 1> computeKaiserKernel()
    {
        std::array<float, 2 * kKaiserRadius + 1> kernel{};
        float sum = 0.0f;
        for (int i = -kKaiserRadius; i <= kKaiserRadius; ++i)
        {
            // Lanczos-like sinc for 0.5 freq cutoff
            float sinc = static_cast<float>((i == 0) ? 1.0f : std::sin(PI * i / 2.0f) / (PI * i / 2.0f));
            kernel[i + kKaiserRadius] = sinc * kaiser(static_cast<float>(i), kKaiserBeta, kKaiserHalfWidth);
            sum += kernel[i + kKaiserRadius];
        }

        for (auto& w : kernel) w /= sum; // Normalize
        return kernel;
    }

    static auto kKaiserKernel = computeKaiserKernel();

    constexpr int   kLanczosRadius = 3;
    constexpr float kLanczosA = 3.0f;

    std::array<float, 2 * kLanczosRadius + 1> computeLanczosKernel()
    {
        std::array<float, 2 * kLanczosRadius + 1> kernel{};
        float sum = 0.0f;
        for (int i = -kLanczosRadius; i <= kLanczosRadius; ++i)
        {
            float x = static_cast<float>(i);
            float sinc_val = static_cast<float>((i == 0) ? 1.0f : std::sin(PI * x / 2.0f) / (PI * x / 2.0f));
            float window = static_cast<float>((x == 0.0f) ? 1.0f : ((std::abs(x) < kLanczosA) ? std::sin(PI * x / kLanczosA) / (PI * x / kLanczosA) : 0.0f));
            kernel[i + kLanczosRadius] = sinc_val * window;
            sum += kernel[i + kLanczosRadius];
        }

        for (auto& w : kernel) w /= sum; // Normalize
        return kernel;
    }

    static auto kLanczosKernel = computeLanczosKernel();

    constexpr int   kTriangleRadius = 1;
    static constexpr std::array<float, 2 * kTriangleRadius + 1> kTriangleKernel = { 0.25f, 0.5f, 0.25f };

    // Helper to get wrapped coordinate
    int wrapCoord(int coord, int size, xbitmap::wrap_mode mode, bool force_clamp = false)
    {
        if (force_clamp) {
            return std::clamp(coord, 0, size - 1);
        }
        if (mode == xbitmap::wrap_mode::WRAP)
        {
            return (coord % size + size) % size;
        }
        else if (mode == xbitmap::wrap_mode::MIRROR)
        {
            int period = size * 2 - 2;
            coord = std::abs((coord % period + period) % period);
            if (coord >= size) coord = period - coord;
            return coord;
        }
        else
        { // CLAMP_TO_EDGE or others
            return std::clamp(coord, 0, size - 1);
        }
    }

    //---------------------------------------------------------------------------------------------
    // Main function: Generate high-quality mipmaps using Kaiser filter
    // - nMips = -1: full chain down to 1x1
    // - renormalizeNormals: If true, treat as tangent-space normal map (renormalize RGB vectors)
    // Based on:
    // - Kaiser for reduced blur/aliasing (Godot proposal, common in engines)
    // - Linear space filtering for sRGB
    // - Renormalization for normals (NVIDIA paper, Godot)
    // - Handles float/half for HDR
    // - Per-face for cubemaps
    // - Respects wrap modes in sampling
    //
    // Variance - based specular adjustment(e.g., Toksvig mapping) in mipmaps modifies roughness values in lower - resolution levels based on the variance(spread)
    // of normal vectors from the higher mip.This simulates lost sub - pixel detail, increasing effective roughness to reduce aliasing and overly shiny artifacts
    // on distant surfaces.For ARME (AO, Roughness, Metalness, Emisive) textures, it targets the roughness channel, making far - away materials appear duller and
    // more realistic, improving visual stability without extra runtime cost.High - variance normals(bumpy surfaces) get stronger adjustment.
    //---------------------------------------------------------------------------------------------
    void GenerateMipMaps(xbitmap& Bitmap, const int MinSize, bool isNormalMap, mips_filter_type filter, bool preserveAlphaCoverage, float alphaCutoff, const xbitmap* pNormalBitmap, int roughnessChannel)
    {
        if (!Bitmap.isValid() || Bitmap.getMipCount() > 1) return; // Assume starts with only mip 0

        const int orig_w = Bitmap.getWidth();
        const int orig_h = Bitmap.getHeight();
        if (orig_w == 0 || orig_h == 0) return;

        // Compute number of mips
        int max_mips = (MinSize <= 0) ? 1 + static_cast<int>(std::floor(std::log2(std::max(orig_w, orig_h)))) : ComputeMaxMips(MinSize, orig_w, orig_h);

        if (max_mips <= 1) return;

        const xbitmap::format fmt = Bitmap.getFormat();
        const bool  is_srgb = Bitmap.getColorSpace() == xbitmap::color_space::SRGB;
        const bool  renormalizeNormals = isNormalMap;

        const bool  is_float32 = (fmt == xbitmap::format::R32G32B32A32_FLOAT
            || fmt == xbitmap::format::R32G32B32_FLOAT
            || fmt == xbitmap::format::R32G32_FLOAT
            || fmt == xbitmap::format::R32_FLOAT
            );
        const bool  is_half = (fmt == xbitmap::format::R16G16B16A16_SFLOAT
            || fmt == xbitmap::format::R16G16_SFLOAT
            || fmt == xbitmap::format::R16_SFLOAT
            );
        const bool  is_int = (!is_float32 && !is_half); // Assume R8G8B8A8 or similar
        const int   channels = (fmt == xbitmap::format::R32G32B32A32_FLOAT || fmt == xbitmap::format::R16G16B16A16_SFLOAT)
            ? 4 : (fmt == xbitmap::format::R32G32B32_FLOAT)
            ? 3 : (fmt == xbitmap::format::R32G32_FLOAT || fmt == xbitmap::format::R16G16_SFLOAT)
            ? 2 : (fmt == xbitmap::format::R32_FLOAT || fmt == xbitmap::format::R16_SFLOAT)
            ? 1 : 4; // Default to 4 for int
        const int   bytes_per_pixel = is_int ? channels : (is_half ? channels * 2 : channels * 4);

        const int face_count = Bitmap.getFaceCount();
        const xbitmap::wrap_mode u_wrap = Bitmap.getUWrapMode();
        const xbitmap::wrap_mode v_wrap = Bitmap.getVWrapMode();
        const bool is_cube = Bitmap.isCubemap();

        bool use_toksvig = (pNormalBitmap != nullptr) && (channels >= 4); // Assume ARME is 4 channels

        // Assume normal bitmap has same properties, but 3 channels effective (RG for xy, reconstruct z)
        xbitmap::format normal_fmt;
        int normal_channels = 3; // Always process as xyz
        int normal_bytes_per_pixel = 0;
        bool is_normal_float32 = false;
        bool is_normal_half = false;
        bool is_normal_int = false;
        if (use_toksvig) {
            normal_fmt = pNormalBitmap->getFormat();
            is_normal_float32 = (normal_fmt == xbitmap::format::R32G32B32A32_FLOAT ||
                normal_fmt == xbitmap::format::R32G32B32_FLOAT ||
                normal_fmt == xbitmap::format::R32G32_FLOAT ||
                normal_fmt == xbitmap::format::R32_FLOAT);
            is_normal_half = (normal_fmt == xbitmap::format::R16G16B16A16_SFLOAT ||
                normal_fmt == xbitmap::format::R16G16_SFLOAT ||
                normal_fmt == xbitmap::format::R16_SFLOAT);
            is_normal_int = !is_normal_float32 && !is_normal_half;
            normal_bytes_per_pixel = static_cast<int>(pNormalBitmap->getMip<std::byte>(0).size() / (orig_w * orig_h));
        }

        // Compute mip sizes
        std::vector<std::uint64_t> mip_sizes(max_mips);
        int w = orig_w, h = orig_h;
        for (int level = 0; level < max_mips; ++level)
        {
            mip_sizes[level] = static_cast<std::uint64_t>(w) * h * bytes_per_pixel;
            w = std::max(1, w / 2);
            h = std::max(1, h / 2);
        }

        // Total face size (sum of all mips for one face)
        std::uint64_t face_size = 0;
        for (int level = 0; level < max_mips; ++level) face_size += mip_sizes[level];

        // Mip table size
        const auto MipTableSize = sizeof(xbitmap::mip) * max_mips;

        // Total data size (mip table + data)
        std::uint64_t total_size = MipTableSize + face_size * face_count;

        // Allocate new data
        auto new_data = std::make_unique<std::byte[]>(total_size);

        // Set mip table
        auto* p_mip_table = reinterpret_cast<xbitmap::mip*>(new_data.get());
        int cur_offset = 0;
        for (int level = 0; level < max_mips; ++level)
        {
            p_mip_table[level].m_Offset = cur_offset;
            cur_offset += static_cast<int>(mip_sizes[level]);
        }

        // Data start
        auto* p_data = new_data.get() + MipTableSize;

        // Copy level 0 for all faces
        for (int face = 0; face < face_count; ++face)
        {
            auto src_span = Bitmap.getMip<std::byte>(0, face);
            std::memcpy(p_data + face * face_size + p_mip_table[0].m_Offset, src_span.data(), mip_sizes[0]);
        }

        // Re-setup bitmap with new data and mip count
        auto old_bitmap_local = std::move(Bitmap); // Ensure old data freed
        Bitmap.setup(orig_w, orig_h, fmt, face_size, { new_data.release(), total_size }, true, max_mips, 1, old_bitmap_local.isCubemap());
        Bitmap.setUWrapMode(u_wrap);
        Bitmap.setVWrapMode(v_wrap);
        Bitmap.setColorSpace(old_bitmap_local.getColorSpace());

        // Select kernel and radius
        std::span<const float> kernel_span;
        int radius = 0;
        bool is_box = false;
        switch (filter)
        {
        case mips_filter_type::BOX:
            is_box = true;
            break;
        case mips_filter_type::TRIANGLE:
            radius = kTriangleRadius;
            kernel_span = kTriangleKernel;
            break;
        case mips_filter_type::LANCZOS:
            radius = kLanczosRadius;
            kernel_span = kLanczosKernel;
            break;
        case mips_filter_type::KAISER:
            radius = kKaiserRadius;
            kernel_span = kKaiserKernel;
            break;
        }

        const bool use_log_filter = (is_float32 || is_half) && !renormalizeNormals && (filter == mips_filter_type::LANCZOS || filter == mips_filter_type::KAISER);
        const int log_channels = (preserveAlphaCoverage && channels >= 4) ? channels - 1 : channels;

        // Generate lower mips
        for (int level = 1; level < max_mips; ++level)
        {
            const int prev_w = std::max(1, orig_w >> (level - 1));
            const int prev_h = std::max(1, orig_h >> (level - 1));
            const int cur_w = std::max(1, prev_w / 2);
            const int cur_h = std::max(1, prev_h / 2);

#pragma omp parallel for
            for (int face = 0; face < face_count; ++face)
            {
                auto prev_span = Bitmap.getMip<std::byte>(level - 1, face);
                auto cur_span = Bitmap.getMip<std::byte>(level, face);

                // Temp buffer for horizontal pass (float accum for quality)
                std::vector<float> horiz(static_cast<std::size_t>(cur_w) * prev_h * channels, 0.0f);

                std::vector<float> horiz_normal;
                if (use_toksvig) {
                    horiz_normal.resize(static_cast<std::size_t>(cur_w) * prev_h * normal_channels, 0.0f);
                }

                // Helper to sample prev mip with wrap, returns array of channels
                auto sample = [&](int x, int y) -> std::array<float, 4>
                    {
                        x = wrapCoord(x, prev_w, u_wrap, is_cube);
                        y = wrapCoord(y, prev_h, v_wrap, is_cube);
                        std::array<float, 4> color{ 0.0f, 0.0f, 0.0f, 1.0f }; // Default A=1

                        const std::byte* ptr = &prev_span[(y * prev_w + x) * bytes_per_pixel];
                        if (is_int)
                        {
                            const xcolori& pixel = *reinterpret_cast<const xcolori*>(ptr);
                            color[0] = pixel.m_R / 255.0f;
                            color[1] = pixel.m_G / 255.0f;
                            color[2] = pixel.m_B / 255.0f;
                            if (channels >= 4) color[3] = pixel.m_A / 255.0f;
                            if (is_srgb && !renormalizeNormals)
                            {
                                color[0] = std::pow(color[0], 2.2f);
                                color[1] = std::pow(color[1], 2.2f);
                                color[2] = std::pow(color[2], 2.2f);
                                // A linear
                            }
                        }
                        else if (is_half)
                        {
                            const uint16_t* hptr = reinterpret_cast<const uint16_t*>(ptr);
                            xmath::half h;
                            for (int c = 0; c < channels; ++c)
                            {
                                h.m_Value = hptr[c];
                                color[c] = static_cast<float>(h);
                            }
                        }
                        else
                        { // float32
                            const float* fptr = reinterpret_cast<const float*>(ptr);
                            for (int c = 0; c < channels; ++c)
                            {
                                color[c] = fptr[c];
                            }
                        }
                        if ((is_float32 || is_half) && !renormalizeNormals) {
                            for (int c = 0; c < channels; ++c) {
                                if (std::isnan(color[c])) {
                                    color[c] = 0.0f;
                                }
                                else if (std::isinf(color[c])) {
                                    color[c] = is_half ? 65504.0f : 3.4028235e+38f; // FLT_MAX
                                }
                                else {
                                    color[c] = std::max(0.0f, color[c]);
                                }
                            }
                        }
                        if (preserveAlphaCoverage && channels >= 4) {
                            color[3] = (color[3] > alphaCutoff) ? 1.0f : 0.0f;
                        }
                        if (use_log_filter) {
                            for (int c = 0; c < log_channels; ++c) {
                                color[c] = std::log1p(color[c]);
                            }
                        }
                        return color;
                    };

                // Helper to sample normal map
                auto sample_normal = [&](int x, int y) -> std::array<float, 3>
                    {
                        x = wrapCoord(x, prev_w, u_wrap, is_cube);
                        y = wrapCoord(y, prev_h, v_wrap, is_cube);
                        std::array<float, 3> n{ 0.0f, 0.0f, 1.0f };

                        auto normal_span = pNormalBitmap->getMip<std::byte>(level - 1, face);
                        const std::byte* ptr = &normal_span[(y * prev_w + x) * normal_bytes_per_pixel];
                        float rx, gy;
                        if (is_normal_int) {
                            const xcolori& pixel = *reinterpret_cast<const xcolori*>(ptr);
                            rx = pixel.m_R / 255.0f;
                            gy = pixel.m_G / 255.0f;
                        }
                        else if (is_normal_half) {
                            const uint16_t* hptr = reinterpret_cast<const uint16_t*>(ptr);
                            xmath::half h;
                            h.m_Value = hptr[0];
                            rx = static_cast<float>(h);
                            h.m_Value = hptr[1];
                            gy = static_cast<float>(h);
                        }
                        else {
                            const float* fptr = reinterpret_cast<const float*>(ptr);
                            rx = fptr[0];
                            gy = fptr[1];
                        }
                        float normal_unpack_scale = is_normal_int ? 2.0f : 1.0f;
                        float normal_unpack_bias = is_normal_int ? -1.0f : 0.0f;
                        float ix = rx * normal_unpack_scale + normal_unpack_bias;
                        float iy = gy * normal_unpack_scale + normal_unpack_bias;
                        float iz = std::sqrt(std::max(0.0f, 1.0f - ix * ix - iy * iy)); // Assume positive z
                        n[0] = ix;
                        n[1] = iy;
                        n[2] = iz;
                        return n;
                    };

                // Horizontal pass
#pragma omp parallel for
                for (int y = 0; y < prev_h; ++y)
                {
                    for (int x = 0; x < cur_w; ++x)
                    {
                        const int sx = x * 2;
                        float* hptr = &horiz[(y * cur_w + x) * channels];
                        if (is_box)
                        {
                            auto c0 = sample(sx, y);
                            auto c1 = sample(sx + 1, y);
                            for (int c = 0; c < channels; ++c)
                            {
                                hptr[c] = (c0[c] + c1[c]) * 0.5f;
                            }
                        }
                        else
                        {
                            std::array<float, 4> accum{ 0.0f, 0.0f, 0.0f, 0.0f };
                            for (int tap = -radius; tap <= radius; ++tap)
                            {
                                float weight = kernel_span[tap + radius];
                                auto scolor = sample(sx + tap, y);
                                for (int c = 0; c < channels; ++c)
                                {
                                    accum[c] += weight * scolor[c];
                                }
                            }
                            for (int c = 0; c < channels; ++c)
                            {
                                hptr[c] = accum[c];
                            }
                        }

                        if (use_toksvig) {
                            float* hnptr = &horiz_normal[(y * cur_w + x) * normal_channels];
                            if (is_box) {
                                auto n0 = sample_normal(sx, y);
                                auto n1 = sample_normal(sx + 1, y);
                                for (int c = 0; c < normal_channels; ++c) {
                                    hnptr[c] = (n0[c] + n1[c]) * 0.5f;
                                }
                            }
                            else {
                                std::array<float, 3> accum_n{ 0.0f, 0.0f, 0.0f };
                                for (int tap = -radius; tap <= radius; ++tap)
                                {
                                    float weight = kernel_span[tap + radius];
                                    auto sn = sample_normal(sx + tap, y);
                                    for (int c = 0; c < normal_channels; ++c) {
                                        accum_n[c] += weight * sn[c];
                                    }
                                }
                                for (int c = 0; c < normal_channels; ++c) {
                                    hnptr[c] = accum_n[c];
                                }
                            }
                        }
                    }
                }

                // Vertical pass + output
#pragma omp parallel for
                for (int y = 0; y < cur_h; ++y)
                {
                    for (int x = 0; x < cur_w; ++x)
                    {
                        const int sy = y * 2;
                        std::array<float, 4> accum{ 0.0f, 0.0f, 0.0f, 0.0f };
                        if (is_box)
                        {
                            auto c0 = std::array<float, 4>{};
                            auto c1 = std::array<float, 4>{};
                            int wrapped_sy = wrapCoord(sy, prev_h, v_wrap, is_cube);
                            int wrapped_sy1 = wrapCoord(sy + 1, prev_h, v_wrap, is_cube);
                            const float* h0 = &horiz[(wrapped_sy * cur_w + x) * channels];
                            const float* h1 = &horiz[(wrapped_sy1 * cur_w + x) * channels];
                            for (int c = 0; c < channels; ++c)
                            {
                                c0[c] = h0[c];
                                c1[c] = h1[c];
                            }
                            for (int c = 0; c < channels; ++c)
                            {
                                accum[c] = (c0[c] + c1[c]) * 0.5f;
                            }
                        }
                        else
                        {
                            for (int tap = -radius; tap <= radius; ++tap)
                            {
                                float weight = kernel_span[tap + radius];
                                int wrapped_y = wrapCoord(sy + tap, prev_h, v_wrap, is_cube);
                                const float* hptr = &horiz[(wrapped_y * cur_w + x) * channels];
                                for (int c = 0; c < channels; ++c)
                                {
                                    accum[c] += weight * hptr[c];
                                }
                            }
                        }

                        if (use_log_filter) {
                            for (int c = 0; c < log_channels; ++c) {
                                accum[c] = std::max(accum[c], 0.0f);
                                accum[c] = std::expm1(accum[c]);
                            }
                        }

                        std::array<float, 3> accum_n{ 0.0f, 0.0f, 0.0f };
                        if (use_toksvig) {
                            if (is_box) {
                                auto n0 = std::array<float, 3>{};
                                auto n1 = std::array<float, 3>{};
                                int wrapped_sy = wrapCoord(sy, prev_h, v_wrap, is_cube);
                                int wrapped_sy1 = wrapCoord(sy + 1, prev_h, v_wrap, is_cube);
                                const float* hn0 = &horiz_normal[(wrapped_sy * cur_w + x) * normal_channels];
                                const float* hn1 = &horiz_normal[(wrapped_sy1 * cur_w + x) * normal_channels];
                                for (int c = 0; c < normal_channels; ++c) {
                                    n0[c] = hn0[c];
                                    n1[c] = hn1[c];
                                }
                                for (int c = 0; c < normal_channels; ++c) {
                                    accum_n[c] = (n0[c] + n1[c]) * 0.5f;
                                }
                            }
                            else {
                                for (int tap = -radius; tap <= radius; ++tap) {
                                    float weight = kernel_span[tap + radius];
                                    int wrapped_y = wrapCoord(sy + tap, prev_h, v_wrap, is_cube);
                                    const float* hnptr = &horiz_normal[(wrapped_y * cur_w + x) * normal_channels];
                                    for (int c = 0; c < normal_channels; ++c) {
                                        accum_n[c] += weight * hnptr[c];
                                    }
                                }
                            }

                            float len_sq = accum_n[0] * accum_n[0] + accum_n[1] * accum_n[1] + accum_n[2] * accum_n[2];
                            float len = std::sqrt(len_sq);
                            if (len > 0.0001f) {
                                float variance = 1.0f - len_sq; // Approx variance = 1 - Ft^2, Ft = len
                                accum[roughnessChannel] = std::sqrt(accum[roughnessChannel] * accum[roughnessChannel] + variance);
                            }
                        }

                        // Renormalize for normals (assume RGB normal, A unused or height)
                        if (renormalizeNormals && channels >= 3)
                        {
                            float main_unpack_scale = is_int ? 2.0f : 1.0f;
                            float main_unpack_bias = is_int ? -1.0f : 0.0f;
                            xmath::fvec3 n(accum[0] * main_unpack_scale + main_unpack_bias,
                                accum[1] * main_unpack_scale + main_unpack_bias,
                                accum[2] * main_unpack_scale + main_unpack_bias);
                            n = n.NormalizeSafeCopy();
                            float main_pack_scale = is_int ? 0.5f : 1.0f;
                            float main_pack_bias = is_int ? 0.5f : 0.0f;
                            accum[0] = n[0] * main_pack_scale + main_pack_bias;
                            accum[1] = n[1] * main_pack_scale + main_pack_bias;
                            accum[2] = n[2] * main_pack_scale + main_pack_bias;
                            // A unchanged
                        }

                        // Clamp for HDR to prevent negatives/overflow from ringing
                        if (!renormalizeNormals && (is_float32 || is_half))
                        {
                            for (int c = 0; c < channels; ++c)
                            {
                                accum[c] = std::max(0.0f, accum[c]);
                                if (is_half)
                                {
                                    accum[c] = std::min(accum[c], 65504.0f);
                                }
                            }
                        }

                        // Store
                        std::byte* optr = &cur_span[(y * cur_w + x) * bytes_per_pixel];
                        if (is_int)
                        {
                            float r = accum[0];
                            float g = accum[1];
                            float b = accum[2];
                            float a = (channels >= 4) ? accum[3] : 1.0f;
                            if (is_srgb && !renormalizeNormals)
                            {
                                r = std::pow(std::clamp(r, 0.0f, 1.0f), 1.0f / 2.2f);
                                g = std::pow(std::clamp(g, 0.0f, 1.0f), 1.0f / 2.2f);
                                b = std::pow(std::clamp(b, 0.0f, 1.0f), 1.0f / 2.2f);
                                // A linear
                            }
                            xcolori pixel;
                            pixel.m_R = static_cast<std::uint8_t>(std::clamp(r * 255.0f + 0.5f, 0.0f, 255.0f));
                            pixel.m_G = static_cast<std::uint8_t>(std::clamp(g * 255.0f + 0.5f, 0.0f, 255.0f));
                            pixel.m_B = static_cast<std::uint8_t>(std::clamp(b * 255.0f + 0.5f, 0.0f, 255.0f));
                            pixel.m_A = static_cast<std::uint8_t>(std::clamp(a * 255.0f + 0.5f, 0.0f, 255.0f));
                            std::memcpy(optr, &pixel, bytes_per_pixel);
                        }
                        else if (is_half)
                        {
                            uint16_t* hptr = reinterpret_cast<uint16_t*>(optr);
                            for (int c = 0; c < channels; ++c)
                            {
                                hptr[c] = xmath::half(accum[c]).m_Value;
                            }
                        }
                        else
                        { // float32
                            float* fptr = reinterpret_cast<float*>(optr);
                            for (int c = 0; c < channels; ++c)
                            {
                                fptr[c] = accum[c];
                            }
                        }
                    }
                }
            }
        }
    }
}








/////////////////////////////////// FIRST FULLY WORKING VERSION AFTER SECOND ROUND....
/*
namespace xbmp::tools::mips
{

    // Modified Bessel functions of the first kind I_nu(x)
    // - bessel_i0(x): I_0(x)
    // - bessel_i1(x): I_1(x)
    // - cyl_bessel_i(nu, x): general order (nu >= 0)
    //
    // Notes:
    //   - For small and medium |x|, these use ascending power series.
    //   - For large |x|, I0 and I1 use simple asymptotic expansions.
    //   - cyl_bessel_i handles nu >= 0. For negative non-integer nu, a full
    //     implementation requires K_nu(x); this header does not provide K_nu.
    inline double eps() noexcept { return std::numeric_limits<double>::epsilon(); }
    inline double inf() noexcept { return std::numeric_limits<double>::infinity(); }
    inline double nan() noexcept { return std::numeric_limits<double>::quiet_NaN(); }

    // Use our own PI to avoid platform macros.
    static constexpr double PI = 3.141592653589793238462643383279502884;

    // I0(x): Modified Bessel of the first kind, order 0.
    // Series: I0(x) = sum_{k=0..inf} (x/2)^(2k) / (k!)^2
    inline double cyl_bessel_i0(double x)
    {
        const double ax = std::fabs(x);

        if (ax <= 8.0) {
            const double y = (x * x) * 0.25; // (x/2)^2
            double term = 1.0;               // k = 0
            double sum = 1.0;
            for (int k = 1; k < 200; ++k) {
                term *= y / (static_cast<double>(k) * static_cast<double>(k));
                sum += term;
                if (std::fabs(term) < std::fabs(sum) * eps() * 4.0) break;
            }
            return sum;
        }

        // Asymptotic: I0(x) ~ exp(x) / sqrt(2*pi*x) * (1 + 1/(8x) + 9/(128x^2) + 225/(3072x^3))
        const double invx = 1.0 / ax;
        const double pre = std::exp(ax) / std::sqrt(2.0 * PI * ax);
        const double p = 1.0 + invx * (1.0 / 8.0) + invx * invx * (9.0 / 128.0) + invx * invx * invx * (225.0 / 3072.0);
        return pre * p;
    }

    //---------------------------------------------------------------------------------------------
    // Kaiser windowed sinc filter implementation based on research
    // Beta ~4.0 for good balance between sharpness and ringing (from Godot proposal and general practices)
    // Window width: 3 lobes (support radius 3), for high quality without excessive cost
    //---------------------------------------------------------------------------------------------
    inline float kaiser(float x, float beta, float half_width)
    {
        if (std::abs(x) > half_width) return 0.0f;
        float t = x / half_width;
        float b = beta * std::sqrt(1.0f - t * t);
        return static_cast<float>(cyl_bessel_i0(b) / cyl_bessel_i0(beta)); // Modified Bessel function of order 0
    }

    constexpr int   kKaiserRadius = 3; // Support: -3 to 3, 7 taps
    constexpr float kKaiserBeta = 4.0f;
    constexpr float kKaiserHalfWidth = static_cast<float>(kKaiserRadius);

    // Precompute kernel weights for downsampling (separable, normalized)
    std::array<float, 2 * kKaiserRadius + 1> computeKaiserKernel()
    {
        std::array<float, 2 * kKaiserRadius + 1> kernel{};
        float sum = 0.0f;
        for (int i = -kKaiserRadius; i <= kKaiserRadius; ++i)
        {
            // Lanczos-like sinc for 0.5 freq cutoff
            float sinc = static_cast<float>((i == 0) ? 1.0f : std::sin(PI * i / 2.0f) / (PI * i / 2.0f));
            kernel[i + kKaiserRadius] = sinc * kaiser(static_cast<float>(i), kKaiserBeta, kKaiserHalfWidth);
            sum += kernel[i + kKaiserRadius];
        }

        for (auto& w : kernel) w /= sum; // Normalize
        return kernel;
    }

    static auto kKaiserKernel = computeKaiserKernel();

    constexpr int   kLanczosRadius = 3;
    constexpr float kLanczosA = 3.0f;

    std::array<float, 2 * kLanczosRadius + 1> computeLanczosKernel()
    {
        std::array<float, 2 * kLanczosRadius + 1> kernel{};
        float sum = 0.0f;
        for (int i = -kLanczosRadius; i <= kLanczosRadius; ++i)
        {
            float x = static_cast<float>(i);
            float sinc_val = static_cast<float>((i == 0) ? 1.0f : std::sin(PI * x / 2.0f) / (PI * x / 2.0f));
            float window = (x == 0.0f) ? 1.0f : ((std::abs(x) < kLanczosA) ? std::sin(PI * x / kLanczosA) / (PI * x / kLanczosA) : 0.0f);
            kernel[i + kLanczosRadius] = sinc_val * window;
            sum += kernel[i + kLanczosRadius];
        }

        for (auto& w : kernel) w /= sum; // Normalize
        return kernel;
    }

    static auto kLanczosKernel = computeLanczosKernel();

    constexpr int   kTriangleRadius = 1;
    static constexpr std::array<float, 2 * kTriangleRadius + 1> kTriangleKernel = { 0.25f, 0.5f, 0.25f };

    // Helper to get wrapped coordinate
    int wrapCoord(int coord, int size, xbitmap::wrap_mode mode, bool force_clamp = false)
    {
        if (force_clamp) {
            return std::clamp(coord, 0, size - 1);
        }
        if (mode == xbitmap::wrap_mode::WRAP)
        {
            return (coord % size + size) % size;
        }
        else if (mode == xbitmap::wrap_mode::MIRROR)
        {
            int period = size * 2 - 2;
            coord = std::abs((coord % period + period) % period);
            if (coord >= size) coord = period - coord;
            return coord;
        }
        else
        { // CLAMP_TO_EDGE or others
            return std::clamp(coord, 0, size - 1);
        }
    }

    //---------------------------------------------------------------------------------------------
    // Main function: Generate high-quality mipmaps using Kaiser filter
    // - nMips = -1: full chain down to 1x1
    // - renormalizeNormals: If true, treat as tangent-space normal map (renormalize RGB vectors)
    // Based on:
    // - Kaiser for reduced blur/aliasing (Godot proposal, common in engines)
    // - Linear space filtering for sRGB
    // - Renormalization for normals (NVIDIA paper, Godot)
    // - Handles float/half for HDR
    // - Per-face for cubemaps
    // - Respects wrap modes in sampling
    //
    // Variance - based specular adjustment(e.g., Toksvig mapping) in mipmaps modifies roughness values in lower - resolution levels based on the variance(spread)
    // of normal vectors from the higher mip.This simulates lost sub - pixel detail, increasing effective roughness to reduce aliasing and overly shiny artifacts
    // on distant surfaces.For ARME (AO, Roughness, Metalness, Emisive) textures, it targets the roughness channel, making far - away materials appear duller and
    // more realistic, improving visual stability without extra runtime cost.High - variance normals(bumpy surfaces) get stronger adjustment.
    //---------------------------------------------------------------------------------------------
    void GenerateMipMaps(xbitmap& Bitmap, const int MinSize, bool isNormalMap, mips_filter_type filter, bool preserveAlphaCoverage, float alphaCutoff, const xbitmap* pNormalBitmap, int roughnessChannel)
    {
        if (!Bitmap.isValid() || Bitmap.getMipCount() > 1) return; // Assume starts with only mip 0

        const int orig_w = Bitmap.getWidth();
        const int orig_h = Bitmap.getHeight();
        if (orig_w == 0 || orig_h == 0) return;

        // Compute number of mips
        int max_mips = (MinSize <= 0) ? 1 + static_cast<int>(std::floor(std::log2(std::max(orig_w, orig_h)))) : ComputeMaxMips(MinSize, orig_w, orig_h);

        if (max_mips <= 1) return;

        const xbitmap::format fmt = Bitmap.getFormat();
        const bool  is_srgb = Bitmap.getColorSpace() == xbitmap::color_space::SRGB;
        const bool  renormalizeNormals = isNormalMap;

        const bool  is_float32 = (fmt == xbitmap::format::R32G32B32A32_FLOAT
            || fmt == xbitmap::format::R32G32B32_FLOAT
            || fmt == xbitmap::format::R32G32_FLOAT
            || fmt == xbitmap::format::R32_FLOAT
            );
        const bool  is_half = (fmt == xbitmap::format::R16G16B16A16_SFLOAT
            || fmt == xbitmap::format::R16G16_SFLOAT
            || fmt == xbitmap::format::R16_SFLOAT
            );
        const bool  is_int = (!is_float32 && !is_half); // Assume R8G8B8A8 or similar
        const int   channels = (fmt == xbitmap::format::R32G32B32A32_FLOAT || fmt == xbitmap::format::R16G16B16A16_SFLOAT)
            ? 4 : (fmt == xbitmap::format::R32G32B32_FLOAT)
            ? 3 : (fmt == xbitmap::format::R32G32_FLOAT || fmt == xbitmap::format::R16G16_SFLOAT)
            ? 2 : (fmt == xbitmap::format::R32_FLOAT || fmt == xbitmap::format::R16_SFLOAT)
            ? 1 : 4; // Default to 4 for int
        const int   bytes_per_pixel = is_int ? channels : (is_half ? channels * 2 : channels * 4);

        const int face_count = Bitmap.getFaceCount();
        const xbitmap::wrap_mode u_wrap = Bitmap.getUWrapMode();
        const xbitmap::wrap_mode v_wrap = Bitmap.getVWrapMode();
        const bool is_cube = Bitmap.isCubemap();

        bool use_toksvig = (pNormalBitmap != nullptr) && (channels >= 4); // Assume ARME is 4 channels

        // Assume normal bitmap has same properties, but 3 channels effective (RG for xy, reconstruct z)
        xbitmap::format normal_fmt;
        int normal_channels = 3; // Always process as xyz
        int normal_bytes_per_pixel = 0;
        bool is_normal_float32 = false;
        bool is_normal_half = false;
        bool is_normal_int = false;
        if (use_toksvig) {
            normal_fmt = pNormalBitmap->getFormat();
            is_normal_float32 = (normal_fmt == xbitmap::format::R32G32B32A32_FLOAT ||
                normal_fmt == xbitmap::format::R32G32B32_FLOAT ||
                normal_fmt == xbitmap::format::R32G32_FLOAT ||
                normal_fmt == xbitmap::format::R32_FLOAT);
            is_normal_half = (normal_fmt == xbitmap::format::R16G16B16A16_SFLOAT ||
                normal_fmt == xbitmap::format::R16G16_SFLOAT ||
                normal_fmt == xbitmap::format::R16_SFLOAT);
            is_normal_int = !is_normal_float32 && !is_normal_half;
            normal_bytes_per_pixel = static_cast<int>(pNormalBitmap->getMip<std::byte>(0).size() / (orig_w * orig_h));
        }

        // Compute mip sizes
        std::vector<std::uint64_t> mip_sizes(max_mips);
        int w = orig_w, h = orig_h;
        for (int level = 0; level < max_mips; ++level)
        {
            mip_sizes[level] = static_cast<std::uint64_t>(w) * h * bytes_per_pixel;
            w = std::max(1, w / 2);
            h = std::max(1, h / 2);
        }

        // Total face size (sum of all mips for one face)
        std::uint64_t face_size = 0;
        for (int level = 0; level < max_mips; ++level) face_size += mip_sizes[level];

        // Mip table size
        const auto MipTableSize = sizeof(xbitmap::mip) * max_mips;

        // Total data size (mip table + data)
        std::uint64_t total_size = MipTableSize + face_size * face_count;

        // Allocate new data
        auto new_data = std::make_unique<std::byte[]>(total_size);

        // Set mip table
        auto* p_mip_table = reinterpret_cast<xbitmap::mip*>(new_data.get());
        int cur_offset = 0;
        for (int level = 0; level < max_mips; ++level)
        {
            p_mip_table[level].m_Offset = cur_offset;
            cur_offset += static_cast<int>(mip_sizes[level]);
        }

        // Data start
        auto* p_data = new_data.get() + MipTableSize;

        // Copy level 0 for all faces
        for (int face = 0; face < face_count; ++face)
        {
            auto src_span = Bitmap.getMip<std::byte>(0, face);
            std::memcpy(p_data + face * face_size + p_mip_table[0].m_Offset, src_span.data(), mip_sizes[0]);
        }

        // Re-setup bitmap with new data and mip count
        auto old_bitmap_local = std::move(Bitmap); // Ensure old data freed
        Bitmap.setup(orig_w, orig_h, fmt, face_size, { new_data.release(), total_size }, true, max_mips, 1, old_bitmap_local.isCubemap());
        Bitmap.setUWrapMode(u_wrap);
        Bitmap.setVWrapMode(v_wrap);
        Bitmap.setColorSpace(old_bitmap_local.getColorSpace());

        // Select kernel and radius
        std::span<const float> kernel_span;
        int radius = 0;
        bool is_box = false;
        switch (filter)
        {
        case mips_filter_type::BOX:
            is_box = true;
            break;
        case mips_filter_type::TRIANGLE:
            radius = kTriangleRadius;
            kernel_span = kTriangleKernel;
            break;
        case mips_filter_type::LANCZOS:
            radius = kLanczosRadius;
            kernel_span = kLanczosKernel;
            break;
        case mips_filter_type::KAISER:
            radius = kKaiserRadius;
            kernel_span = kKaiserKernel;
            break;
        }

        const bool use_log_filter = (is_float32 || is_half) && !renormalizeNormals && (filter == mips_filter_type::LANCZOS || filter == mips_filter_type::KAISER);
        const int log_channels = (preserveAlphaCoverage && channels >= 4) ? channels - 1 : channels;

        // Generate lower mips
        for (int level = 1; level < max_mips; ++level)
        {
            const int prev_w = std::max(1, orig_w >> (level - 1));
            const int prev_h = std::max(1, orig_h >> (level - 1));
            const int cur_w = std::max(1, prev_w / 2);
            const int cur_h = std::max(1, prev_h / 2);

            for (int face = 0; face < face_count; ++face)
            {
                auto prev_span = Bitmap.getMip<std::byte>(level - 1, face);
                auto cur_span = Bitmap.getMip<std::byte>(level, face);

                // Temp buffer for horizontal pass (float accum for quality)
                std::vector<float> horiz(static_cast<std::size_t>(cur_w) * prev_h * channels, 0.0f);

                std::vector<float> horiz_normal;
                if (use_toksvig) {
                    horiz_normal.resize(static_cast<std::size_t>(cur_w) * prev_h * normal_channels, 0.0f);
                }

                // Helper to sample prev mip with wrap, returns array of channels
                auto sample = [&](int x, int y) -> std::array<float, 4>
                    {
                        x = wrapCoord(x, prev_w, u_wrap, is_cube);
                        y = wrapCoord(y, prev_h, v_wrap, is_cube);
                        std::array<float, 4> color{ 0.0f, 0.0f, 0.0f, 1.0f }; // Default A=1

                        const std::byte* ptr = &prev_span[(y * prev_w + x) * bytes_per_pixel];
                        if (is_int)
                        {
                            const xcolori& pixel = *reinterpret_cast<const xcolori*>(ptr);
                            color[0] = pixel.m_R / 255.0f;
                            color[1] = pixel.m_G / 255.0f;
                            color[2] = pixel.m_B / 255.0f;
                            if (channels >= 4) color[3] = pixel.m_A / 255.0f;
                            if (is_srgb && !renormalizeNormals)
                            {
                                color[0] = std::pow(color[0], 2.2f);
                                color[1] = std::pow(color[1], 2.2f);
                                color[2] = std::pow(color[2], 2.2f);
                                // A linear
                            }
                        }
                        else if (is_half)
                        {
                            const uint16_t* hptr = reinterpret_cast<const uint16_t*>(ptr);
                            xmath::half h;
                            for (int c = 0; c < channels; ++c)
                            {
                                h.m_Value = hptr[c];
                                color[c] = static_cast<float>(h);
                            }
                        }
                        else
                        { // float32
                            const float* fptr = reinterpret_cast<const float*>(ptr);
                            for (int c = 0; c < channels; ++c)
                            {
                                color[c] = fptr[c];
                            }
                        }
                        if ((is_float32 || is_half) && !renormalizeNormals) {
                            for (int c = 0; c < channels; ++c) {
                                if (std::isnan(color[c])) {
                                    color[c] = 0.0f;
                                }
                                else if (std::isinf(color[c])) {
                                    color[c] = is_half ? 65504.0f : 3.4028235e+38f; // FLT_MAX
                                }
                                else {
                                    color[c] = std::max(0.0f, color[c]);
                                }
                            }
                        }
                        if (preserveAlphaCoverage && channels >= 4) {
                            color[3] = (color[3] > alphaCutoff) ? 1.0f : 0.0f;
                        }
                        if (use_log_filter) {
                            for (int c = 0; c < log_channels; ++c) {
                                color[c] = std::log1p(color[c]);
                            }
                        }
                        return color;
                    };

                // Helper to sample normal map
                auto sample_normal = [&](int x, int y) -> std::array<float, 3>
                    {
                        x = wrapCoord(x, prev_w, u_wrap, is_cube);
                        y = wrapCoord(y, prev_h, v_wrap, is_cube);
                        std::array<float, 3> n{ 0.0f, 0.0f, 1.0f };

                        auto normal_span = pNormalBitmap->getMip<std::byte>(level - 1, face);
                        const std::byte* ptr = &normal_span[(y * prev_w + x) * normal_bytes_per_pixel];
                        float rx, gy;
                        if (is_normal_int) {
                            const xcolori& pixel = *reinterpret_cast<const xcolori*>(ptr);
                            rx = pixel.m_R / 255.0f;
                            gy = pixel.m_G / 255.0f;
                        }
                        else if (is_normal_half) {
                            const uint16_t* hptr = reinterpret_cast<const uint16_t*>(ptr);
                            xmath::half h;
                            h.m_Value = hptr[0];
                            rx = static_cast<float>(h);
                            h.m_Value = hptr[1];
                            gy = static_cast<float>(h);
                        }
                        else {
                            const float* fptr = reinterpret_cast<const float*>(ptr);
                            rx = fptr[0];
                            gy = fptr[1];
                        }
                        float normal_unpack_scale = is_normal_int ? 2.0f : 1.0f;
                        float normal_unpack_bias = is_normal_int ? -1.0f : 0.0f;
                        float ix = rx * normal_unpack_scale + normal_unpack_bias;
                        float iy = gy * normal_unpack_scale + normal_unpack_bias;
                        float iz = std::sqrt(std::max(0.0f, 1.0f - ix * ix - iy * iy)); // Assume positive z
                        n[0] = ix;
                        n[1] = iy;
                        n[2] = iz;
                        return n;
                    };

                // Horizontal pass
                for (int y = 0; y < prev_h; ++y)
                {
                    for (int x = 0; x < cur_w; ++x)
                    {
                        const int sx = x * 2;
                        float* hptr = &horiz[(y * cur_w + x) * channels];
                        if (is_box)
                        {
                            auto c0 = sample(sx, y);
                            auto c1 = sample(sx + 1, y);
                            for (int c = 0; c < channels; ++c)
                            {
                                hptr[c] = (c0[c] + c1[c]) * 0.5f;
                            }
                        }
                        else
                        {
                            std::array<float, 4> accum{ 0.0f, 0.0f, 0.0f, 0.0f };
                            for (int tap = -radius; tap <= radius; ++tap)
                            {
                                float weight = kernel_span[tap + radius];
                                auto scolor = sample(sx + tap, y);
                                for (int c = 0; c < channels; ++c)
                                {
                                    accum[c] += weight * scolor[c];
                                }
                            }
                            for (int c = 0; c < channels; ++c)
                            {
                                hptr[c] = accum[c];
                            }
                        }

                        if (use_toksvig) {
                            float* hnptr = &horiz_normal[(y * cur_w + x) * normal_channels];
                            if (is_box) {
                                auto n0 = sample_normal(sx, y);
                                auto n1 = sample_normal(sx + 1, y);
                                for (int c = 0; c < normal_channels; ++c) {
                                    hnptr[c] = (n0[c] + n1[c]) * 0.5f;
                                }
                            }
                            else {
                                std::array<float, 3> accum_n{ 0.0f, 0.0f, 0.0f };
                                for (int tap = -radius; tap <= radius; ++tap)
                                {
                                    float weight = kernel_span[tap + radius];
                                    auto sn = sample_normal(sx + tap, y);
                                    for (int c = 0; c < normal_channels; ++c) {
                                        accum_n[c] += weight * sn[c];
                                    }
                                }
                                for (int c = 0; c < normal_channels; ++c) {
                                    hnptr[c] = accum_n[c];
                                }
                            }
                        }
                    }
                }

                // Vertical pass + output
                for (int y = 0; y < cur_h; ++y)
                {
                    for (int x = 0; x < cur_w; ++x)
                    {
                        const int sy = y * 2;
                        std::array<float, 4> accum{ 0.0f, 0.0f, 0.0f, 0.0f };
                        if (is_box)
                        {
                            auto c0 = std::array<float, 4>{};
                            auto c1 = std::array<float, 4>{};
                            int wrapped_sy = wrapCoord(sy, prev_h, v_wrap, is_cube);
                            int wrapped_sy1 = wrapCoord(sy + 1, prev_h, v_wrap, is_cube);
                            const float* h0 = &horiz[(wrapped_sy * cur_w + x) * channels];
                            const float* h1 = &horiz[(wrapped_sy1 * cur_w + x) * channels];
                            for (int c = 0; c < channels; ++c)
                            {
                                c0[c] = h0[c];
                                c1[c] = h1[c];
                            }
                            for (int c = 0; c < channels; ++c)
                            {
                                accum[c] = (c0[c] + c1[c]) * 0.5f;
                            }
                        }
                        else
                        {
                            for (int tap = -radius; tap <= radius; ++tap)
                            {
                                float weight = kernel_span[tap + radius];
                                int wrapped_y = wrapCoord(sy + tap, prev_h, v_wrap, is_cube);
                                const float* hptr = &horiz[(wrapped_y * cur_w + x) * channels];
                                for (int c = 0; c < channels; ++c)
                                {
                                    accum[c] += weight * hptr[c];
                                }
                            }
                        }

                        if (use_log_filter) {
                            for (int c = 0; c < log_channels; ++c) {
                                accum[c] = std::max(accum[c], 0.0f);
                                accum[c] = std::expm1(accum[c]);
                            }
                        }

                        std::array<float, 3> accum_n{ 0.0f, 0.0f, 0.0f };
                        if (use_toksvig) {
                            if (is_box) {
                                auto n0 = std::array<float, 3>{};
                                auto n1 = std::array<float, 3>{};
                                int wrapped_sy = wrapCoord(sy, prev_h, v_wrap, is_cube);
                                int wrapped_sy1 = wrapCoord(sy + 1, prev_h, v_wrap, is_cube);
                                const float* hn0 = &horiz_normal[(wrapped_sy * cur_w + x) * normal_channels];
                                const float* hn1 = &horiz_normal[(wrapped_sy1 * cur_w + x) * normal_channels];
                                for (int c = 0; c < normal_channels; ++c) {
                                    n0[c] = hn0[c];
                                    n1[c] = hn1[c];
                                }
                                for (int c = 0; c < normal_channels; ++c) {
                                    accum_n[c] = (n0[c] + n1[c]) * 0.5f;
                                }
                            }
                            else {
                                for (int tap = -radius; tap <= radius; ++tap) {
                                    float weight = kernel_span[tap + radius];
                                    int wrapped_y = wrapCoord(sy + tap, prev_h, v_wrap, is_cube);
                                    const float* hnptr = &horiz_normal[(wrapped_y * cur_w + x) * normal_channels];
                                    for (int c = 0; c < normal_channels; ++c) {
                                        accum_n[c] += weight * hnptr[c];
                                    }
                                }
                            }

                            float len_sq = accum_n[0] * accum_n[0] + accum_n[1] * accum_n[1] + accum_n[2] * accum_n[2];
                            float len = std::sqrt(len_sq);
                            if (len > 0.0001f) {
                                float variance = 1.0f - len_sq; // Approx variance = 1 - Ft^2, Ft = len
                                accum[roughnessChannel] = std::sqrt(accum[roughnessChannel] * accum[roughnessChannel] + variance);
                            }
                        }

                        // Renormalize for normals (assume RGB normal, A unused or height)
                        if (renormalizeNormals && channels >= 3)
                        {
                            float main_unpack_scale = is_int ? 2.0f : 1.0f;
                            float main_unpack_bias = is_int ? -1.0f : 0.0f;
                            xmath::fvec3 n(accum[0] * main_unpack_scale + main_unpack_bias,
                                accum[1] * main_unpack_scale + main_unpack_bias,
                                accum[2] * main_unpack_scale + main_unpack_bias);
                            n = n.NormalizeSafeCopy();
                            float main_pack_scale = is_int ? 0.5f : 1.0f;
                            float main_pack_bias = is_int ? 0.5f : 0.0f;
                            accum[0] = n[0] * main_pack_scale + main_pack_bias;
                            accum[1] = n[1] * main_pack_scale + main_pack_bias;
                            accum[2] = n[2] * main_pack_scale + main_pack_bias;
                            // A unchanged
                        }

                        // Clamp for HDR to prevent negatives/overflow from ringing
                        if (!renormalizeNormals && (is_float32 || is_half))
                        {
                            for (int c = 0; c < channels; ++c)
                            {
                                accum[c] = std::max(0.0f, accum[c]);
                                if (is_half)
                                {
                                    accum[c] = std::min(accum[c], 65504.0f);
                                }
                            }
                        }

                        // Store
                        std::byte* optr = &cur_span[(y * cur_w + x) * bytes_per_pixel];
                        if (is_int)
                        {
                            float r = accum[0];
                            float g = accum[1];
                            float b = accum[2];
                            float a = (channels >= 4) ? accum[3] : 1.0f;
                            if (is_srgb && !renormalizeNormals)
                            {
                                r = std::pow(std::clamp(r, 0.0f, 1.0f), 1.0f / 2.2f);
                                g = std::pow(std::clamp(g, 0.0f, 1.0f), 1.0f / 2.2f);
                                b = std::pow(std::clamp(b, 0.0f, 1.0f), 1.0f / 2.2f);
                                // A linear
                            }
                            xcolori pixel;
                            pixel.m_R = static_cast<std::uint8_t>(std::clamp(r * 255.0f + 0.5f, 0.0f, 255.0f));
                            pixel.m_G = static_cast<std::uint8_t>(std::clamp(g * 255.0f + 0.5f, 0.0f, 255.0f));
                            pixel.m_B = static_cast<std::uint8_t>(std::clamp(b * 255.0f + 0.5f, 0.0f, 255.0f));
                            pixel.m_A = static_cast<std::uint8_t>(std::clamp(a * 255.0f + 0.5f, 0.0f, 255.0f));
                            std::memcpy(optr, &pixel, bytes_per_pixel);
                        }
                        else if (is_half)
                        {
                            uint16_t* hptr = reinterpret_cast<uint16_t*>(optr);
                            for (int c = 0; c < channels; ++c)
                            {
                                hptr[c] = xmath::half(accum[c]).m_Value;
                            }
                        }
                        else
                        { // float32
                            float* fptr = reinterpret_cast<float*>(optr);
                            for (int c = 0; c < channels; ++c)
                            {
                                fptr[c] = accum[c];
                            }
                        }
                    }
                }
            }
        }
    }
}
*/
