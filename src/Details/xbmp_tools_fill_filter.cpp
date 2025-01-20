namespace xbmp::tools::filters
{
    //--------------------------------------------------------------------------------

    void ForcePunchThroughAlpha(xcore::bitmap& Bitmap, int AlphaThreshold ) noexcept
    {
        assert( Bitmap.getFormat() == xcore::bitmap::format::R8G8B8A8 );
        assert(Bitmap.getMipCount() == 1);

        auto ColorView = Bitmap.getMip<xcore::icolor>(0);

        for (auto& E : ColorView )
        {
            if (E.m_A <= AlphaThreshold) E.m_A = 0x00;
            else                         E.m_A = 0xFF;
        }
    }

    //--------------------------------------------------------------------------------

    static
    xcore::icolor ComputetAvgSurroundingOpaqueColor(const xcore::bitmap& Bitmap, int TexelX, int TexelY, std::uint8_t AlphaThreshold ) noexcept
    {
        assert(Bitmap.getFormat() == xcore::bitmap::format::R8G8B8A8);
        assert(Bitmap.getMipCount() == 1);

        const auto      Color   = Bitmap.getMip<xcore::icolor>(0);
        const int       Width   = static_cast<int>(Bitmap.getWidth());
        const int       Height  = static_cast<int>(Bitmap.getHeight());

        std::uint32_t   SumR        = 0;
        std::uint32_t   SumG        = 0;
        std::uint32_t   SumB        = 0;
        std::uint32_t   SumTotal    = 0;

        xcore::icolor   FinalColor;

        bool            YWrapping   = false;
        bool            XWrapping   = false;
        bool            YMirror     = false;
        bool            XMirror     = false;

        if( Bitmap.getWrapMode() != xcore::bitmap::wrap_mode::UV_BOTH_CLAMP_TO_EDGE )
        {
            YWrapping = Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_BOTH_WRAP
                     || Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_UCLAMP_VWRAP;

            if (!YWrapping) YMirror = Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_BOTH_MIRROR
                                   || Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_UCLAMP_VMIRROR;

            XWrapping = Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_BOTH_WRAP
                     || Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_UWRAP_VCLAMP;

            if (!XWrapping) XMirror = Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_BOTH_MIRROR
                                   || Bitmap.getWrapMode() == xcore::bitmap::wrap_mode::UV_UMIRROR_VCLAMP;

        }

        constexpr auto ComputeWithWrapMode = []( int& v, const bool Wrapping, const bool Mirror, const int Edge ) constexpr ->bool
        {
            if (Wrapping)
            {
                if (v < 0) v = Edge - 1;
                else if (v >= Edge) v = 0;
            }
            else if (Mirror)
            {
                if (v < 0) v = 0;
                else if (v >= Edge) v = Edge - 1;
            }
            else
            {
                if ((v >= Edge) || (v < 0))
                    return true; // skip this line - out of range
            }

            return false;
        };

        for( int yi = TexelY - 1; yi <= (TexelY + 1); ++yi)
        {
            int Y = yi;

            if (ComputeWithWrapMode( Y, YWrapping, YMirror, Height))
                continue;
            
            for( int xi = TexelX - 1; xi <= (TexelX + 1); ++xi)
            {
                int X = xi;

                if (ComputeWithWrapMode(X, XWrapping, XMirror, Width))
                    continue;

                auto& C = Color[X + Y * Width];

                if( xi == TexelX && yi == TexelY )
                {
                    FinalColor = C;
                    continue;
                }

                // transparent pixel
                if (C.m_A < AlphaThreshold)
                    continue; 

                SumR += C.m_R;
                SumG += C.m_G;
                SumB += C.m_B;
                ++SumTotal;
            }
        }

        if (SumTotal > 0)
        {
            FinalColor.m_R = static_cast<std::uint8_t>(SumR / SumTotal);
            FinalColor.m_G = static_cast<std::uint8_t>(SumG / SumTotal);
            FinalColor.m_B = static_cast<std::uint8_t>(SumB / SumTotal);
        }

        return FinalColor;
    }

    //--------------------------------------------------------------------------------

    void FillAvrColorBaseOnAlpha(xcore::bitmap& Bitmap, const std::uint8_t AlphaThreshold, std::uint32_t Depth) noexcept
    {
        assert(Bitmap.getFormat() == xcore::bitmap::format::R8G8B8A8);
        assert(Bitmap.getMipCount() == 1);

        // Allocate our working buffers
        std::array<xcore::bitmap, 2> Bitmaps;
        for (auto& E : Bitmaps)
            E.CreateBitmap(Bitmap.getWidth(), Bitmap.getHeight());

        // Copy the data from the original buffer
        assert(Bitmap.getMip<std::byte>(0).size() == Bitmaps[0].getMip<std::byte>(0).size());
        std::memcpy( Bitmaps[0].getMip<std::byte>(0).data(), Bitmap.getMip<std::byte>(0).data(), Bitmap.getMip<std::byte>(0).size() );

        const auto  Width       = Bitmap.getWidth();
        const auto  Height      = Bitmap.getHeight();
        int         DepthIndex  = 0;
        int         DepthCount  = 0;

        while(Depth--)
        {
            int         PixelChanged    = 0;
            const auto  CS              = Bitmaps[    DepthIndex].getMip<xcore::icolor>(0);
            auto        CD              = Bitmaps[1 - DepthIndex].getMip<xcore::icolor>(0);

            for( auto y=0u ; y<Height; ++y )
            {
                for (auto x = 0u; x < Width; ++x)
                {
                    auto C = CS[ x + y * Width ];
                    if (C.m_A < AlphaThreshold)
                    {
                        C = ComputetAvgSurroundingOpaqueColor(Bitmaps[DepthIndex], x, y, AlphaThreshold);

                        // Check if any thing changed
                        if(reinterpret_cast<std::uint32_t&>(C) != reinterpret_cast<std::uint32_t&>(CS[x + y * Width]) ) 
                        {
                            PixelChanged++;

                            // mark this pixel as computed
                            C.m_A = 0xff;
                        }
                    }

                    // Write the final pixel
                    CD[x + y * Width] = C;
                }
            }

            // Switch buffers
            DepthIndex = 1 - DepthIndex;

            // If we have no changed any pixels we can call it done...
            if (PixelChanged == 0) break;
        }

        //
        // Copy the results to the original bitmap
        // Note that we don't overrite the alpha channel since we used it as a flag
        //
        {
            const auto CS = Bitmaps[DepthIndex].getMip<xcore::icolor>(0);
            auto       CD = Bitmap.getMip<xcore::icolor>(0);

            for ( auto i = 0; i< CS.size(); ++i )
            {
                CD[i].m_R = CS[i].m_R;
                CD[i].m_G = CS[i].m_G;
                CD[i].m_B = CS[i].m_B;
            }
        }
    }

    //--------------------------------------------------------------------------------

    void MakeBitmapTilable(xcore::bitmap& Bitmap, float WidthOverlapPercentage, float HeightOverlapPercentage ) noexcept
    {
        assert(WidthOverlapPercentage  >= 0);
        assert(WidthOverlapPercentage  <= 1);
        assert(HeightOverlapPercentage >= 0);
        assert(HeightOverlapPercentage <= 1);

        auto       Dest = Bitmap.getMip<xcore::icolor>(0);
        const auto W    = Bitmap.getWidth();
        const auto H    = Bitmap.getHeight();

        std::uint32_t MixSize = static_cast<std::uint32_t>(W * WidthOverlapPercentage);

        // Blend right and left edges
        for (auto y = 0u; y < H; ++y)
        {
            for (auto x = 0u; x < MixSize; ++x)
            {
                const float alpha = float(x) / MixSize;

                auto& L = Dest[y * W + x];
                auto& R = Dest[y * W + (W - x - 1)];

                L.m_R = static_cast<std::uint8_t>((1 - alpha) * R.m_R + alpha * L.m_R);
                L.m_G = static_cast<std::uint8_t>((1 - alpha) * R.m_G + alpha * L.m_G);
                L.m_B = static_cast<std::uint8_t>((1 - alpha) * R.m_B + alpha * L.m_B);
                L.m_A = static_cast<std::uint8_t>((1 - alpha) * R.m_A + alpha * L.m_A);

                R.m_R = static_cast<std::uint8_t>((1 - alpha) * L.m_R + alpha * R.m_R);
                R.m_G = static_cast<std::uint8_t>((1 - alpha) * L.m_G + alpha * R.m_G);
                R.m_B = static_cast<std::uint8_t>((1 - alpha) * L.m_B + alpha * R.m_B);
                R.m_A = static_cast<std::uint8_t>((1 - alpha) * L.m_A + alpha * R.m_A);
            }
        }

        // Blend top and bottom edges
        MixSize = static_cast<std::uint32_t>(H * HeightOverlapPercentage);

        for (auto x = 0u; x < W; ++x)
        {
            for (auto y = 0u; y < MixSize; ++y)
            {
                const float     alpha = float(y) / MixSize;

                auto& T = Dest[y * W + x];
                auto& B = Dest[x + (H - y - 1) * W];

                T.m_R = static_cast<std::uint8_t>((1 - alpha) * B.m_R + alpha * T.m_R);
                T.m_G = static_cast<std::uint8_t>((1 - alpha) * B.m_G + alpha * T.m_G);
                T.m_B = static_cast<std::uint8_t>((1 - alpha) * B.m_B + alpha * T.m_B);
                T.m_A = static_cast<std::uint8_t>((1 - alpha) * B.m_A + alpha * T.m_A);

                B.m_R = static_cast<std::uint8_t>((1 - alpha) * T.m_R + alpha * B.m_R);
                B.m_G = static_cast<std::uint8_t>((1 - alpha) * T.m_G + alpha * B.m_G);
                B.m_B = static_cast<std::uint8_t>((1 - alpha) * T.m_B + alpha * B.m_B);
                B.m_A = static_cast<std::uint8_t>((1 - alpha) * T.m_A + alpha * B.m_A);
            }
        }
    }
}
