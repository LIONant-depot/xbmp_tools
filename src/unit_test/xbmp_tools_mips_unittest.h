#pragma once
//
// Mip generation (xbmp::tools::mips::GenerateMipMaps): what each kind of texture must keep across its mips.
// Each check measures against the meaning of the source, not against a previous output.
//
#include <cmath>
#include <cstdio>
#include <vector>

namespace xbmp::tools::mips::unit_test
{
    inline int s_Failures = 0;

    inline void Check(bool bOk, const char* pWhat, double Got, double Expected)
    {
        std::printf("  %s %-64s got %9.4f expected %9.4f\n", bOk ? "ok  " : "FAIL", pWhat, Got, Expected);
        if (!bOk) ++s_Failures;
    }

    inline xbitmap MakeRGBA(int W, int H, xbitmap::color_space Space, auto&& Texel)
    {
        xbitmap B;
        B.CreateBitmap(W, H);
        B.setColorSpace(Space);
        B.setUWrapMode(xbitmap::wrap_mode::WRAP);
        B.setVWrapMode(xbitmap::wrap_mode::WRAP);
        auto Px = B.getMip<xcolori>(0);
        for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) Px[x + y * W] = Texel(x, y);
        return B;
    }

    //------------------------------------------------------------------------------
    // sRGB color is filtered in linear light with the sRGB curve: a 0 / 64 checkerboard averages to linear 0.02566, sRGB 44
    // (a pow 2.2 approximation gives 47; 0 / 255 would not tell them apart: 188 vs 187)

    inline void TestSRGBLinearFiltering()
    {
        std::printf("sRGB filtered in linear light\n");
        auto B = MakeRGBA(64, 64, xbitmap::color_space::SRGB, [](int x, int y) { const std::uint8_t v = ((x + y) & 1) ? 64 : 0; return xcolori{ v, v, v, 255 }; });
        GenerateMipMaps(B, 1, false, mips_filter_type::BOX);
        const auto Px = B.getMip<xcolori>(1);
        double Sum = 0; for (auto& P : Px) Sum += P.m_R;
        Check(std::abs(Sum / Px.size() - 44.0) <= 1.0, "mip 1 of a 0/64 sRGB checkerboard (box)", Sum / Px.size(), 44.0);
    }

    //------------------------------------------------------------------------------
    // A cutout keeps the fraction of texels that pass the alpha test (Castano 2010)

    inline double CoverageOf(std::span<const xcolori> Px, float Cutoff)
    {
        int Pass = 0; for (auto& P : Px) Pass += P.m_A > Cutoff * 255.0f;
        return double(Pass) / Px.size();
    }

    inline std::uint32_t Hash(std::uint32_t x) { x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16; return x; }

    inline void TestAlphaCoverage()
    {
        constexpr float Cutoff = 0.5f;

        // Foliage-like: soft-edged blobs (value noise, 16 texel cells) plus thin strands one texel wide in random columns
        std::printf("Cutout alpha keeps its coverage at the cutoff (blobs + thin strands)\n");
        {
            auto Noise = [](int x, int y)
            {
                const int cx = x / 16, cy = y / 16; const float fx = (x % 16) / 16.0f, fy = (y % 16) / 16.0f;
                auto V = [](int i, int j) { return (Hash(std::uint32_t(i & 15) * 31u + std::uint32_t(j & 15) * 977u) & 0xffff) / 65535.0f; };
                const float a = V(cx, cy) + (V(cx + 1, cy) - V(cx, cy)) * fx, b = V(cx, cy + 1) + (V(cx + 1, cy + 1) - V(cx, cy + 1)) * fx;
                return a + (b - a) * fy;
            };
            auto B = MakeRGBA(256, 256, xbitmap::color_space::SRGB, [&](int x, int y)
            {
                const float Blob   = std::clamp((Noise(x, y) - 0.55f) / 0.1f, 0.0f, 1.0f);    // soft edges over ~2 texels
                const auto  h      = Hash(std::uint32_t(x));                                       // a thin vertical strand in ~1 of 9 columns,
                const float Strand = (h % 9) == 0 ? 0.55f + 0.45f * ((h >> 8) & 255) / 255.0f : 0.0f;   // each its own opacity (as real ones)
                return xcolori{ 40, 160, 40, std::uint8_t(std::max(Blob, Strand) * 255.0f + 0.5f) };
            });
            const double Target = CoverageOf(B.getMip<xcolori>(0), Cutoff);
            GenerateMipMaps(B, 1, false, mips_filter_type::KAISER, true, Cutoff);
            for (int m = 1; m < B.getMipCount() && (256 >> m) >= 8; ++m)
            {
                const double Got = CoverageOf(B.getMip<xcolori>(m), Cutoff);
                char What[64]; std::snprintf(What, sizeof(What), "coverage of mip %d (%dx%d)", m, 256 >> m, 256 >> m);
                Check(std::abs(Got - Target) <= 0.05, What, Got, Target);
            }
        }

        // A regular fence: wires one texel wide every 4. From mip 2 every texel has the same alpha, so an alpha test can only pass all
        // or nothing (25% needs hashed alpha / alpha-to-coverage at run time); what it must never do is vanish (the old result: 0%)
        std::printf("A regular fence never vanishes\n");
        {
            auto B = MakeRGBA(256, 256, xbitmap::color_space::SRGB, [](int x, int) { return xcolori{ 90, 90, 90, std::uint8_t((x % 4) == 0 ? 255 : 0) }; });
            GenerateMipMaps(B, 1, false, mips_filter_type::KAISER, true, Cutoff);
            for (int m = 1; m < B.getMipCount() && (256 >> m) >= 8; ++m)
            {
                const double Got = CoverageOf(B.getMip<xcolori>(m), Cutoff);
                char What[64]; std::snprintf(What, sizeof(What), "fence coverage of mip %d is above 0", m);
                Check(Got > 0.0, What, Got, 0.25);
            }
        }
    }

    //------------------------------------------------------------------------------
    // A tangent-space normal keeps its direction: a constant normal (0.8, 0, 0.6) (R = x, G = y, B = z) stays x = 0.8 at every mip

    inline void TestNormalKeepsTilt()
    {
        std::printf("Normal map keeps its tilt\n");
        auto B = MakeRGBA(64, 64, xbitmap::color_space::LINEAR, [](int, int) { return xcolori{ std::uint8_t(0.8f * 127.5f + 127.5f), 128, std::uint8_t(0.6f * 127.5f + 127.5f), 255 }; });
        GenerateMipMaps(B, 1, true, mips_filter_type::KAISER);
        const int Last = B.getMipCount() - 1;
        const auto P = B.getMip<xcolori>(Last)[0];
        Check(std::abs((P.m_R / 127.5 - 1.0) - 0.8) <= 0.02, "x of the last mip of a constant normal (0.8, 0, 0.6)", P.m_R / 127.5 - 1.0, 0.8);
    }

    inline int Test()
    {
        TestSRGBLinearFiltering();
        TestAlphaCoverage();
        TestNormalKeepsTilt();
        std::printf("%s (%d failed)\n", s_Failures ? "FAILED" : "ALL OK", s_Failures);
        return s_Failures;
    }
}
