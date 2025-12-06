#include <cmath>
#include <algorithm>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>


namespace xbmp::tools::lighting
{
    static constexpr double PI      = 3.141592653589793238462643383279502884;
    static constexpr float  F_PI    = 3.14159265358979323846f;

    inline
    float pow5(float x)
    {
        const float x2 = x * x;
        return x2 * x2 * x;
    }

    inline
    float saturate(float x)
    {
        return std::min(std::max(x, 0.0f), 1.0f);
    }

    xmath::fvec2 hammersley(uint32_t i, float invN)
    {
        uint32_t bits = i;
        bits = (bits << 16u) | (bits >> 16u);
        bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
        bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
        bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
        bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
        float rdi = static_cast<float>(bits) * 2.3283064365386963e-10f; // / 0x100000000
        return {static_cast<float>(i) * invN, rdi};
    }

    xmath::fvec3 hemisphereImportanceSampleDggx(const xmath::fvec2& u, float a)
    {
        const float phi         = 2.0f * F_PI * u.m_X;
        const float cosTheta2   = (1 - u.m_Y) / (1 + (a + 1) * ((a - 1) * u.m_Y));
        const float cosTheta    = std::sqrt(cosTheta2);
        const float sinTheta    = std::sqrt(1 - cosTheta2);
        return {sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
    }

    float Visibility(float NoV, float NoL, float linearRoughness)
    {
        const float a       = linearRoughness;
        const float a2      = a * a;
        const float GGXV    = NoL * std::sqrt(NoV * NoV * (1.0f - a2) + a2);
        const float GGXL    = NoV * std::sqrt(NoL * NoL * (1.0f - a2) + a2);
        return 0.5f / (GGXV + GGXL);
    }

    xmath::fvec2 integrateDFG(float NoV, float linearRoughness, size_t numSamples)
    {
        float sumX = 0.0f;
        float sumY = 0.0f;
        const xmath::fvec3 V = { std::sqrt(1.0f - NoV * NoV), 0.0f, NoV };
        const float invN = 1.0f / static_cast<float>(numSamples);

#pragma omp parallel for reduction(+:sumX,sumY)
        for (int i = 0; i < static_cast<int>(numSamples); ++i)
        {
            const xmath::fvec2 u = hammersley(static_cast<uint32_t>(i), invN);
            const xmath::fvec3 H = hemisphereImportanceSampleDggx(u, linearRoughness);
            const xmath::fvec3 L = 2.0f * xmath::fvec3::Dot(V, H) * H - V;
            const float VoH = saturate(xmath::fvec3::Dot(V, H));
            const float NoL = saturate(L.m_Z);
            const float NoH = saturate(H.m_Z);

            if (NoL > 0.0f)
            {
                const float v = Visibility(NoV, NoL, linearRoughness) * NoL * (VoH / NoH);
                const float Fc = pow5(1.0f - VoH);
                sumX += v * (1.0f - Fc);
                sumY += v * Fc;
            }
        }

        return { sumX * (4.0f * invN), sumY * (4.0f * invN) };
    }

    void GenerateGGX_BRDF_RG_32FLOAT(xbitmap& Bitmap, int Width, int Height) noexcept
    {
        const int width = Width;
        const int height = Height;
        const size_t numSamples = 1024;
        const size_t data_size = static_cast<size_t>(width) * height * 2;
        const size_t data_byte_size = (data_size + 1) * sizeof(float);

        auto data_owner = std::make_unique<float[]>(data_byte_size);
        auto data = std::span{ data_owner.get() + 1, data_size };

        data_owner[0] = 0.0f;

#pragma omp parallel for
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                const float NoV = saturate((static_cast<float>(x) + 0.5f) / static_cast<float>(width));
                const float roughness = saturate((static_cast<float>(height - y) + 0.5f) / static_cast<float>(height));
                const float linearRoughness = roughness * roughness;

                const xmath::fvec2 dfg = integrateDFG(NoV, linearRoughness, numSamples);

                const size_t index = (static_cast<size_t>(y) * width + x) * 2;
                data[index]     = dfg.m_X;
                data[index + 1] = dfg.m_Y;
            }
        }

        Bitmap.setup
        (Width
            , Height
            , xbitmap::format::R32G32_FLOAT
            , data_size * sizeof(float)
            , std::span<std::byte>{reinterpret_cast<std::byte*>(data_owner.release()), data_byte_size}
        , true
            , 1
            , 1
            , false
            );
    }


    void GenerateGGX_BRDF_RG_16SFLOAT(xbitmap& Bitmap, int Width, int Height) noexcept
    {
        const size_t width          = Width;
        const size_t height         = Height;
        const size_t numSamples     = 1024;
        const size_t data_size      = width * height * 2;   // (2 = r,g)
        const size_t data_byte_size = (data_size + 2) * sizeof(xmath::half);        // +2 because xbitmap requires the first 4 bytes as an offset

        auto data_owner = std::make_unique<xmath::half[]>(data_byte_size);
        auto data       = std::span{ data_owner.get()+2, data_size };

        // set first 4 bytes to zero
        data_owner[0].m_Value = 0;
        data_owner[1].m_Value = 0;

#pragma omp parallel for
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                const float         NoV             = saturate((static_cast<float>(x) + 0.5f) / static_cast<float>(width));
                const float         roughness       = saturate((static_cast<float>(height - y) + 0.5f) / static_cast<float>(height));
                const float         linearRoughness = roughness * roughness;

                const xmath::fvec2  dfg             = integrateDFG(NoV, linearRoughness, numSamples);

                const size_t index = (y * width + x) * 2;
                data[index]     = dfg.m_X;  // scale / red
                data[index + 1] = dfg.m_Y;  // bias / green
            }
        }

        Bitmap.setup
        ( Width
        , Height
        , xbitmap::format::R16G16_SFLOAT
        , data_size * sizeof(xmath::half)
        , std::span<std::byte>{ (std::byte*)data_owner.release(), data_byte_size }
        , true
        , 1
        , 1
        , false
        );
    }

    // New function for cubemap prefiltering
    xmath::fvec3 importanceSampleGGX(const xmath::fvec2& Xi, float roughness, const xmath::fvec3& N)
    {
        float a = roughness * roughness;

        float phi = 2.0f * F_PI * Xi[0];
        float cosTheta = std::sqrt((1.0f - Xi[1]) / (1.0f + (a * a - 1.0f) * Xi[1]));
        float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);

        xmath::fvec3 H(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);

        xmath::fvec3 up = (std::abs(N[2]) < 0.999f) ? xmath::fvec3(0.0f, 0.0f, 1.0f) : xmath::fvec3(1.0f, 0.0f, 0.0f);
        xmath::fvec3 tangent = up.Cross(N).NormalizeSafeCopy();
        xmath::fvec3 bitangent = N.Cross(tangent);

        return tangent * H[0] + bitangent * H[1] + N * H[2];
    }

    xmath::fvec3 importanceSampleCosine(const xmath::fvec2& Xi, const xmath::fvec3& N)
    {
        float phi = 2.0f * F_PI * Xi[0];
        float cosTheta = std::sqrt(1.0f - Xi[1]);
        float sinTheta = std::sqrt(Xi[1]);

        xmath::fvec3 H(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);

        xmath::fvec3 up = (std::abs(N[2]) < 0.999f) ? xmath::fvec3(0.0f, 0.0f, 1.0f) : xmath::fvec3(1.0f, 0.0f, 0.0f);
        xmath::fvec3 tangent = up.Cross(N).NormalizeSafeCopy();
        xmath::fvec3 bitangent = N.Cross(tangent);

        return tangent * H[0] + bitangent * H[1] + N * H[2];
    }

    struct CubemapFaceInfo
    {
        int face;
        xmath::fvec2 uv;
    };

    CubemapFaceInfo getCubemapFace(const xmath::fvec3& dir)
    {
        xmath::fvec3 ad = xmath::fvec3(std::abs(dir[0]), std::abs(dir[1]), std::abs(dir[2]));
        float maxAxis = std::max(ad[0], std::max(ad[1], ad[2]));

        CubemapFaceInfo info;

        if (maxAxis == ad[0])
        {
            if (dir[0] > 0.0f)
            {
                info.face = 0; // +X
                info.uv = xmath::fvec2(dir[2], -dir[1]) / ad[0] * 0.5f + 0.5f;
            }
            else
            {
                info.face = 1; // -X
                info.uv = xmath::fvec2(-dir[2], -dir[1]) / ad[0] * 0.5f + 0.5f;
            }
        }
        else if (maxAxis == ad[1])
        {
            if (dir[1] > 0.0f)
            {
                info.face = 2; // +Y
                info.uv = xmath::fvec2(dir[0], dir[2]) / ad[1] * 0.5f + 0.5f;
            }
            else
            {
                info.face = 3; // -Y
                info.uv = xmath::fvec2(dir[0], -dir[2]) / ad[1] * 0.5f + 0.5f;
            }
        }
        else
        {
            if (dir[2] > 0.0f)
            {
                info.face = 4; // +Z
                info.uv = xmath::fvec2(-dir[0], -dir[1]) / ad[2] * 0.5f + 0.5f;
            }
            else
            {
                info.face = 5; // -Z
                info.uv = xmath::fvec2(dir[0], -dir[1]) / ad[2] * 0.5f + 0.5f;
            }
        }

        return info;
    }

    xmath::fvec3 sampleCubemapBilinear(const xbitmap& cubemap, const xmath::fvec3& dir, int mip = 0)
    {
        CubemapFaceInfo info = getCubemapFace(dir.NormalizeSafeCopy());

        int size = cubemap.getWidth() >> mip;
        xmath::fvec2 pixel = info.uv * static_cast<float>(size) - 0.5f;
        xmath::fvec2 frac = xmath::fvec2(pixel[0] - std::floor(pixel[0]), pixel[1] - std::floor(pixel[1]));
        int x0 = static_cast<int>(std::floor(pixel[0]));
        int y0 = static_cast<int>(std::floor(pixel[1]));
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        auto getColor = [&](int x, int y) -> xmath::fvec3
            {
                x = std::clamp(x, 0, size - 1);
                y = std::clamp(y, 0, size - 1);

                const xbitmap::format   fmt = cubemap.getFormat();
                const bool              is_float32 = (fmt == xbitmap::format::R32G32B32_FLOAT     || fmt == xbitmap::format::R32G32B32A32_FLOAT);
                const bool              is_half    = (fmt == xbitmap::format::R16G16B16A16_SFLOAT);
                const bool              is_int     = (!is_float32 && !is_half); // Assume R8G8B8 or R8G8B8A8
                xmath::fvec3            color(0.0f);

                if (is_int)
                {
                    assert(fmt == xbitmap::format::R8G8B8A8 || fmt == xbitmap::format::R8G8B8U8 );
                    auto            span    = cubemap.getMip<xcolori>(mip, info.face);
                    const xcolori&  pixel   = span[y * size + x];
                    color[0] = pixel.m_R / 255.0f;
                    color[1] = pixel.m_G / 255.0f;
                    color[2] = pixel.m_B / 255.0f;
                }
                else if (is_half)
                {
                    assert(fmt == xbitmap::format::R16G16B16A16_SFLOAT);

                    auto            span = cubemap.getMip<uint16_t>(mip, info.face);
                    const uint16_t* hptr = &span[(y * size + x) * 4];
                    xmath::half h;
                    h.m_Value = hptr[0]; color[0] = static_cast<float>(h);
                    h.m_Value = hptr[1]; color[1] = static_cast<float>(h);
                    h.m_Value = hptr[2]; color[2] = static_cast<float>(h);
                }
                else // float32
                {
                    assert(fmt == xbitmap::format::R32G32B32A32_FLOAT || fmt == xbitmap::format::R32G32B32_FLOAT);

                    auto         span       = cubemap.getMip<float>(mip, info.face);
                    const int    dimensions = fmt == xbitmap::format::R32G32B32A32_FLOAT ? 4 : 3;
                    const float* fptr       = &span[(y * size + x) * dimensions];
                    color[0] = fptr[0];
                    color[1] = fptr[1];
                    color[2] = fptr[2];
                }
                return color;
            };

        xmath::fvec3 c00 = getColor(x0, y0);
        xmath::fvec3 c10 = getColor(x1, y0);
        xmath::fvec3 c01 = getColor(x0, y1);
        xmath::fvec3 c11 = getColor(x1, y1);

        xmath::fvec3 c0 = c00 * (1 - frac[0]) + c10 * frac[0];
        xmath::fvec3 c1 = c01 * (1 - frac[0]) + c11 * frac[0];

        return c0 * (1 - frac[1]) + c1 * frac[1];
    }

    xmath::fvec3 getCubemapDirection(int face, float u, float v)
    {
        xmath::fvec3 dir;
        switch (face)
        {
        case 0: // +X
            dir = xmath::fvec3(1.0f, -v, -u);
            break;
        case 1: // -X
            dir = xmath::fvec3(-1.0f, -v, u);
            break;
        case 2: // +Y
            dir = xmath::fvec3(u, 1.0f, v);
            break;
        case 3: // -Y
            dir = xmath::fvec3(u, -1.0f, -v);
            break;
        case 4: // +Z
            dir = xmath::fvec3(u, -v, 1.0f);
            break;
        case 5: // -Z
            dir = xmath::fvec3(-u, -v, -1.0f);
            break;
        }
        return dir.NormalizeSafeCopy();
    }

    void PrefilterCubemap
    ( const xbitmap&                inputCubemap
    , xbitmap&                      specularOutput
    , const int                     nMips
    , const int                     specularSamples
    , const diffuse_type            diffuse_type
    , xbitmap*                      diffuseCubemap
    , std::vector<xmath::fvec3>*    shCoeffs
    , const int                     diffuseRes
    , const int                     diffuseSamples
    )
    {
        if (!inputCubemap.isCubemap() || inputCubemap.getFaceCount() != 6 || inputCubemap.getWidth() != inputCubemap.getHeight()) return;

        int baseSize = inputCubemap.getWidth();
        xbitmap::format outputFmt = xbitmap::format::R16G16B16A16_SFLOAT; // Use 4 channels, A=1
        int channels = 4;
        int bytes_per_pixel = 8; // half *4

        const int numMips = nMips < 0 ? xbmp::tools::mips::ComputeMaxMips(1, baseSize, baseSize) : nMips;

        // For specular
        std::vector<std::uint64_t> mip_sizes(numMips);
        int w = baseSize;
        for (int level = 0; level < numMips; ++level)
        {
            mip_sizes[level] = static_cast<std::uint64_t>(w) * w * bytes_per_pixel;
            w = std::max(1, w / 2);
        }

        std::uint64_t face_size = 0;
        for (auto s : mip_sizes) face_size += s;

        auto MipTableSize = sizeof(xbitmap::mip) * numMips;

        std::uint64_t total_size = MipTableSize + face_size * 6;

        auto new_data = std::make_unique<std::byte[]>(total_size);

        auto* p_mip_table = reinterpret_cast<xbitmap::mip*>(new_data.get());
        int cur_offset = 0;
        for (int level = 0; level < numMips; ++level)
        {
            p_mip_table[level].m_Offset = cur_offset;
            cur_offset += static_cast<int>(mip_sizes[level]);
        }

        auto* p_data = new_data.get() + MipTableSize;

        specularOutput.setup(baseSize, baseSize, outputFmt, face_size, { new_data.release(), total_size }, true, numMips, 1, true);

        // Copy mip 0 from input, convert to half
        const xbitmap::format input_fmt = inputCubemap.getFormat();
        const bool input_float32 = (input_fmt == xbitmap::format::R32G32B32_FLOAT || input_fmt == xbitmap::format::R32G32B32A32_FLOAT);
        const bool input_half = (input_fmt == xbitmap::format::R16G16B16A16_SFLOAT);
        const bool input_int = (!input_float32 && !input_half);
        const int input_channels = (input_fmt == xbitmap::format::R32G32B32A32_FLOAT || input_fmt == xbitmap::format::R16G16B16A16_SFLOAT || input_fmt == xbitmap::format::R8G8B8A8) ? 4 : 3;

        for (int face = 0; face < 6; ++face)
        {
            auto dst_span_typed = specularOutput.getMip<uint16_t>(0, face);
            xmath::fvec3 color(0.0f);

            if (input_int)
            {
                assert(input_fmt == xbitmap::format::R8G8B8A8 || input_fmt == xbitmap::format::R8G8B8U8 );
                auto src_span_typed = inputCubemap.getMip<xcolori>(0, face);
                for (int y = 0; y < baseSize; ++y)
                {
                    for (int x = 0; x < baseSize; ++x)
                    {
                        const xcolori& pixel = src_span_typed[y * baseSize + x];
                        color[0] = pixel.m_R / 255.0f;
                        color[1] = pixel.m_G / 255.0f;
                        color[2] = pixel.m_B / 255.0f;

                        uint16_t* dst_ptr = &dst_span_typed[(y * baseSize + x) * 4];
                        dst_ptr[0] = xmath::half(color[0]).m_Value;
                        dst_ptr[1] = xmath::half(color[1]).m_Value;
                        dst_ptr[2] = xmath::half(color[2]).m_Value;
                        dst_ptr[3] = xmath::half(1.0f).m_Value;
                    }
                }
            }
            else if (input_half)
            {
                assert(input_fmt == xbitmap::format::R16G16B16A16_SFLOAT);

                auto src_span_typed = inputCubemap.getMip<uint16_t>(0, face);
                for (int y = 0; y < baseSize; ++y)
                {
                    for (int x = 0; x < baseSize; ++x)
                    {
                        const uint16_t* hptr = &src_span_typed[(y * baseSize + x) * input_channels];
                        xmath::half h;
                        h.m_Value = hptr[0]; color[0] = static_cast<float>(h);
                        h.m_Value = hptr[1]; color[1] = static_cast<float>(h);
                        h.m_Value = hptr[2]; color[2] = static_cast<float>(h);

                        uint16_t* dst_ptr = &dst_span_typed[(y * baseSize + x) * 4];
                        dst_ptr[0] = xmath::half(color[0]).m_Value;
                        dst_ptr[1] = xmath::half(color[1]).m_Value;
                        dst_ptr[2] = xmath::half(color[2]).m_Value;
                        dst_ptr[3] = xmath::half(1.0f).m_Value;
                    }
                }
            }
            else // float32
            {
                assert(input_fmt == xbitmap::format::R32G32B32A32_FLOAT || input_fmt == xbitmap::format::R32G32B32_FLOAT);

                auto src_span_typed = inputCubemap.getMip<float>(0, face);
                for (int y = 0; y < baseSize; ++y)
                {
                    for (int x = 0; x < baseSize; ++x)
                    {
                        const float* fptr = &src_span_typed[(y * baseSize + x) * input_channels];
                        color[0] = fptr[0];
                        color[1] = fptr[1];
                        color[2] = fptr[2];

                        uint16_t* dst_ptr = &dst_span_typed[(y * baseSize + x) * 4];
                        dst_ptr[0] = xmath::half(color[0]).m_Value;
                        dst_ptr[1] = xmath::half(color[1]).m_Value;
                        dst_ptr[2] = xmath::half(color[2]).m_Value;
                        dst_ptr[3] = xmath::half(1.0f).m_Value;
                    }
                }
            }
        }

        // Generate specular mips
        for (int level = 1; level < numMips; ++level)
        {
            float roughness = static_cast<float>(level) / static_cast<float>(numMips - 1);
            int cur_size = baseSize >> level;

#pragma omp parallel for
            for (int face = 0; face < 6; ++face)
            {
                auto dst_span = specularOutput.getMip<std::byte>(level, face);

#pragma omp parallel for
                for (int y = 0; y < cur_size; ++y)
                {
                    for (int x = 0; x < cur_size; ++x)
                    {
                        float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(cur_size) * 2.0f - 1.0f;
                        float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(cur_size) * 2.0f - 1.0f;
                        xmath::fvec3 N = getCubemapDirection(face, u, v);

                        xmath::fvec3 color(0.0f);
                        float totalWeight = 0.0f;

                        for (int s = 0; s < specularSamples; ++s)
                        {
                            xmath::fvec2 Xi = hammersley(s, 1/static_cast<float>(specularSamples));
                            xmath::fvec3 H = importanceSampleGGX(Xi, roughness, N);
                            xmath::fvec3 L = (2.0f * N.Dot(H) * H - N).NormalizeSafeCopy();
                            float NoL = std::max(0.0f, N.Dot(L));

                            if (NoL > 0.0f)
                            {
                                xmath::fvec3 sampleColor = sampleCubemapBilinear(inputCubemap, L);
                                color += sampleColor * NoL;
                                totalWeight += NoL;
                            }
                        }

                        if (totalWeight > 0.0f)
                        {
                            color /= totalWeight;
                        }

                        uint16_t* dst_ptr = reinterpret_cast<uint16_t*>(&dst_span[(y * cur_size + x) * bytes_per_pixel]);
                        dst_ptr[0] = xmath::half(color[0]).m_Value;
                        dst_ptr[1] = xmath::half(color[1]).m_Value;
                        dst_ptr[2] = xmath::half(color[2]).m_Value;
                        dst_ptr[3] = xmath::half(1.0f).m_Value;
                    }
                }
            }
        }

        // Diffuse
        if (diffuse_type == diffuse_type::CUBEMAP && diffuseCubemap != nullptr)
        {
            // Setup diffuse cubemap, no mips, small res
            std::uint64_t diffuse_face_size = static_cast<std::uint64_t>(diffuseRes) * diffuseRes * bytes_per_pixel;
            std::uint64_t diffuse_total_size = sizeof(xbitmap::mip) + diffuse_face_size * 6;

            auto diffuse_data = std::make_unique<std::byte[]>(diffuse_total_size);
            auto* diffuse_mip = reinterpret_cast<xbitmap::mip*>(diffuse_data.get());
            diffuse_mip[0].m_Offset = 0;

            auto* diffuse_p_data = diffuse_data.get() + sizeof(xbitmap::mip);

            diffuseCubemap->setup(diffuseRes, diffuseRes, outputFmt, diffuse_face_size, { diffuse_data.release(), diffuse_total_size }, true, 1, 1, true);

#pragma omp parallel for
            for (int face = 0; face < 6; ++face)
            {
                auto dst_span = diffuseCubemap->getMip<std::byte>(0, face);

#pragma omp parallel for
                for (int y = 0; y < diffuseRes; ++y)
                {
                    for (int x = 0; x < diffuseRes; ++x)
                    {
                        float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(diffuseRes) * 2.0f - 1.0f;
                        float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(diffuseRes) * 2.0f - 1.0f;
                        xmath::fvec3 N = getCubemapDirection(face, u, v);

                        xmath::fvec3 color(0.0f);

                        for (int s = 0; s < diffuseSamples; ++s)
                        {
                            xmath::fvec2 Xi = hammersley(s, 1/static_cast<float>(diffuseSamples));
                            xmath::fvec3 L = importanceSampleCosine(Xi, N);
                            color += sampleCubemapBilinear(inputCubemap, L);
                        }

                        color = color * static_cast<float>(PI) / static_cast<float>(diffuseSamples);

                        uint16_t* dst_ptr = reinterpret_cast<uint16_t*>(&dst_span[(y * diffuseRes + x) * bytes_per_pixel]);
                        dst_ptr[0] = xmath::half(color[0]).m_Value;
                        dst_ptr[1] = xmath::half(color[1]).m_Value;
                        dst_ptr[2] = xmath::half(color[2]).m_Value;
                        dst_ptr[3] = xmath::half(1.0f).m_Value;
                    }
                }
            }
        }
        else if (diffuse_type == diffuse_type::SH && shCoeffs != nullptr)
        {
            shCoeffs->resize(9, xmath::fvec3(0.0f));

            // Simple projection, loop over input cubemap pixels
            int res = inputCubemap.getWidth();
            float solidAngleAvg = 4.0f * F_PI / (6.0f * static_cast<float>(res * res));

            for (int face = 0; face < 6; ++face)
            {
                xmath::fvec3 color(0.0f);

                if (input_int)
                {
                    assert(input_fmt == xbitmap::format::R8G8B8A8 || input_fmt == xbitmap::format::R8G8B8U8 );
                    auto src_span_typed = inputCubemap.getMip<xcolori>(0, face);
                    for (int y = 0; y < res; ++y)
                    {
                        for (int x = 0; x < res; ++x)
                        {
                            float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(res) * 2.0f - 1.0f;
                            float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(res) * 2.0f - 1.0f;
                            xmath::fvec3 dir = getCubemapDirection(face, u, v);

                            const xcolori& pixel = src_span_typed[y * res + x];
                            color[0] = pixel.m_R / 255.0f;
                            color[1] = pixel.m_G / 255.0f;
                            color[2] = pixel.m_B / 255.0f;

                            // SH basis (order 2, 9 coeffs)
                            float sqrtPI = static_cast<float>(std::sqrt(PI));

                            float y00 = 0.282095f; // 1/2 sqrt(1/PI)

                            float y1n1 = 0.488603f * dir[1]; // sqrt(3/(4PI)) y

                            float y10 = 0.488603f * dir[2]; // z

                            float y1p1 = 0.488603f * dir[0]; // x

                            float y2n2 = 1.092548f * dir[0] * dir[1]; // sqrt(15/(4PI)) x y

                            float y2n1 = 1.092548f * dir[1] * dir[2]; // y z

                            float y20 = 0.315392f * (3.0f * dir[2] * dir[2] - 1.0f); // 1/4 sqrt(5/PI) (3z^2 -1)

                            float y2p1 = 1.092548f * dir[0] * dir[2]; // x z

                            float y2p2 = 0.546274f * (dir[0] * dir[0] - dir[1] * dir[1]); // 1/4 sqrt(15/PI) (x^2 - y^2)

                            std::array<float, 9> basis = { y00, y1n1, y10, y1p1, y2n2, y2n1, y20, y2p1, y2p2 };

                            for (int k = 0; k < 9; ++k)
                            {
                                (*shCoeffs)[k] += color * basis[k] * solidAngleAvg;
                            }
                        }
                    }
                }
                else if (input_half)
                {
                    assert(input_fmt == xbitmap::format::R16G16B16A16_SFLOAT );

                    auto src_span_typed = inputCubemap.getMip<uint16_t>(0, face);
                    for (int y = 0; y < res; ++y)
                    {
                        for (int x = 0; x < res; ++x)
                        {
                            float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(res) * 2.0f - 1.0f;
                            float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(res) * 2.0f - 1.0f;
                            xmath::fvec3 dir = getCubemapDirection(face, u, v);

                            const uint16_t* hptr = &src_span_typed[(y * res + x) * input_channels];
                            xmath::half h;
                            h.m_Value = hptr[0]; color[0] = static_cast<float>(h);
                            h.m_Value = hptr[1]; color[1] = static_cast<float>(h);
                            h.m_Value = hptr[2]; color[2] = static_cast<float>(h);

                            // SH basis (order 2, 9 coeffs)
                            float sqrtPI = static_cast<float>(std::sqrt(PI));

                            float y00 = 0.282095f; // 1/2 sqrt(1/PI)

                            float y1n1 = 0.488603f * dir[1]; // sqrt(3/(4PI)) y

                            float y10 = 0.488603f * dir[2]; // z

                            float y1p1 = 0.488603f * dir[0]; // x

                            float y2n2 = 1.092548f * dir[0] * dir[1]; // sqrt(15/(4PI)) x y

                            float y2n1 = 1.092548f * dir[1] * dir[2]; // y z

                            float y20 = 0.315392f * (3.0f * dir[2] * dir[2] - 1.0f); // 1/4 sqrt(5/PI) (3z^2 -1)

                            float y2p1 = 1.092548f * dir[0] * dir[2]; // x z

                            float y2p2 = 0.546274f * (dir[0] * dir[0] - dir[1] * dir[1]); // 1/4 sqrt(15/PI) (x^2 - y^2)

                            std::array<float, 9> basis = { y00, y1n1, y10, y1p1, y2n2, y2n1, y20, y2p1, y2p2 };

                            for (int k = 0; k < 9; ++k)
                            {
                                (*shCoeffs)[k] += color * basis[k] * solidAngleAvg;
                            }
                        }
                    }
                }
                else // float32
                {
                    assert(input_fmt == xbitmap::format::R32G32B32A32_FLOAT || input_fmt == xbitmap::format::R32G32B32_FLOAT);

                    auto src_span_typed = inputCubemap.getMip<float>(0, face);
                    for (int y = 0; y < res; ++y)
                    {
                        for (int x = 0; x < res; ++x)
                        {
                            float u = (static_cast<float>(x) + 0.5f) / static_cast<float>(res) * 2.0f - 1.0f;
                            float v = (static_cast<float>(y) + 0.5f) / static_cast<float>(res) * 2.0f - 1.0f;
                            xmath::fvec3 dir = getCubemapDirection(face, u, v);

                            const float* fptr = &src_span_typed[(y * res + x) * input_channels];
                            color[0] = fptr[0];
                            color[1] = fptr[1];
                            color[2] = fptr[2];

                            // SH basis (order 2, 9 coeffs)
                            float sqrtPI = static_cast<float>(std::sqrt(PI));

                            float y00 = 0.282095f; // 1/2 sqrt(1/PI)

                            float y1n1 = 0.488603f * dir[1]; // sqrt(3/(4PI)) y

                            float y10 = 0.488603f * dir[2]; // z

                            float y1p1 = 0.488603f * dir[0]; // x

                            float y2n2 = 1.092548f * dir[0] * dir[1]; // sqrt(15/(4PI)) x y

                            float y2n1 = 1.092548f * dir[1] * dir[2]; // y z

                            float y20 = 0.315392f * (3.0f * dir[2] * dir[2] - 1.0f); // 1/4 sqrt(5/PI) (3z^2 -1)

                            float y2p1 = 1.092548f * dir[0] * dir[2]; // x z

                            float y2p2 = 0.546274f * (dir[0] * dir[0] - dir[1] * dir[1]); // 1/4 sqrt(15/PI) (x^2 - y^2)

                            std::array<float, 9> basis = { y00, y1n1, y10, y1p1, y2n2, y2n1, y20, y2p1, y2p2 };

                            for (int k = 0; k < 9; ++k)
                            {
                                (*shCoeffs)[k] += color * basis[k] * solidAngleAvg;
                            }
                        }
                    }
                }
            }
        }
    }
}