//
// TIFF and DNG through tinydng (compiled in xbmp_tools_tinydng.c).
// 8, 16 and 32 bit integer and 32 bit float samples, 1 to 4 channels.
// A DNG's raw sensor image (CFA, not debayered) is skipped: its largest developed image (the preview) is loaded.
//
#include <algorithm>
#include <format>
#include <functional>
#include "dependencies/tinydng/tinydng.h"
#include "xbmp_tools_file.h"

namespace xbmp::tools::loader
{
    namespace details_tiff
    {
        // Decodes the largest non CFA image of the file to RGBA floats. Integers are mapped to 0..1, floats are kept.
        static xerr Decode( std::vector<float>& RGBA, int& W, int& H, std::wstring_view FileName ) noexcept
        {
            std::vector<std::uint8_t> FileData;
            {
                auto File = details::OpenFile(FileName, L"rb");
                if (File == nullptr || details::ReadFile(FileData, *File) == false)
                    return xerr::create_f<xbmp::tools::state, "Fail to read the TIFF/DNG file">();
            }

            tinydng_error Err{};
            std::unique_ptr<tinydng_context, decltype(&tinydng_context_destroy)> Ctx{ tinydng_context_create(nullptr, &Err), &tinydng_context_destroy };
            if (Ctx == nullptr) return xerr::create_f<xbmp::tools::state, "Fail to create the tinydng context">();

            tinydng_open_options Opts{};
            Opts.flags = TINYDNG_OPEN_PARSE_SUBIFDS;

            tinydng_document* pDoc = nullptr;
            if (tinydng_open_memory(Ctx.get(), FileData.data(), FileData.size(), &Opts, &pDoc, &Err) != TINYDNG_OK)
            {
                xerr::LogMessage<state::FAILURE>(std::format("[Error]: Trying to load TIFF/DNG: {} but found error {}", xstrtool::UTF16ToUTF8(FileName), Err.message));
                return xerr::create_f<xbmp::tools::state, "Fail to open the TIFF/DNG file">();
            }
            auto DocClose = std::unique_ptr<tinydng_document, std::function<void(tinydng_document*)>>{ pDoc, [&](tinydng_document* p) { tinydng_document_destroy(Ctx.get(), p); } };

            // Pick the largest image that is not a raw sensor (CFA) image
            std::size_t     iBest = ~std::size_t{0};
            std::uint64_t   Best  = 0;
            for (std::size_t i = 0, n = tinydng_image_count(pDoc); i < n; ++i)
            {
                const auto& Info = *tinydng_image_get(pDoc, i);
                if (Info.cfa.present) continue;
                if (const std::uint64_t Area = std::uint64_t{ Info.width } * Info.height; Area > Best) { Best = Area; iBest = i; }
            }
            if (iBest == ~std::size_t{0})
                return xerr::create_f<xbmp::tools::state, "The TIFF/DNG file only has raw sensor (CFA) images, they need debayering which is not supported">();

            tinydng_pixels Px{};
            if (tinydng_decode_image(Ctx.get(), pDoc, iBest, nullptr, &Px, &Err) != TINYDNG_OK)
            {
                xerr::LogMessage<state::FAILURE>(std::format("[Error]: Trying to decode TIFF/DNG: {} but found error {}", xstrtool::UTF16ToUTF8(FileName), Err.message));
                return xerr::create_f<xbmp::tools::state, "Fail to decode the TIFF/DNG image">();
            }
            auto PxFree = std::unique_ptr<tinydng_pixels, std::function<void(tinydng_pixels*)>>{ &Px, [&](tinydng_pixels* p) { tinydng_pixels_free(Ctx.get(), p); } };

            const int nChannels = Px.samples_per_pixel;
            const int Bits      = Px.bits_per_sample;
            const bool bFloat   = Px.sample_format == TINYDNG_SAMPLEFORMAT_IEEEFP;
            if (nChannels < 1 || nChannels > 4 || (Bits != 8 && Bits != 16 && Bits != 32) || (bFloat && Bits != 32))
                return xerr::create_f<xbmp::tools::state, "Unsupported TIFF/DNG pixel layout (channels, bits per sample or sample format)">();

            W = static_cast<int>(Px.width);
            H = static_cast<int>(Px.height);
            const std::size_t nPixels = std::size_t(W) * H;
            if (Px.size < nPixels * nChannels * (Bits / 8))
                return xerr::create_f<xbmp::tools::state, "The decoded TIFF/DNG image is smaller than its size says">();

            const bool bSigned = Px.sample_format == TINYDNG_SAMPLEFORMAT_INT;
            auto Sample = [&](std::size_t i) -> float
            {
                switch (Bits)
                {
                case 8:  return bSigned ? std::max(reinterpret_cast<const std::int8_t*>(Px.data)[i] / 127.0f, -1.0f) : Px.data[i] / 255.0f;
                case 16: return bSigned ? std::max(reinterpret_cast<const std::int16_t*>(Px.data)[i] / 32767.0f, -1.0f) : reinterpret_cast<const std::uint16_t*>(Px.data)[i] / 65535.0f;
                default: return bFloat  ? reinterpret_cast<const float*>(Px.data)[i]
                              : bSigned ? std::max(static_cast<float>(reinterpret_cast<const std::int32_t*>(Px.data)[i] / 2147483647.0), -1.0f)
                                        : static_cast<float>(reinterpret_cast<const std::uint32_t*>(Px.data)[i] / 4294967295.0);
                }
            };

            RGBA.resize(nPixels * 4);
            for (std::size_t p = 0; p < nPixels; ++p)
            {
                float C[4] = { 0, 0, 0, 1 };
                for (int c = 0; c < nChannels; ++c)
                    C[c] = Sample(p * nChannels + c);     // tinydng interleaves planar images on decode

                if (nChannels <= 2) { C[3] = nChannels == 2 ? C[1] : 1.0f; C[1] = C[2] = C[0]; }   // gray (+ alpha)

                std::copy_n(C, 4, &RGBA[p * 4]);
            }

            return {};
        }
    }

    //------------------------------------------------------------------------------

    xerr LoadTIFFImage(xbitmap& Bitmap, std::wstring_view FileName) noexcept
    {
        std::vector<float> RGBA;
        int W, H;
        if (auto Err = details_tiff::Decode(RGBA, W, H, FileName); Err) return Err;

        const auto DataSize      = std::size_t(W) * H * 4;
        const auto TotalDataSize = DataSize + sizeof(int);
        auto       FinalData     = std::make_unique<std::byte[]>(TotalDataSize);

        FinalData[0] = FinalData[1] = FinalData[2] = FinalData[3] = std::byte{ 0 };
        for (std::size_t i = 0; i < DataSize; ++i)
            FinalData[4 + i] = static_cast<std::byte>(static_cast<std::uint8_t>(std::clamp(RGBA[i], 0.0f, 1.0f) * 255.0f + 0.5f));

        Bitmap.setup(W, H, xbitmap::format::R8G8B8A8, DataSize, { FinalData.release(), TotalDataSize }, true, 1, 1);
        return {};
    }

    //------------------------------------------------------------------------------

    xerr LoadHDRTIFFImage(xbitmap& Bitmap, std::wstring_view FileName) noexcept
    {
        std::vector<float> RGBA;
        int W, H;
        if (auto Err = details_tiff::Decode(RGBA, W, H, FileName); Err) return Err;

        const auto DataSize      = RGBA.size() * sizeof(float);
        const auto TotalDataSize = DataSize + sizeof(int);
        auto       FinalData     = std::make_unique<std::byte[]>(TotalDataSize);

        FinalData[0] = FinalData[1] = FinalData[2] = FinalData[3] = std::byte{ 0 };
        std::memcpy(&FinalData[4], RGBA.data(), DataSize);

        Bitmap.setup(W, H, xbitmap::format::R32G32B32A32_FLOAT, DataSize, { FinalData.release(), TotalDataSize }, true, 1, 1);
        return {};
    }
}
