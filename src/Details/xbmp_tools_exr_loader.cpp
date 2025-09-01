// Must add to your includes: xcore\dependencies\zstd\lib

#define TINYEXR_USE_MINIZ 0
#define TINYEXR_USE_STB_ZLIB 1
#include "dependencies/zstd/lib/zstd.h"
#define TINYEXR_IMPLEMENTATION
#include "dependencies/tinyexr/tinyexr.h"
#include <format>

#include "dependencies/xstrtool/source/xstrtool.h"

namespace xbmp::tools::loader
{
    xerr LoadEXRImage(xbitmap& Bitmap, std::wstring_view FileName) noexcept
    {
        float*          pOut; // width * height * RGBA
        int             Width;
        int             Height;
        const char*     err     = nullptr;

        //
        // Load the image
        //
        if( auto ret = LoadEXR(&pOut, &Width, &Height, xstrtool::To(FileName).c_str(), &err);  ret != TINYEXR_SUCCESS)
        {
            // Display error message
            if (err) 
            {
                xerr::LogMessage<state::FAILURE>( std::format("[Error]: Trying to load EXR: {} but found error {}", xstrtool::To(FileName), err) );
                FreeEXRErrorMessage(err); // release memory of error message.
            }

            return xerr::create_f<state, "Unable to load the EXR image">();
        }

        //
        // Define some basic structures
        //
        struct rgb_f
        {
            float m_R;
            float m_G;
            float m_B;
            float m_A;
        };

        struct rgb_i
        {
            std::uint8_t m_R;
            std::uint8_t m_G;
            std::uint8_t m_B;
            std::uint8_t m_A;
        };

        //
        // Setup xbmp
        //
        const auto nPixels       = Width * Height;
        const auto FrameSize     = nPixels * sizeof(rgb_i);
        const auto TotalDataSize = FrameSize + sizeof(int);
        auto*      pData         = new std::byte[TotalDataSize];

        pData[0] = pData[1] = pData[2] = pData[3] = std::byte{0};

        // First set the data
        Bitmap.setup
        ( Width
        , Height
        , xbitmap::format::XCOLOR
        , FrameSize
        , { pData, TotalDataSize }
        , true
        , 1
        , 1
        );

        //
        // Now lets convert the data
        //
        auto pFOut = reinterpret_cast<const rgb_f*>(pOut);
        constexpr auto inv_gamma = 1.0 / 2.2;
        for( auto& E : Bitmap.getMip<rgb_i>(0) )
        {
            E.m_R = std::clamp(static_cast<std::uint32_t>(std::pow(pFOut->m_R, inv_gamma) * 255.0f), 0u, 0xffu );
            E.m_G = std::clamp(static_cast<std::uint32_t>(std::pow(pFOut->m_G, inv_gamma) * 255.0f), 0u, 0xffu );
            E.m_B = std::clamp(static_cast<std::uint32_t>(std::pow(pFOut->m_B, inv_gamma) * 255.0f), 0u, 0xffu );
            E.m_A = std::clamp(static_cast<std::uint32_t>(pFOut->m_A * 255.0f), 0u, 0xffu );
            pFOut++;
        }

        //
        // release memory of image data
        //
        free(pOut);

        return {};
    }

    xerr LoadHDREXRImage(xbitmap& Bitmap, std::wstring_view FileName) noexcept
    {
        float*          pOut; // width * height * RGBA
        int             Width;
        int             Height;
        const char*     err     = nullptr;

        //
        // Load the image
        //
        if( auto ret = LoadEXR(&pOut, &Width, &Height, xstrtool::To(FileName).c_str(), &err);  ret != TINYEXR_SUCCESS)
        {
            // Display error message
            if (err) 
            {
                xerr::LogMessage<state::FAILURE>(std::format("[Error]: Trying to load HDR EXR: {} but found error: {}", xstrtool::To(FileName), err));
                FreeEXRErrorMessage(err); // release memory of error message.
            }

            return xerr::create_f<state, "Unable to load the HDR EXR image">();
        }

        //
        // Setup xbmp
        //
        const auto nPixels       = Width * Height;
        const auto FrameSize     = nPixels * 4 * sizeof(float);
        const auto TotalDataSize = FrameSize + sizeof(int);
        auto*      pData         = new std::byte[TotalDataSize];

        pData[0] = pData[1] = pData[2] = pData[3] = std::byte{0};

        // copy the image data to the correct place
        memcpy(&pData[4], pOut, FrameSize);
        
        Bitmap.setup
        ( Width
        , Height
        , xbitmap::format::R32G32B32A32_FLOAT
        , FrameSize
        , { pData, TotalDataSize }
        , true
        , 1
        , 1
        );

        //
        // release memory of image data
        //
        free(pOut);

        return {};
    }
}



