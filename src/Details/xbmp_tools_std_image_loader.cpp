#define STB_IMAGE_IMPLEMENTATION
#include "../dependencies/stb/stb_image.h"

namespace xbmp::tools::loader
{
    error* LoadSTDImage(xcore::bitmap& Bitmap, const char* pFileName) noexcept
    {
        int nChannels;
        int W, H;
        unsigned char* pData;

        if( pData = stbi_load(pFileName, &W, &H, &nChannels, 0 ); nullptr == pData )
            return error_code< error::FAILURE, error_str("Fail to load image") >;

        const auto DataSize      = (W * H * nChannels);
        const auto TotalDataSize = DataSize + sizeof(int);

        auto FinalData = std::make_unique<std::byte[]>(TotalDataSize);

        FinalData[0] = std::byte{ 0 };
        FinalData[1] = std::byte{ 0 };
        FinalData[2] = std::byte{ 0 };
        FinalData[3] = std::byte{ 0 };
        std::memcpy( &FinalData[4], pData, DataSize );

        stbi_image_free(pData);


        xcore::bitmap::format Fmt;
        switch(nChannels*8)
        {
        case 32: Fmt = xcore::bitmap::format::R8G8B8A8; break;
        case 24: Fmt = xcore::bitmap::format::R8G8B8; break;
        case 16: Fmt = xcore::bitmap::format::R5G6B5; break;
        default:
            return error_code< error::FAILURE, error_str("Unkown pixel depth (bits per pixel) in loaded image") >;
        };

        Bitmap.setup
        ( W
        , H
        , Fmt
        , DataSize
        , { FinalData.release(), TotalDataSize }
        , true
        , 1
        , 1
        );

        return nullptr;
    }

}