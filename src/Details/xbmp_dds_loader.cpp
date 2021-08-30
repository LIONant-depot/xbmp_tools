#define TINYDDSLOADER_IMPLEMENTATION
#include "../dependencies/tinyddsloader/tinyddsloader.h"

namespace xbmp::tools::loader {

    namespace dds
    {
        using namespace tinyddsloader;

        constexpr static auto error_codes_v = std::array<error*,8>
        { error_code< error::SUCCESS, error_str("Success") >
        , error_code< error::FAILURE, error_str("Error File Open") >
        , error_code< error::FAILURE, error_str("Error Read") >
        , error_code< error::FAILURE, error_str("Error Magic Word") >
        , error_code< error::FAILURE, error_str("Error Size") >
        , error_code< error::FAILURE, error_str("Error Verify") >
        , error_code< error::FAILURE, error_str("Error Not Supported") >
        , error_code< error::FAILURE, error_str("Error Invalid Data") >
        };

        //-------------------------------------------------------------------------------------

        std::optional<std::tuple<xcore::bitmap::format, xcore::bitmap::color_space, bool>> ConvertFormat( DDSFile::DXGIFormat Format ) noexcept
        {
            switch ( Format )
            {
                case DDSFile::DXGIFormat::BC1_UNorm:
                    return std::tuple{ xcore::bitmap::format::BC1_4RGBA1, xcore::bitmap::color_space::LINEAR, false };
                case DDSFile::DXGIFormat::BC1_UNorm_SRGB:
                    return std::tuple{ xcore::bitmap::format::BC1_4RGBA1, xcore::bitmap::color_space::SRGB, false };
                case DDSFile::DXGIFormat::BC2_UNorm:
                    return std::tuple{ xcore::bitmap::format::BC2_8RGBA, xcore::bitmap::color_space::LINEAR, false };
                case DDSFile::DXGIFormat::BC2_UNorm_SRGB:
                    return std::tuple{ xcore::bitmap::format::BC2_8RGBA, xcore::bitmap::color_space::SRGB, false };
                case DDSFile::DXGIFormat::BC3_UNorm:
                    return std::tuple{ xcore::bitmap::format::BC3_8RGBA, xcore::bitmap::color_space::LINEAR, false };
                case DDSFile::DXGIFormat::BC3_UNorm_SRGB:
                    return std::tuple{ xcore::bitmap::format::BC3_8RGBA, xcore::bitmap::color_space::SRGB, false };
                case DDSFile::DXGIFormat::R8G8B8A8_UNorm:
                    return std::tuple{ xcore::bitmap::format::R8G8B8A8, xcore::bitmap::color_space::LINEAR, false };
                case DDSFile::DXGIFormat::R8G8B8A8_UNorm_SRGB:
                    return std::tuple{ xcore::bitmap::format::R8G8B8A8, xcore::bitmap::color_space::SRGB, false };
                case DDSFile::DXGIFormat::R8G8B8A8_SNorm:
                    return std::tuple{ xcore::bitmap::format::R8G8B8A8, xcore::bitmap::color_space::LINEAR, true };
                case DDSFile::DXGIFormat::B8G8R8A8_UNorm:
                    return std::tuple{ xcore::bitmap::format::B8G8R8A8, xcore::bitmap::color_space::LINEAR, false };
                case DDSFile::DXGIFormat::B8G8R8A8_UNorm_SRGB:
                    return std::tuple{ xcore::bitmap::format::B8G8R8A8, xcore::bitmap::color_space::SRGB, false };
                case DDSFile::DXGIFormat::B8G8R8A8_Typeless:
                    return std::tuple{ xcore::bitmap::format::B8G8R8A8, xcore::bitmap::color_space::LINEAR, true };
                case DDSFile::DXGIFormat::B8G8R8X8_UNorm:
                    return std::tuple{ xcore::bitmap::format::B8G8R8U8, xcore::bitmap::color_space::LINEAR, false };
                case DDSFile::DXGIFormat::B8G8R8X8_UNorm_SRGB:
                    return std::tuple{ xcore::bitmap::format::B8G8R8U8, xcore::bitmap::color_space::SRGB, false };
                case DDSFile::DXGIFormat::B8G8R8X8_Typeless:
                    return std::tuple{ xcore::bitmap::format::B8G8R8U8, xcore::bitmap::color_space::LINEAR, true };
            }

            return {};
        }
    } // xbmp::tools::loader::dds

    //-------------------------------------------------------------------------------

    error* LoadDSS(xcore::bitmap& Bitmap, const char* pFileName) noexcept
    {
        using namespace tinyddsloader;

        DDSFile Image;

        if( auto Err = Image.Load(pFileName); Err )
        {
            return dds::error_codes_v[Err];
        }

        // Make sure we can handle the texture
        if( Image.GetTextureDimension() == DDSFile::TextureDimension::Unknown 
         || Image.GetTextureDimension() == DDSFile::TextureDimension::Texture3D )
        {
            return error_code< error::FAILURE, error_str("Unsupported dimension of texture") >;
        }

        // get basic information about the texture
        xcore::bitmap::format       Format;
        xcore::bitmap::color_space  ColorSpace;
        bool                        isSigned;
        if( auto E = dds::ConvertFormat( Image.GetFormat() ); E )
        {
            Format     = std::get<xcore::bitmap::format>(*E);
            ColorSpace = std::get<xcore::bitmap::color_space>(*E);
            isSigned   = std::get<bool>(*E);
        }
        else
        {
            return error_code< error::FAILURE, error_str("Unsupported texture format") >;
        }

        //
        // Prepare the memory
        //
        const auto MipTableBytes = Image.GetMipCount() * sizeof(xcore::bitmap::mip);
        const auto FaceByteSize  = [&]
        {
            auto FaceByteSize = 0;

            // Mips should be organized from biggest to smallest
            std::uint32_t PrevW = 0xffffffff;
            for( std::uint32_t i=0; i<Image.GetMipCount(); i++ )
            {
                auto View = Image.GetImageData(i, 0);
                FaceByteSize += View->m_memSlicePitch;
                if( View->m_width >= PrevW)
                {
                    return -1; 
                }
                PrevW = View->m_width;
            }

            return FaceByteSize;
        }();
        if(FaceByteSize == -1 ) return error_code< error::FAILURE, error_str("mips are organized incorrectly") >;

        const auto nSubFaces     = Image.IsCubemap() ? 6u : 1u;
        const auto FrameByteSize = FaceByteSize * nSubFaces;
        const auto nFrames       = Image.GetArraySize();
        const auto TotalByteSize = MipTableBytes + FrameByteSize * nFrames;

        //
        // Copy memory
        //
        auto Memory     = std::make_unique<std::byte[]>(TotalByteSize);
        auto pMipOffset = reinterpret_cast<xcore::bitmap::mip*>(Memory.get());
        auto pFrame     = reinterpret_cast<std::byte*>(&pMipOffset[Image.GetMipCount()]);

        // Set the very first offset
        pMipOffset[0].m_Offset = 0;

        for( std::uint32_t iFrame=0; iFrame < nFrames; ++iFrame)
        {
            for(std::uint32_t iSubFace=0; iSubFace < nSubFaces; ++iSubFace )
            {
                auto TopMipView = Image.GetImageData(0, iFrame * iSubFace);
                for (std::uint32_t iMip = 0; iMip < Image.GetMipCount(); ++iMip)
                {
                    auto View      = Image.GetImageData(iMip, iFrame * iSubFace );
                    auto ByteSize  = View->m_memSlicePitch;

                    // Set the offset of the next mip
                    if( iFrame == 0 && iSubFace == 0 ) 
                    {
                        if( (iMip +1) < Image.GetMipCount() ) pMipOffset[iMip+1].m_Offset = pMipOffset[iMip].m_Offset + static_cast<int>(ByteSize);
                    }
                    else
                    {
                        if( pMipOffset[iMip].m_Offset != static_cast<int>(ByteSize) )
                            return error_code< error::FAILURE, error_str("Unexcepted mipmap offset") >;
                    }

                    // Make sure that the size formula follows what we expect
                    if( std::max( 1u, (TopMipView->m_height>>iMip)) != View->m_height 
                     || std::max( 1u, (TopMipView->m_width >>iMip)) != View->m_width )
                    {
                        return error_code< error::FAILURE, error_str("Unexcepted mipmap size formulation") >;
                    }

                    // Copy the mip data
                    std::memcpy(&pFrame[pMipOffset[iMip].m_Offset + iSubFace * FaceByteSize + FrameByteSize * iFrame], View->m_mem, ByteSize);
                }
            }
        }

        //
        // Ready to setup the bitmap
        //
        Bitmap.setup
        ( Image.GetImageData(0,0)->m_width
        , Image.GetImageData(0,0)->m_height
        , Format
        , FrameByteSize
        , { Memory.release(), TotalByteSize }
        , true
        , Image.GetMipCount()
        , nFrames
        );

        return nullptr;
    }

//----------------------------------------------------------
// END
//----------------------------------------------------------
} // xbmp::tools::loader