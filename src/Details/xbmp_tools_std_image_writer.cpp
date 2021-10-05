
#define STB_IMAGE_WRITE_IMPLEMENTATION
#pragma warning( push )
#pragma warning( disable : 4996 )  // error C4996: 'sprintf': This function or variable may be unsafe.
#include "../dependencies/stb/stb_image_write.h"
#pragma warning( pop )

namespace xbmp::tools::writers
{
    error* SaveSTDImage( const char* pFileName, int W, int H, int Bpp, const std::byte* pData ) noexcept
    {
        if( xcore::string::FindStrI( pFileName, ".png") != -1 )
        {
            if( 0 == stbi_write_png(pFileName, W, H, Bpp/8, pData, W * Bpp / 8 ) )
                return error_code< error::FAILURE, error_str("Fail to write a PNG") >;
        }
        else if (xcore::string::FindStrI(pFileName, ".bmp") != -1)
        {
            if (0 == stbi_write_bmp(pFileName, W, H, Bpp / 8, pData ) )
                return error_code< error::FAILURE, error_str("Fail to write a BMP") >;
        }
        else if (xcore::string::FindStrI(pFileName, ".tga") != -1)
        {
            if (0 == stbi_write_tga(pFileName, W, H, Bpp / 8, pData ) )
                return error_code< error::FAILURE, error_str("Fail to write a TGA") >;
        }
        else if (xcore::string::FindStrI(pFileName, ".jpg") != -1)
        {
            if (0 == stbi_write_jpg(pFileName, W, H, Bpp / 8, pData, 90 ) )
                return error_code< error::FAILURE, error_str("Fail to write a JPG") >;
        }
        else if (xcore::string::FindStrI(pFileName, ".hdr") != -1)
        {
            // int stbi_write_hdr(pFileName, Bitmap.getWidth(), Bitmap.getHeight(), int comp, const float* data);
            return error_code< error::FAILURE, error_str("Unsupported file format") >;
        }
        else
        {
            return error_code< error::FAILURE, error_str("Unsupported file format") >;
        }

        return nullptr;
    }


    error* SaveSTDImage( const char* pFileName, const xcore::bitmap& Bitmap ) noexcept
    {
        xcore::color::format    ColorFmt = xcore::color::format{};
        switch( Bitmap.getFormat() )
        {
        case xcore::bitmap::format::R8G8B8A8:   ColorFmt = xcore::color::format{xcore::color::format::type::UINT_32_RGBA_8888}; break;
        case xcore::bitmap::format::R8G8B8:     ColorFmt = xcore::color::format{xcore::color::format::type::UINT_24_RGB_888  };   break;
        case xcore::bitmap::format::R5G6B5:     ColorFmt = xcore::color::format{xcore::color::format::type::UINT_16_RGB_565  };   break;
        }

        if( ColorFmt.m_Value == xcore::color::format::type::INVALID )
        {
            if( Bitmap.getFormat() >= xcore::bitmap::format::XCOLOR_END )
                return error_code< error::FAILURE, error_str("Image has the wrong format to the writters. Make sure it has a non-compress format.") >;

            auto  Data              = std::make_unique<std::byte[]>(Bitmap.getFrameSize());
            auto  ColorFmt          = xcore::color::format{ static_cast<xcore::color::format::type>( Bitmap.getFormat() ) };
            auto& Descriptor        = ColorFmt.getDescriptor();
            auto  pBitmapData       = Bitmap.getMip<std::byte>(0).data();
            const auto BypePerPixel = Descriptor.m_TB / 8;

            //
            // Default to 32bits
            //
            auto pData = reinterpret_cast<std::uint32_t*>(Data.get());
            for( int y=0, end_y = Bitmap.getHeight(); y < end_y; ++y )
            for( int x=0, end_x = Bitmap.getWidth();  x < end_x; ++x )
            {
                const std::uint32_t D = *reinterpret_cast<const std::uint32_t*>(pBitmapData);
                xcore::icolor       C = xcore::icolor{ D, ColorFmt };

                *pData = C.getDataFromColor({xcore::color::format::type::UINT_32_ARGB_8888});

                pData++;
                pBitmapData += BypePerPixel;
            }

            //
            // Save it
            //
            if( auto Err = SaveSTDImage(pFileName, Bitmap.getWidth(), Bitmap.getHeight(), 32, Data.get() ); Err )
                return Err;
        }
        else
        {
            auto  ColorFmt   = xcore::color::format{ static_cast<xcore::color::format::type>(Bitmap.getFormat()) };
            auto& Descriptor = ColorFmt.getDescriptor();

            //
            // Save it
            //
            if (auto Err = SaveSTDImage(pFileName, Bitmap.getWidth(), Bitmap.getHeight(), Descriptor.m_TB, Bitmap.getMip<std::byte>(0).data() ); Err)
                return Err;
        }

        return nullptr;
    }
}