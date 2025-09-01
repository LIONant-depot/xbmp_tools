
#define STB_IMAGE_WRITE_IMPLEMENTATION
#pragma warning( push )
#pragma warning( disable : 4996 )  // error C4996: 'sprintf': This function or variable may be unsafe.
#include "stb_image_write.h"
#pragma warning( pop )

#include "dependencies/xstrtool/source/xstrtool.h"

namespace xbmp::tools::writers
{
    //------------------------------------------------------------------------------

    xerr SaveSTDImage(std::wstring_view FileName, int W, int H, int Bpp, const std::byte* pData ) noexcept
    {
        if( xstrtool::findI(FileName, L".png") != std::string::npos )
        {
            if( 0 == stbi_write_png( xstrtool::To(FileName).c_str(), W, H, Bpp / 8, pData, W * Bpp / 8))
                return xerr::create_f<xbmp::tools::state, "Fail to write a PNG" >(); 
        }
        else if (xstrtool::findI(FileName, L".bmp") != std::string::npos )
        {
            if (0 == stbi_write_bmp(xstrtool::To(FileName).c_str(), W, H, Bpp / 8, pData ) )
                return xerr::create_f<xbmp::tools::state, "Fail to write a BMP" >();
        }
        else if (xstrtool::findI(FileName, L".tga") != std::string::npos )
        {
            if (0 == stbi_write_tga(xstrtool::To(FileName).c_str(), W, H, Bpp / 8, pData ) )
                return xerr::create_f<xbmp::tools::state, "Fail to write a TGA" >();
        }
        else if (xstrtool::findI(FileName, L".jpg") != std::string::npos )
        {
            if (0 == stbi_write_jpg(xstrtool::To(FileName).c_str(), W, H, Bpp / 8, pData, 90 ) )
                return xerr::create_f<xbmp::tools::state, "Fail to write a JPG" >();
        }
        else if (xstrtool::findI(FileName, L".hdr") != std::string::npos )
        {
            // int stbi_write_hdr(pFileName, Bitmap.getWidth(), Bitmap.getHeight(), int comp, const float* data);
            return xerr::create_f<xbmp::tools::state, "Unsupported file format" >(); 
        }
        else
        {
            return xerr::create_f<xbmp::tools::state, "Unsupported file format" >();
        }

        return {};
    }


    xerr SaveSTDImage( std::wstring_view FileName, const xbitmap& Bitmap ) noexcept
    {
        xcolor::format    ColorFmt = xcolor::format{};
        switch( Bitmap.getFormat() )
        {
        case xbitmap::format::R8G8B8A8:   ColorFmt = xcolor::format{xcolor::format::type::UINT_32_RGBA_8888}; break;
        case xbitmap::format::R8G8B8:     ColorFmt = xcolor::format{xcolor::format::type::UINT_24_RGB_888  };   break;
        case xbitmap::format::R5G6B5:     ColorFmt = xcolor::format{xcolor::format::type::UINT_16_RGB_565  };   break;
        }

        if( ColorFmt.m_Value == xcolor::format::type::INVALID )
        {
            if( Bitmap.getFormat() >= xbitmap::format::XCOLOR_END )
                return xerr::create_f<xbmp::tools::state, "Image has the wrong format to the writters. Make sure it has a non-compress format." >();

            auto  Data              = std::make_unique<std::byte[]>(Bitmap.getFrameSize());
            auto  ColorFmt          = xcolor::format{ static_cast<xcolor::format::type>( Bitmap.getFormat() ) };
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
                xcolori       C = xcolori{ D, ColorFmt };

                *pData = C.getDataFromColor({xcolor::format::type::UINT_32_ARGB_8888});

                pData++;
                pBitmapData += BypePerPixel;
            }

            //
            // Save it
            //
            if( auto Err = SaveSTDImage(FileName, Bitmap.getWidth(), Bitmap.getHeight(), 32, Data.get() ); Err )
                return Err;
        }
        else
        {
            auto  ColorFmt   = xcolor::format{ static_cast<xcolor::format::type>(Bitmap.getFormat()) };
            auto& Descriptor = ColorFmt.getDescriptor();

            //
            // Save it
            //
            if (auto Err = SaveSTDImage(FileName, Bitmap.getWidth(), Bitmap.getHeight(), Descriptor.m_TB, Bitmap.getMip<std::byte>(0).data() ); Err)
                return Err;
        }

        return {};
    }
}