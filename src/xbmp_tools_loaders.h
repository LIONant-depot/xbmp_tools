namespace xbmp::tools::loader
{
    xerr LoadDSS          ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;
    xerr LoadDSS          ( xbitmap& Bitmap, std::span<const std::byte> Buffer ) noexcept;
    xerr LoadSTDImage     ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;
    xerr LoadEXRImage     ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;
    xerr LoadHDRSTDImage  ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;
    xerr LoadHDREXRImage  ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;
    xerr LoadTIFFImage    ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;   // .tif .tiff .dng (8 bits RGBA)
    xerr LoadHDRTIFFImage ( xbitmap& Bitmap, std::wstring_view FileName ) noexcept;   // .tif .tiff .dng (32 bits float RGBA)
}