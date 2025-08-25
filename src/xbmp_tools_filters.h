namespace xbmp::tools::filters
{
    // Sets the alpha to 0xff or 0xff based on the AlphaThreshold. 
    void ForcePunchThroughAlpha(xbitmap& Bitmap, int AlphaThreshold) noexcept;

    // Computes the average color of the surrounding pixels that are not transparent
    // This is a way to fill the transparent pixels with the average color of the surrounding pixels
    // which can solved the black border issue when using alpha blending
    void FillAvrColorBaseOnAlpha(xbitmap& Bitmap, const std::uint8_t AlphaThreshold = 128, std::uint32_t Depth = 0xffffffff) noexcept;

    // Makes the bitmap tilable by blending the edges based on a percentage
    void MakeBitmapTilable( xbitmap& Bitmap, float WidthOverlapPercentage = 0.3f, float HeightOverlapPercentage = 0.3f ) noexcept;
    void MakeBitmapTilableHDR(xbitmap& Bitmap, float WidthOverlapPercentage = 0.3f, float HeightOverlapPercentage = 0.3f) noexcept;

    void MakeBitmapTilableHEPDP(xbitmap& seamless, const xbitmap& src_image, float WidthOverlapPercentage = 0.15f, float HeightOverlapPercentage = 0.15f) noexcept;


    // Covert a 2D Image to a cubemap
    bool ConvertToCubeMapHDR( xbitmap& CubeMap, const xbitmap& Bitmap, int CubemapResolution, bool bUseBilinear = true ) noexcept;
    bool ConvertToCubeMap   ( xbitmap& CubeMap, const xbitmap& Bitmap, int CubemapResolution, bool bUseBilinear = true ) noexcept;
}    