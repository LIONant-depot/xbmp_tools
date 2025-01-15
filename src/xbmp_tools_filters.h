namespace xbmp::tools::filters
{
    // Sets the alpha to 0xff or 0xff based on the AlphaThreshold. 
    void ForcePunchThroughAlpha(xcore::bitmap& Bitmap, int AlphaThreshold) noexcept;

    // Computes the average color of the surrounding pixels that are not transparent
    // This is a way to fill the transparent pixels with the average color of the surrounding pixels
    // which can solved the black border issue when using alpha blending
    void FillAvrColorBaseOnAlpha(xcore::bitmap& Bitmap, const std::uint8_t AlphaThreshold = 128, std::uint32_t Depth = 0xffffffff) noexcept;

}    