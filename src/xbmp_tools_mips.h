namespace xbmp::tools::mips
{
    inline
    int CalcMaxMipLevels(int width, int height, int min_size) 
    {
        int levels = 1;
        while (width > min_size || height > min_size) {
            width = std::max(1, width / 2);
            height = std::max(1, height / 2);
            levels++;
        }
        return levels;
    }

    inline
    int ComputeMaxMips(int MinSize, int orig_w, int orig_h)
    {
        int max_mips = 1;  // Always include level 0
        int prev_w = orig_w;
        int prev_h = orig_h;
        for (int l = 1; l < 32; ++l)
        {  // Safety limit
            int lw = std::max(1, orig_w >> l);
            int lh = std::max(1, orig_h >> l);
            if (lw == prev_w && lh == prev_h) break;  // No further change (e.g., at 1x1)
            if (std::min(lw, lh) < MinSize) break;
            max_mips++;
            prev_w = lw;
            prev_h = lh;
        }

        return max_mips;
    }

    enum class mips_filter_type
    { BOX
    , TRIANGLE
    , LANCZOS
    , KAISER
    };

    void GenerateMipMaps( xbitmap&          Bitmap
                        , const int         MinSize                 = 1
                        , const bool        isNormapMap             = false
                        , mips_filter_type  filter                  = mips_filter_type::KAISER
                        , bool              preserveAlphaCoverage   = false
                        , float             alphaCutoff             = 0.5f
                        , const xbitmap*    pNormalBitmap           = nullptr
                        , int               roughnessChannel        = 1
                        );
}