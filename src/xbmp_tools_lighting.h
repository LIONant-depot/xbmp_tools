namespace xbmp::tools::lighting
{
    void GenerateGGX_BRDF_RG_32FLOAT ( xbitmap& Bitmap, int Width = 512, int Height = 512) noexcept;
    void GenerateGGX_BRDF_RG_16SFLOAT( xbitmap& Bitmap, int Width = 512, int Height = 512) noexcept;


    enum class diffuse_type { NONE, CUBEMAP, SH };
    void PrefilterCubemap
    ( const xbitmap&                inputCubemap
    , xbitmap&                      specularOutput
    , std::function<void(float)>    progressCallback    = {}
    , int                           nMips               = -1
    , int                           specularSamples     = 1024
    , diffuse_type                  diffuse_type        = diffuse_type::NONE
    , xbitmap*                      diffuseCubemap      = nullptr
    , std::vector<xmath::fvec3>*    shCoeffs            = nullptr
    , int                           diffuseRes          = 32
    , int                           diffuseSamples      = 4096
    );

}