//
// tinydng (TIFF / DNG reader, https://github.com/syoyo/tinydng) compiled as one C translation unit.
// tinydng is C11, so it can not go into xbmp_tools.cpp with the other loaders; this file is compiled as C.
// Only the reading side is built (no writers, no PSD: stb_image reads those). Its private stb_image copy (JPEG inside TIFF/DNG) is
// made static, so it does not clash with the stb_image that xbmp_tools.cpp builds.
// tinydng's miniz is not built: other libraries in the same programs carry their own miniz (crunch in the texture compiler) and
// miniz's functions can not be made static, so they would be defined twice. tinydng only calls mz_uncompress (Deflate TIFF), which
// is given here by stb_image's zlib decoder, under a name of its own.
//
#define _CRT_SECURE_NO_WARNINGS
#define TINYDNG_NO_PSD
#define STB_IMAGE_STATIC
#define STBI_ONLY_PNG                       // keeps stb_image's zlib decoder (STBI_ONLY_JPEG alone compiles it out)
#define mz_uncompress xbmp_tools_tinydng_uncompress

#include "dependencies/tinydng/tinydng_api.c"
#include "dependencies/tinydng/tinydng_io.c"
#include "dependencies/tinydng/tinydng_tiff.c"
#include "dependencies/tinydng/tinydng_dng.c"
#include "dependencies/tinydng/tinydng_codec.c"
#include "dependencies/tinydng/tinydng_stb_image.c"
#include "dependencies/tinydng/tiny_dng_ljpeg92_v2.c"

int mz_uncompress(unsigned char* pDest, mz_ulong* pDest_len, const unsigned char* pSource, mz_ulong source_len)
{
    const int n = stbi_zlib_decode_buffer((char*)pDest, (int)*pDest_len, (const char*)pSource, (int)source_len);
    if (n < 0) return MZ_DATA_ERROR;
    *pDest_len = (mz_ulong)n;
    return MZ_OK;
}
