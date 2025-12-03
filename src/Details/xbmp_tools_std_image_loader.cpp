#define STB_IMAGE_IMPLEMENTATION
#include "dependencies/stb/stb_image.h"
#include "dependencies/xstrtool/source/xstrtool.h"

namespace xbmp::tools::loader
{
    /*
    constexpr uint16_t floatToHalf(float f) noexcept
    {
        uint32_t    x       = *reinterpret_cast<uint32_t*>(&f);
        uint32_t    sign    = (x >> 16) & 0x8000;
        int32_t     exp     = ((x >> 23) & 0xFF) - 127 + 15;
        uint32_t    mant    = x & 0x007FFFFF;

        if (exp <= 0) 
        {
            if (exp < -10) 
            {
                return sign;
            }
            mant = (mant | 0x00800000) >> (1 - exp);
            return sign | (mant >> 13);
        }
        else if (exp == 0xFF - (127 - 15)) 
        {
            if (mant == 0) 
            {
                return sign | 0x7C00; // Inf
            }
            else 
            {
                mant >>= 13;
                return sign | 0x7C00 | mant | (mant == 0); // NaN
            }
        }
        else if (exp > 30) 
        {
            return sign | 0x7C00; // Overflow
        }

        return sign | (exp << 10) | (mant >> 13);
    }

    static constexpr auto float_to_half_v = []() consteval
    {
        std::array<uint16_t, 65536> lut = {};

        for (uint32_t i = 0; i < 65536; ++i)
        {
            uint32_t floatBits = i << 16;
            float f = *reinterpret_cast<float*>(&floatBits);
            lut[i] = floatToHalf(f);
        }

        return lut;
    }();

    constexpr std::uint16_t FastFloatToHalf( float x ) noexcept
    {
        uint32_t    i = *reinterpret_cast<uint32_t*>(&x);
        return float_to_half_v[i >> 16];
    }
    */

    constexpr float overflow()
    {
        assert(false);
        return 0;
        /*
        volatile float f = 1e10;

        for (int i = 0; i < 10; i++)
        {
            // this will overflow before
            // the for loop terminates
            f *= f; 
        }
        
        return f;
        */
    }

    consteval std::uint16_t FloatToHalf(const int i) noexcept
    {
        // Our floating point number, f, is represented by the bit
        // pattern in integer i.  Disassemble that bit pattern into
        // the sign, s, the exponent, e, and the significand, m.
        // Shift s into the position where it will go in in the
        // resulting half number.
        // Adjust e, accounting for the different exponent bias
        // of float and half (127 versus 15).
        int s = (i >> 16) & 0x00008000;
        int e = ((i >> 23) & 0x000000ff) - (127 - 15);
        int m = i & 0x007fffff;

        //
        // Now reassemble s, e and m into a half:
        //
        if (e <= 0) 
        {
            if (e < -10) 
            {
                // E is less than -10.  The absolute value of f is
                // less than HALF_MIN (f may be a small normalized
                // float, a denormalized float or a zero).
                //
                // We convert f to a half zero with the same sign as f.
                return s;
            }

            // E is between -10 and 0.  F is a normalized float
            // whose magnitude is less than HALF_NRM_MIN.
            //
            // We convert f to a denormalized half.
            //
            // Add an explicit leading 1 to the significand.
            m = m | 0x00800000;

            //
            // Round to m to the nearest (10+e)-bit value (with e between
            // -10 and 0); in case of a tie, round to the nearest even value.
            //
            // Rounding may cause the significand to overflow and make
            // our number normalized.  Because of the way a half's bits
            // are laid out, we don't have to treat this case separately;
            // the code below will handle it correctly.
            //
            int t = 14 - e;
            int a = (1 << (t - 1)) - 1;
            int b = (m >> t) & 1;

            m = (m + a + b) >> t;

            //
            // Assemble the half from s, e (zero) and m.
            //

            return s | m;
        }

        if (e == 0xff - (127 - 15)) 
        {
            if (m == 0) 
            {
                //
                // F is an infinity; convert f to a half
                // infinity with the same sign as f.
                //
                return s | 0x7c00;
            }

            // F is a NAN; we produce a half NAN that preserves
            // the sign bit and the 10 leftmost bits of the
            // significand of f, with one exception: If the 10
            // leftmost bits are all zero, the NAN would turn
            // into an infinity, so we have to set at least one
            // bit in the significand.

            m >>= 13;
            return s | 0x7c00 | m | (m == 0);
        }

        // E is greater than zero.  F is a normalized float.
        // We try to convert f to a normalized half.
        //
        // Round to m to the nearest 10-bit value.  In case of
        // a tie, round to the nearest even value.
        m = m + 0x00000fff + ((m >> 13) & 1);

        if (m & 0x00800000) 
        {
            // overflow in significand,
            m = 0;

            // adjust exponent
            e += 1;        
        }

        //
        // Handle exponent overflow
        //
        if (e > 30) 
        {
            // Cause a hardware floating point overflow;
            overflow();

            // if this returns, the half becomes an
            return s | 0x7c00;

            // infinity with the same sign as f.
        }

        //
        // Assemble the half from s, e and m.
        //
        return s | (e << 10) | (m >> 13);
    }


    xerr LoadSTDImage(xbitmap& Bitmap, std::wstring_view FileName) noexcept
    {
        int nChannels;
        int W, H;
        unsigned char* pData;

        if( pData = stbi_load(xstrtool::To(FileName).c_str(), &W, &H, &nChannels, 0); nullptr == pData)
            return xerr::create_f<xbmp::tools::state, "Fail to load image">();

        const auto DataSize      = (W * H * nChannels * sizeof(std::uint8_t));
        const auto TotalDataSize = DataSize + sizeof(int);

        auto FinalData = std::make_unique<std::byte[]>(TotalDataSize);

        FinalData[0] = std::byte{ 0 };
        FinalData[1] = std::byte{ 0 };
        FinalData[2] = std::byte{ 0 };
        FinalData[3] = std::byte{ 0 };
        std::memcpy( &FinalData[4], pData, DataSize );

        stbi_image_free(pData);


        xbitmap::format Fmt;
        switch(nChannels*8)
        {
        case 32: Fmt = xbitmap::format::R8G8B8A8; break;
        case 24: Fmt = xbitmap::format::R8G8B8; break;
        case 16: Fmt = xbitmap::format::R5G6B5; break;
        case  8: Fmt = xbitmap::format::R8; break;
        default:
            return xerr::create_f<xbmp::tools::state, "Unkown pixel depth (bits per pixel) in loaded image" >();
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

        return {};
    }

    xerr LoadHDRSTDImage(xbitmap& Bitmap, std::wstring_view FileName) noexcept
    {
        int         nChannels;
        int         W, H;
        float*      pData;

        if( pData = stbi_loadf( xstrtool::To(FileName).c_str(), &W, &H, &nChannels, 0); nullptr == pData)
            return xerr::create_f<xbmp::tools::state, "Fail to hdr load image" >(); 

        const auto DataSize      = (W * H * nChannels * sizeof(float));
        const auto TotalDataSize = DataSize + sizeof(int);

        auto FinalData = std::make_unique<std::byte[]>(TotalDataSize);

        FinalData[0] = std::byte{ 0 };
        FinalData[1] = std::byte{ 0 };
        FinalData[2] = std::byte{ 0 };
        FinalData[3] = std::byte{ 0 };
        std::memcpy( &FinalData[4], pData, DataSize );

        stbi_image_free(pData);

        xbitmap::format Fmt;
        switch(nChannels)
        {
        case 4: Fmt = xbitmap::format::R32G32B32A32_FLOAT; break;
        case 3: Fmt = xbitmap::format::R32G32B32_FLOAT;    break;
        default:
            return xerr::create_f<xbmp::tools::state, "Unkown pixel depth (bits per pixel) in loaded image" >();
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

        return {};
    }

}