#ifndef _XBMP_TOOLS_H
#define _XBMP_TOOLS_H
#pragma once

#include "source/xerr.h"
#include"source/xbitmap.h"
#include <memory>
#include <vector>
#include <cwctype> // for std::towlower


#define XBMP_TOOLS_INLINE __forceinline

namespace xbmp::tools
{
    enum class state : std::uint8_t
    { OK
    , FAILURE
    };

    //-----------------------------------------------------------------------------

    bool iequals(std::wstring_view lhs, std::wstring_view rhs)
    {
        if (lhs.size() != rhs.size()) return false;
        for (size_t i = 0; i < lhs.size(); ++i)
        {
            assert(std::towlower(rhs[i]) == rhs[i]);
            if (std::towlower(lhs[i]) != std::towlower(rhs[i]))
            {
                return false;
            }
        }
        return true;
    }

    //-----------------------------------------------------------------------------

    inline 
    std::string wstring_view_to_char(std::wstring_view wsv) noexcept
    {
        // Worst-case: each wchar_t could take up to 4 bytes in UTF-8
        std::string result;

        // Best guess reservation
        result.reserve(wsv.size() * 4); 

        for (wchar_t wc : wsv) 
        {
            // 1-byte UTF-8: ASCII range (U+0000 to U+007F)
            if (wc <= 0x7F) 
            {
                result += static_cast<char>(wc);
            }
            // 2-byte UTF-8: U+0080 to U+07FF
            else if (wc <= 0x7FF) 
            {
                result += static_cast<char>(0xC0 | ((wc >> 6) & 0x1F));
                result += static_cast<char>(0x80 | (wc & 0x3F));
            }
            // 3-byte UTF-8: U+0800 to U+FFFF (excluding surrogate pairs)
            else if (wc <= 0xFFFF) 
            {
                result += static_cast<char>(0xE0 | ((wc >> 12) & 0x0F));
                result += static_cast<char>(0x80 | ((wc >> 6) & 0x3F));
                result += static_cast<char>(0x80 | (wc & 0x3F));
            }
            // 4-byte UTF-8: U+10000 to U+10FFFF (for characters outside BMP)
            else if constexpr(sizeof(wchar_t) == 4)
            {
                result += static_cast<char>(0xF0 | ((unsigned int(wc) >> 18) & 0x07));
                result += static_cast<char>(0x80 | ((unsigned int(wc) >> 12) & 0x3F));
                result += static_cast<char>(0x80 | ((unsigned int(wc) >>  6) & 0x3F));
                result += static_cast<char>(0x80 | (unsigned int(wc) & 0x3F));
            }
        }

        return result;
    }

}

#include "xbmp_tools_loaders.h"
#include "xbmp_tools_writers.h"
#include "xbmp_tools_atlas.h"
#include "xbmp_tools_filters.h"

//---------------------------------------------------------------------------
// END
//---------------------------------------------------------------------------
#endif