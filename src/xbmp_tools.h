#ifndef _XBMP_TOOLS_H
#define _XBMP_TOOLS_H
#pragma once

#include"xcore.h"

#define XBMP_TOOLS_INLINE __forceinline

namespace xbmp::tools
{
    namespace details
    {
        template< typename T_ERR, std::size_t T_SIZE_V >
        struct error_info
        {
            static_assert(sizeof(T_ERR) == 1);
            T_ERR                           m_Code;
            std::array<char, T_SIZE_V>      m_Message;
            constexpr operator T_ERR*() const noexcept
            {
                return const_cast<T_ERR*>(&m_Code);
            }
            constexpr T_ERR* get( void ) const noexcept
            {
                return const_cast<T_ERR*>(&m_Code);
            }
        };
    }

    template< auto T_ERR_V, auto MSG >
    static constexpr auto error_code = details::error_info<decltype(T_ERR_V), MSG.size()>{ T_ERR_V, MSG };

    template< std::size_t T_SIZE_V >
    consteval auto error_str(const char(&Message)[T_SIZE_V] ) noexcept
    {
        std::array<char, T_SIZE_V> M;
        for( int i=0; M[i] = Message[i]; ++i );
        return M;
    }

    template<typename T> constexpr
    const char* getErrorMsg(T pErr) noexcept
    {
        assert(pErr != nullptr);
        return reinterpret_cast<const char*>(&pErr[1]);
    }

    template<typename T> constexpr
    int getErrorInt(T pErr) noexcept
    {
        return (pErr == nullptr) ? 0 : -(1 + static_cast<int>(*pErr));
    }

    template<typename ... Args> XBMP_TOOLS_INLINE
    std::string FormatString(const char* pFmt, Args ... args)
    {
        // Extra space for '\0' so + 1
        int size_s = std::snprintf(nullptr, 0, pFmt, args ...) + 1;
        if (size_s <= 0) { return "Error formatting string...."; }
        auto size = static_cast<size_t>(size_s);
        auto buf = std::make_unique<char[]>(size);
        std::snprintf(buf.get(), size, pFmt, args ...);
        // We don't want the '\0' inside so size -1
        return std::string(buf.get(), buf.get() + size - 1);
    }
}

#include "xbmp_tools_loaders.h"
#include "xbmp_tools_writers.h"

//---------------------------------------------------------------------------
// END
//---------------------------------------------------------------------------
#endif