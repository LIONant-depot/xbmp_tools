#ifndef _XBMP_TOOLS_H
#define _XBMP_TOOLS_H
#pragma once

#include "dependencies/xerr/source/xerr.h"
#include "dependencies/xbitmap/source/xbitmap.h"
#include <memory>
#include <vector>
#include <cwctype> // for std::towlower
#include <functional>

#include "dependencies/xmath/source/xmath.h"

#define XBMP_TOOLS_INLINE __forceinline
namespace xbmp::tools
{
    enum class state : std::uint8_t
    { OK
    , FAILURE
    };
}

#include "xbmp_tools_loaders.h"
#include "xbmp_tools_writers.h"
#include "xbmp_tools_atlas.h"
#include "xbmp_tools_filters.h"
#include "xbmp_tools_lighting.h"
#include "xbmp_tools_mips.h"

//---------------------------------------------------------------------------
// END
//---------------------------------------------------------------------------
#endif