// xbmp_tools unit tests: python -I D:\LIONant\xLION_Docs\tools\xbmp_tools_test.py (builds this with xbmp_tools.cpp and runs it)
#include "dependencies/xbmp_tools/src/xbmp_tools.h"
#include "dependencies/xbitmap/source/xbitmap.cpp"
#include "dependencies/xstrtool/source/xstrtool.cpp"
#include "dependencies/xbmp_tools/src/xbmp_tools.cpp"
#include "xbmp_tools_mips_unittest.h"

int main()
{
    return xbmp::tools::mips::unit_test::Test();
}
