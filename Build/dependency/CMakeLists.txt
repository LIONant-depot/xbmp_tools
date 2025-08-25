include(${CMAKE_BINARY_DIR}/_deps/xcmake_tools/Common.cmake)

# Fetch dependencies
FetchAndPopulate("https://github.com/benikabocha/tinyddsloader.git" "master")
FetchAndPopulate("https://github.com/nothings/stb.git" "master")
FetchAndPopulate("https://github.com/syoyo/tinyexr.git" "release")
FetchAndPopulate("https://github.com/LIONant-depot/xbitmap.git")

DefineInterfaceComponent(xbmp_tools "dependencies/xcore"
  "src/xbmp_tools.cpp"
  "src/xbmp_tools.h"
  "src/xbmp_tools_atlas.h"
  "src/xbmp_tools_filters.h"
  "src/xbmp_tools_loaders.h"
  "src/xbmp_tools_writers.h"
)
if(xbmp_tools_CREATED)
  AppendComponentIncludes(xbmp_tools
    "${CMAKE_SOURCE_DIR}/dependencies/stb"
    "${CMAKE_SOURCE_DIR}/dependencies/tinyddsloader"
    "${CMAKE_SOURCE_DIR}/dependencies/tinyexr"
  )
  endif()

