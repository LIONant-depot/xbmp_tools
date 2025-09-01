include(${CMAKE_BINARY_DIR}/_deps/xcmake_tools/Common.cmake)

# Fetch dependencies
FetchAndPopulate("https://github.com/benikabocha/tinyddsloader.git" "master")
FetchAndPopulate("https://github.com/nothings/stb.git"              "master")
FetchAndPopulate("https://github.com/syoyo/tinyexr.git"             "release")
FetchAndPopulate("https://github.com/LIONant-depot/xbitmap.git"     "main")
FetchAndPopulate("https://github.com/LIONant-depot/xmath.git"       "main")
FetchAndPopulate("https://github.com/LIONant-depot/xstrtool.git"    "main")

DefineInterfaceComponent(xbmp_tools "dependencies/xcore"
  "src/xbmp_tools.cpp"
  "src/xbmp_tools.h"
  "src/xbmp_tools_atlas.h"
  "src/xbmp_tools_filters.h"
  "src/xbmp_tools_loaders.h"
  "src/xbmp_tools_writers.h"
)
# We could add these here but unless it is required by other complex dependencies 
# it is better to let the user to add the necesary pathing by themselves... 
# if(xbmp_tools_CREATED)
#   AppendComponentIncludes(xbmp_tools
#     "${CMAKE_SOURCE_DIR}/dependencies/stb"
#     "${CMAKE_SOURCE_DIR}/dependencies/tinyddsloader"
#     "${CMAKE_SOURCE_DIR}/dependencies/tinyexr"
#   )
#   endif()

