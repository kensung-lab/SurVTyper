cmake_minimum_required(VERSION 3.5)
project(SurVTyper)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.13/include)
include_directories(libs)

add_executable(genotype_svs libs/ssw.c libs/ssw_cpp.cpp genotype_svs.cpp)

find_library(HTS_LIB hts PATHS htslib-1.13/lib NO_DEFAULT_PATH)
target_link_libraries(genotype_svs "${HTS_LIB}")
