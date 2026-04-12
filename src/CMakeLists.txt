cmake_minimum_required(VERSION 3.16)
project(StumpDemo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------- Eigen ----------------
find_package(Eigen3 REQUIRED)

# ---------------- OpenMP ----------------
find_package(OpenMP REQUIRED)

# ---------------- Library ----------------
add_library(stump
    Stump.cpp
)

target_include_directories(stump
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ---------------- METIS ----------------
find_path(METIS_INCLUDE_DIR metis.h
    PATHS
        /usr/include
        /usr/local/include
)

find_library(METIS_LIBRARY metis
    PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
)

if (NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
    message(FATAL_ERROR "METIS not found. Please install libmetis-dev.")
endif()

target_include_directories(stump PUBLIC ${METIS_INCLUDE_DIR})
target_link_libraries(stump PUBLIC ${METIS_LIBRARY})

# ---------------- AMD (SuiteSparse) ----------------
find_path(AMD_INCLUDE_DIR amd.h
    PATHS
        /usr/include
        /usr/include/suitesparse
        /usr/local/include
        /usr/local/include/suitesparse
)

find_library(AMD_LIBRARY NAMES amd libamd
    PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
)

if (NOT AMD_INCLUDE_DIR OR NOT AMD_LIBRARY)
    message(FATAL_ERROR "AMD not found. Please install libsuitesparse-dev.")
endif()

target_include_directories(stump PUBLIC ${AMD_INCLUDE_DIR})
target_link_libraries(stump PUBLIC ${AMD_LIBRARY})

# stump uses Eigen types in public interfaces (or source), keep it linked
target_link_libraries(stump PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

# ---------------- MATIO (.mat file I/O) ----------------
find_path(MATIO_INCLUDE_DIR matio.h
    PATHS
        /usr/include
        /usr/local/include
)

find_library(MATIO_LIBRARY NAMES matio libmatio
    PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/local/lib
)

if (NOT MATIO_INCLUDE_DIR OR NOT MATIO_LIBRARY)
    message(FATAL_ERROR "MATIO not found. Please install libmatio-dev.")
endif()

# ---------------- Demo executable ----------------
add_executable(stump_demo main.cpp)
target_link_libraries(stump_demo PRIVATE stump ${MATIO_LIBRARY})
target_include_directories(stump_demo PRIVATE ${MATIO_INCLUDE_DIR})

# 将矩阵文件目录以宏传入，运行时无需关心 build 目录位置
target_compile_definitions(stump_demo PRIVATE
    DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../matrices"
)