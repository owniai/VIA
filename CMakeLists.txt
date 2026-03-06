cmake_minimum_required(VERSION 3.13)
project(VIA_PIR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")


add_executable(VIA_PIR ${SOURCES})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

target_compile_options(VIA_PIR PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-O3 -DNDEBUG -march=native>
    $<$<CXX_COMPILER_ID:MSVC>:/O2 /DNDEBUG>
)

target_compile_options(VIA_PIR PRIVATE -mavx512dq -mavx512ifma)


target_include_directories(VIA_PIR PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(VIA_PIR PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/libhexl.a)
