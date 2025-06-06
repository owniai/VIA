cmake_minimum_required(VERSION 3.13)
project(VIA_PIR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")


add_executable(VIA_PIR ${SOURCES})


target_compile_options(VIA_PIR PRIVATE -mavx512dq -mavx512ifma)


# 指定头文件目录
target_include_directories(VIA_PIR PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# 指定库文件
target_link_libraries(VIA_PIR PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/libhexl.a)