cmake_minimum_required(VERSION 3.15)
project(hsimPost LANGUAGES C CXX)


# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 包含头文件目录
include_directories(${PROJECT_SOURCE_DIR})

# 添加dt库头文件路径
include_directories(${PROJECT_SOURCE_DIR}/extern/dt/include)

# 配置可执行文件
add_executable(${PROJECT_NAME}
    main.cpp
    common_type.h
    CLI11.hpp
	kdtree.c
	kdtree.h
	kdtree_m.h
	kdtree_m.cpp
)

# 添加头文件路径
include_directories(
    ./extern/eigen
)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# 添加头文件路径
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/extern/dt/include)


# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
