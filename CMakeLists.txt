cmake_minimum_required(VERSION 3.15)
project(hBLgenPre LANGUAGES CXX)

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
    boundingbox.h
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

# 设置dt库的路径（Debug和Release版本）
set(DT_LIB_PATH_DEBUG ${PROJECT_SOURCE_DIR}/extern/dt/Debug)
set(DT_LIB_PATH_RELEASE ${PROJECT_SOURCE_DIR}/extern/dt/Release)

# 添加头文件路径
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/extern/dt/include)

# 根据构建类型选择相应的库文件
#if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#    message(STATUS "Building in Debug mode")
#    target_link_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/extern/dt/Debug)
#    target_link_libraries(${PROJECT_NAME} dt)
#elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
#    message(STATUS "Building in Release mode")
#    target_link_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/extern/dt/Release)
#    target_link_libraries(${PROJECT_NAME} dt)
#endif()
target_link_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/extern/dt/Release)
target_link_libraries(${PROJECT_NAME} dt)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
