cmake_minimum_required(VERSION 3.10)

project(ReactorModel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined -g")

add_executable(ReactorModel
    src/ReactorEvents.cpp
)

target_include_directories(ReactorModel
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

add_subdirectory(libs/geometry_module)
target_link_libraries(ReactorModel PUBLIC geometry_module)

