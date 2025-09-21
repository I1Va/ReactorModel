cmake_minimum_required(VERSION 3.10)
project(ReactorModel LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_library(ReactorModel STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ReactorEvents.cpp
)

target_include_directories(ReactorModel
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

target_link_libraries(ReactorModel PUBLIC
    geometry_module
)

