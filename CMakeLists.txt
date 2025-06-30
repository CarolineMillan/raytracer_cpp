cmake_minimum_required(VERSION 3.16)
project(RayTracer LANGUAGES CXX)

# 1) Require C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 2) Declare the executable target
add_executable(RayTracer)

# 3) Tell the compiler where to find headers
target_include_directories(RayTracer
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry
    ${CMAKE_CURRENT_SOURCE_DIR}/src/materials
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lights
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core
)

# 4) Automatically gather all .cpp files under src/
file(GLOB_RECURSE ALL_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# 5) Exclude none (we only globbed .cpp, so meshes/ .ply and test.ppm are safe)
# 6) Hook those .cpp files into your target
target_sources(RayTracer PRIVATE ${ALL_SOURCES})

target_compile_definitions(RayTracer
  PRIVATE
    _USE_MATH_DEFINES
)
target_compile_definitions(RayTracer
  PRIVATE
    RESOURCE_DIR=\"${CMAKE_SOURCE_DIR}/src/meshes\"
)


# 7) Optional: print out what’s being compiled
message(STATUS "Sources:")
foreach(src IN LISTS ALL_SOURCES)
  message(STATUS "  ${src}")
endforeach()
