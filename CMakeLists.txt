cmake_minimum_required(VERSION 3.16)

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE INTERNAL FORCE)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(Interactive_ARAP)

# Options
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Libs CACHE PATH "Path to lib folder")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fsanitize=address")

# Executables
LIST(APPEND HEADER_FILES "include/cotangent_weight_matrix.h" "include/rotation_matrix.h")
LIST(APPEND SOURCE_FILES "src/main.cpp" "src/cotangent_weight_matrix.cpp" "src/rotation_matrix.cpp")

set(Eigen3_DIR ${LIBRARY_DIR}/Eigen/share/eigen3/cmake CACHE PATH "Path to installed Eigen")

include(libigl)

# Include optional libigl modules
#igl_include(module)
igl_include(glfw)

# Create a build target
add_executable(app ${SOURCE_FILES} ${HEADER_FILES})

#target_link_libraries(app igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
target_link_libraries(app Eigen3::Eigen igl::glfw)
target_include_directories(app PUBLIC ${EIGEN3_INCLUDE_DIR} "include")

# Visual Studio properties
#set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT app)
#set_property(TARGET app PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/)
