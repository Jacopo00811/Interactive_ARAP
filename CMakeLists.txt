cmake_minimum_required(VERSION 3.16)

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE INTERNAL FORCE)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(Interactive_ARAP)

# Options
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Libs CACHE PATH "Path to lib folder")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fsanitize=address")

# Executables
# LIST(APPEND HEADER_FILES -)
LIST(APPEND SOURCE_FILES "main.cpp")

# Libraries
set(Eigen3_DIR ${LIBRARY_DIR}/Eigen/share/eigen3/cmake CACHE PATH "Path to installed Eigen")
#set(Libigl_INCLUDE_DIR ${LIBRARY_DIR}/libigl/include/ CACHE PATH "Path to libigl source folder")
#set(libigl_DIR ${LIBRARY_DIR}/libigl-2.5.0/cmake/igl CACHE PATH "Path to libigl")

# Libigl options
#option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
#option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
#option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
#option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)

# Finds and loads settings from the installed libraries
find_package(Eigen3 REQUIRED)
#find_package(LIBIGL REQUIRED)

include(libigl)

# Include optional libigl modules

#igl_include(module)
igl_include(glfw)

# Create a build target
#add_executable(app ${SOURCE_FILES} ${HEADER_FILES})
add_executable(app ${SOURCE_FILES})

#target_link_libraries(app igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
target_link_libraries(app Eigen3::Eigen igl::glfw)
#target_include_directories(app PUBLIC ${EIGEN3_INCLUDE_DIR} ${Libigl_INCLUDE_DIR})
target_include_directories(app PUBLIC ${EIGEN3_INCLUDE_DIR})

# Visual Studio properties
#set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT app)
#set_property(TARGET app PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/)
