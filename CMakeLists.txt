cmake_minimum_required(VERSION 4.0)
project(Ising_Model)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(srcs
  src/isingmodel.cpp
  src/utils.cpp
  src/main.cpp)
set(headers
  include/isingmodel.hpp)

add_executable(isingmodel)
target_sources(isingmodel PRIVATE ${srcs})
target_sources(isingmodel PUBLIC
  FILE_SET isingmodel_headers
  TYPE HEADERS
  BASE_DIRS ${PROJECT_SOURCE_DIR}/include
  FILES ${headers})

if (UNIX)
  set(CMAKE_EXPORT_COMPILE_COMMANDS on)
  target_compile_options(isingmodel PRIVATE -Wall -Wextra -Wpedantic)
else()
  target_compile_options(isingmodel PRIVATE /W4)
endif()

option(ISING_ENABLE_VISUALS "Enable the .show_simulation() method, which uses Raylib to render the simulation" OFF)
if (ISING_ENABLE_VISUALS)
  find_package(raylib REQUIRED)

  target_sources(isingmodel src/raylib.cpp)
  target_link_libraries(isingmodel PRIVATE raylib)
  target_compile_definitions(isingmodel PRIVATE ENABLE_RAYLIB=1)
endif()

function(copy_to_build_dir base)
  foreach(file IN LISTS ARGN)
    file(RELATIVE_PATH relative_path ${base} ${file})
    configure_file(${file} ${CMAKE_CURRENT_BINARY_DIR}/${relative_path} COPYONLY)
  endforeach()
endfunction()
file(GLOB plotscripts CONFIGURE_DEPENDS ./*.gplt)
copy_to_build_dir(${CMAKE_CURRENT_SOURCE_DIR} ${plotscripts})
