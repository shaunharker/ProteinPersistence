cmake_minimum_required(VERSION 2.8.12)
project(ProteinPersistence)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})
add_subdirectory(pybind11)

include_directories ( 
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/phat/include
  /usr/local/include 
  ${USER_INCLUDE_PATH} )

link_directories ( 
  ${USER_LIBRARY_PATH} 
  /usr/local/lib )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

pybind11_add_module(ProteinPersistence src/pdb2persistence.cpp)
add_subdirectory (src)
