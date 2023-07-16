cmake_minimum_required(VERSION 3.10)
project(ursa)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

set(SOURCES_LIST
    src/tools/variable_printer.cpp
    src/tools/print_bootstrap.cpp)

file(GLOB SOURCES_LIST_BASIC_TEST src/pcr3bp_basic_test/*.cpp)
list(APPEND SOURCES_LIST ${SOURCES_LIST_BASIC_TEST})

file(GLOB SOURCES_LIST_PROOF src/proof/*.cpp)
list(APPEND SOURCES_LIST ${SOURCES_LIST_PROOF})

add_executable(${PROJECT_NAME} ${SOURCES_LIST})

add_dependencies(${PROJECT_NAME} gtest)
add_dependencies(${PROJECT_NAME} carina)

add_subdirectory(src/googletest)
add_subdirectory(src/carina)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PRIVATE -lstdc++)
target_link_libraries(${PROJECT_NAME} PRIVATE -lm)
target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)

link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_link_libraries(${PROJECT_NAME} PRIVATE gtest)
target_link_libraries(${PROJECT_NAME} PRIVATE gtest_main)

link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/src/carina)
target_link_libraries(${PROJECT_NAME} PRIVATE carina)
