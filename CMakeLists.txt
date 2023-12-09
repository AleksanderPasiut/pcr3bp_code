cmake_minimum_required(VERSION 3.10)
project(pcr3bp_code)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

################################################################################
# primary component build
################################################################################
set(SOURCES_LIST
    src/tools/variable_printer.cpp
    src/tools/print_bootstrap.cpp)

file(GLOB SOURCES_LIST_PROOF src/proof/*.cpp)
list(APPEND SOURCES_LIST ${SOURCES_LIST_PROOF})

add_executable(${PROJECT_NAME} ${SOURCES_LIST})

add_dependencies(${PROJECT_NAME} gtest)
add_dependencies(${PROJECT_NAME} capd_utils)

add_subdirectory(src/googletest)
add_subdirectory(src/capd_utils)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PRIVATE -lstdc++)
target_link_libraries(${PROJECT_NAME} PRIVATE -lm)
target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)

link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_link_libraries(${PROJECT_NAME} PRIVATE gtest)
target_link_libraries(${PROJECT_NAME} PRIVATE gtest_main)

# link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/src/capd_utils)
target_link_libraries(${PROJECT_NAME} PRIVATE capd_utils)


