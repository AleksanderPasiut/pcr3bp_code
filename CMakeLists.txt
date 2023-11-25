cmake_minimum_required(VERSION 3.10)
project(ursa)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

set(SOURCES_LIST
    src/tools/variable_printer.cpp
    src/tools/print_bootstrap.cpp
    src/main.cpp)

file(GLOB SOURCES_LIST_BASIC_TEST src/pcr3bp_basic_test/*.cpp)
list(APPEND SOURCES_LIST ${SOURCES_LIST_BASIC_TEST})

file(GLOB SOURCES_LIST_PROOF src/proof/*.cpp)
list(APPEND SOURCES_LIST ${SOURCES_LIST_PROOF})

add_executable(${PROJECT_NAME} ${SOURCES_LIST})

add_dependencies(${PROJECT_NAME} gtest)
add_dependencies(${PROJECT_NAME} carina)
add_dependencies(${PROJECT_NAME} taurus)

add_subdirectory(src/googletest)
add_subdirectory(src/carina)
add_subdirectory(src/bluesky)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PRIVATE -lstdc++)
target_link_libraries(${PROJECT_NAME} PRIVATE -lm)
target_link_libraries(${PROJECT_NAME} PRIVATE -lpthread)

link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/lib)
target_link_libraries(${PROJECT_NAME} PRIVATE gtest)
target_link_libraries(${PROJECT_NAME} PRIVATE gtest_main)

# link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/src/carina)
target_link_libraries(${PROJECT_NAME} PRIVATE carina)

target_link_libraries(${PROJECT_NAME} PRIVATE taurus)

################################################################################
# plot components build
################################################################################
function(plot_component SLAVE_NAME SOURCE_FILE)

    add_executable(${SLAVE_NAME} ${SOURCE_FILE})
    target_include_directories(${SLAVE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(${SLAVE_NAME} PRIVATE taurus)
    target_link_libraries(${SLAVE_NAME} PRIVATE carina)

endfunction()

if(TRUE)

plot_component( plot_reg_evolution_raw
    src/plot_1/reg_evolution_raw_plot.cpp)

plot_component( plot_reg_evolution_lyapunov_multi
    src/plot_1/reg_evolution_lyapunov_plot.multi.cpp)

endif()
