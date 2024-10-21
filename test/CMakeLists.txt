project(test_my_vegas)

CMAKE_MINIMUM_REQUIRED(VERSION 3.22
)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})


ADD_EXECUTABLE(RandomTest RandomGeneratorTest.cpp
)
ADD_EXECUTABLE(VEGASMAPTEST VEGAS_map_test.cpp)
ADD_EXECUTABLE(VEGASINTERTEST VEGAS_Inter_test.cpp)
#TARGET_LINK_LIBRARIES(RandomTest ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES(VEGASMAPTEST ${ABS_VEGAS_LIBRARIES})
#TARGET_LINK_LIBRARIES(VEGASINTERTEST ${ABS_VEGAS_LIBRARIES})

add_executable(test_my_vegas test_simple_functions.cpp)
add_executable(benchmark benchmark.cpp)
target_link_libraries(test_my_vegas  ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)
target_link_libraries(benchmark benchmark::benchmark)
