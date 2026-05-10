add_executable(test_simulator test_simulator.cpp)
target_link_libraries(test_simulator PRIVATE planets_lib GTest::gtest_main GTest::gmock)

add_executable(test_plummer test_plummer.cpp)
target_link_libraries(test_plummer PRIVATE planets_lib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_simulator)
gtest_discover_tests(test_plummer)
