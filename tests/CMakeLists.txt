add_executable(test_simulator test_simulator.cpp)
target_link_libraries(test_simulator PRIVATE planets_lib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_simulator)
