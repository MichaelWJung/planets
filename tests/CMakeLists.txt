add_executable(test_simulator test_simulator.cpp)
target_link_libraries(test_simulator PRIVATE planets_lib)

add_test(NAME SimulatorTests COMMAND test_simulator)
