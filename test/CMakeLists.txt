set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

add_executable(io_test io.cpp)
target_link_libraries(io_test Boost::unit_test_framework Threads::Threads EMPTool QuadSquad)

add_executable(sharing_test sharing.cpp)
target_link_libraries(sharing_test Boost::unit_test_framework Threads::Threads QuadSquad)

add_executable(utils_test utils.cpp)
target_link_libraries(utils_test Boost::unit_test_framework Threads::Threads QuadSquad)

add_executable(offline_test offline.cpp)
target_link_libraries(offline_test Boost::unit_test_framework Threads::Threads QuadSquad NTL GMP)

add_executable(jump_test jump.cpp)
target_link_libraries(jump_test Boost::unit_test_framework Threads::Threads QuadSquad)

add_executable(online_test online.cpp)
target_link_libraries(online_test Boost::unit_test_framework Threads::Threads QuadSquad)

add_custom_target(tests)
add_dependencies(tests
    io_test sharing_test utils_test offline_test jump_test online_test)

add_test(NAME io_test COMMAND io_test)
add_test(NAME sharing_test COMMAND sharing_test)
add_test(NAME utils_test COMMAND utils_test)
add_test(NAME offline_test COMMAND offline_test)
add_test(NAME jump_test COMMAND jump_test)
add_test(NAME online_test COMMAND online_test)
