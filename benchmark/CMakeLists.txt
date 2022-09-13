set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarks)

set(benchbin "")

macro(add_benchmark source_name)
    add_executable(${source_name} ${source_name}.cpp utils.cpp)
    target_link_libraries(${source_name} Boost::system Boost::program_options nlohmann_json::nlohmann_json QuadSquad Threads::Threads NTL GMP EMPTool)
    list(APPEND benchbin ${source_name})
endmacro()

add_benchmark(online_mpc)
add_benchmark(online_nn)
add_benchmark(offline_mpc_tp)
add_benchmark(offline_mpc_sub)
add_benchmark(sodo_gridlock_iter)

add_custom_target(benchmarks)
add_dependencies(benchmarks ${benchbin})
