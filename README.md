# QuadSquad

This directory contains the implementation of the QuadSquad fair protocol.
The protocol is implemented in C++17 and [CMake](https://cmake.org/) is used as the build system.

## External Dependencies
The following libraries need to be installed separately and should be available to the build system and compiler.

- [GMP](https://gmplib.org/)
- [NTL](https://www.shoup.net/ntl/) (11.0.0 or later)
- [Boost](https://www.boost.org/) (1.72.0 or later)
- [Nlohmann JSON](https://github.com/nlohmann/json)
- [EMP Tool](https://github.com/emp-toolkit/emp-tool)
- [EMP OT](https://github.com/emp-toolkit/emp-ot/)

### Docker
All required dependencies to compile and run the project are available through the docker image.
To build and run the docker image, execute the following commands from the root directory of the repository:

```sh
# Build the QuadSquad Docker image.
#
# Building the Docker image requires at least 4GB RAM. This needs to be set 
# explicitly in case of Windows and MacOS.
docker build -t quadsquad .

# Create and run a container.
#
# This should start the shell from within the container.
docker run -it -v $PWD:/code quadsquad

# The following command changes the working directory to the one containing the 
# source code and should be run on the shell started using the previous command.
cd /code
```

## Compilation
The project uses [CMake](https://cmake.org/) for building the source code. 
To compile, run the following commands from the root directory of the repository:

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..

# The two main targets are 'benchmarks' and 'tests' corresponding to
# binaries used to run benchmarks and unit tests respectively.
make <target>
```

## Usage
A short description of the compiled programs is given below.
All of them provide detailed usage description on using the `--help` option.

- `benchmarks/online_mpc`: Benchmark the performance of the QuadSquad online phase by evaluating a circuit with a given depth and number of multiplication gates at each depth.
- `benchmarks/online_nn`: Benchmark the performance of the QuadSquad online phase for neural network inference on the FCN and LeNet models.
- `benchmarks/offline_mpc_tp`: Benchmark the performance of the QuadSquad offline phase for a circuit with given number of multiplication gates.
- `benchmarks/offline_mpc_sub`: Benchmark the performance of the subprotocols used in the QuadSquad offline phase.
- `benchmarks/sodo_gridlock_iter`: Benchmark the performance of an iteration of the sodoGR protocol of AST22 for liquidity matching.
- `tests/*`: These programs contain unit tests for various parts of the codebase. The test coverage is currently incomplete. However, the protocols have been manually verified for correctness.

Execute the following commands from the `build` directory created during compilation to run the programs:
```sh
# Ferret OT requires the 'ot_data' directory to store pre-OT data.
mkdir ot_data

# Run unit tests. Can be skipped.
#
# 'offline_test' is known to fail sometimes (especially on MacOS) because of 
# running multiple instances of Ferret OT in different threads. However, this
# is not an issue for benchmark programs.
ctest

# Benchmark online phase for MPC.
#
# The command below should be run on four different terminals with $PID set to
# 0, 1, 2, and 3 i.e., one instance corresponding to each party.
#
# The number of threads can be set using the '-t' option. '-g' denotes the 
# number of gates at each level and '-d' denotes the depth of the circuit.
#
# The program can be run on different machines by replacing the `--localhost`
# option with '--net-config <net_config.json>' where 'net_config.json' is a
# JSON file containing the IPs of the parties. A template is given in the
# repository root.
./benchmarks/online_mpc -p $PID --localhost -g 100 -d 10

# The `run.sh` script in the repository root can be used to run the programs 
# for all parties from the same terminal.
# For example, the previous benchmark can be run using the script as shown
# below.
../run.sh ./benchmarks/online_mpc -g 100 -d 10

# All other benchmark programs have similar options and behaviour. The '-h'
# option can be used for detailed usage information.

# Benchmark online phase for neural network inference.
../run.sh ./benchmarks/online_nn -n fcn

# Benchmark offline phase for MPC.
../run.sh ./benchmarks/offline_mpc_tp -g 1024

# Benchmark subprotocols of offline phase for MPC.
../run.sh ./benchmarks/offline_mpc_sub -g 1024

# Benchmark an iteration of the sodoGR protocol of AST22 for liquidity matching.
../run.sh ./benchmarks/sodo_gridlock_iter -b 256 -x 1000
```
