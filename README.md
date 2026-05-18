# BlackScholes

Simple European call Black-Scholes pricing project.

Implements the Euler forward finite difference method to solve the Black-Scholes
PDE directly with example comparison against the analytical price.

## Building

[CMake] >=3.21 and a C++17 compiler is required to build the project.

For Unix:

<!-- Debug or default config by not setting CMAKE_BUILD_TYPE -->

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
```

For Windows:

<!-- Debug config by using Debug instead -->

```shell
cmake -S . -B build_win64 -A x64 && cmake --build build_win64 --config Release -j
```

[CTest] tests can be run to validate the build.

For Unix:

```shell
ctest --test-dir build -j$(nproc)
```

For Windows:

```shell
ctest --test-dir build_win64 -C Release -j%NUMBER_OF_PROCESSORS%
```
