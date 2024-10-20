//
// Created by andreas on 20.10.24.
//
#include "benchmark/benchmark.h"
#include "test_integrands.h"
static void BM_StringCreation(benchmark::State& state) {
    for (auto _ : state)
        std::string empty_string;
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation);

// Define another benchmark
static void BM_StringCopy(benchmark::State& state) {
    std::string x = "hello";
    for (auto _ : state)
        std::string copy(x);
}
BENCHMARK(BM_StringCopy);

static void BM_sin_cos_tan(benchmark::State& state) {
    constexpr int dimension = 3; // Set the number of dimensions
    std::array<double, dimension> x = {0.1, 0.5, 0.8}; // Example input array
    void* param = nullptr; // Example parameter (can be nullptr if unused)

    for (auto _ : state) {
        // Call the function you're benchmarking
        benchmark::DoNotOptimize(sin_cos_tan<dimension>(x, param));
    }
}

BENCHMARK(BM_sin_cos_tan);

BENCHMARK_MAIN();