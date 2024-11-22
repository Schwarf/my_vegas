//
// Created by andreas on 20.10.24.
//
#include "benchmark/benchmark.h"
#include "integrands.h"
#include "VEGAS_Integrator.h"
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


static void BM_VegasPolynom2Integration(benchmark::State& state) {
    constexpr int dimensions = 5;

    // VegasNumericalIntegration object setup
    VegasNumericalIntegration<dimensions> integrator;
    integrator.set_verbosity(VegasVerbosity::None);
    integrator.set_integrand(std::move(polynom2<dimensions>), nullptr);

    // Benchmark loop
    for (auto _ : state) {
        // Improve the integration grid and run the integration
        integrator.improve_grid();
        integrator.integrate();

        // Prevent compiler optimizations from affecting the timing
        benchmark::DoNotOptimize(integrator.get_result());

        // Optionally, print the result and error if needed (you might not want this for every benchmark iteration)
        std::cout << integrator.get_result() << " +/- " << integrator.get_error()
                  << " with chi-square: " << integrator.get_chisquare() << std::endl;
    }
}

// Register the benchmark
BENCHMARK(BM_VegasPolynom2Integration) ->Iterations(10);

BENCHMARK(BM_sin_cos_tan);

BENCHMARK_MAIN();