#include "dy4.h"
#include <benchmark/benchmark.h>
#include "fourier.h"
#include <vector>
#include "iofunc.h"
#include "utils_bench.h"
#include "filter.h"

#define RANGE_MULTIPLIER 2
#define MIN_INPUT_SIZE 128
#define MAX_INPUT_SIZE (8 * MIN_INPUT_SIZE)

const int Fs = 1000;
const int lower_bound = -1;
const int upper_bound = 1;

static void Bench_matrixPSD(benchmark::State& state) {
    int N = state.range(0);
    std::vector<real> samples(N);
    generate_random_values(samples, lower_bound, upper_bound);
    std::vector<real> freq, psd;

    for (auto _ : state) {
        matrixPSD(samples, Fs, freq, psd);
    }
}

BENCHMARK(Bench_matrixPSD)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);

static void Bench_estimatePSD(benchmark::State& state) {
    int N = state.range(0);
    std::vector<real> samples(N);
    generate_random_values(samples, lower_bound, upper_bound);
    std::vector<real> freq, psd;

    for (auto _ : state) {
        estimatePSD(samples, Fs, freq, psd);
    }
}

BENCHMARK(Bench_estimatePSD)->RangeMultiplier(RANGE_MULTIPLIER)->Range(MIN_INPUT_SIZE, MAX_INPUT_SIZE);
