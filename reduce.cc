#include "GaudiKernel/detected.h"
#include <benchmark/benchmark.h>
#include <boost/random/taus88.hpp>
#include <execution>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

boost::random::taus88 random_gen;

template <typename T>
using reduce_t = decltype(
    std::reduce(std::execution::unseq, std::declval<T>(), std::declval<T>()));

namespace detail {
template <typename T> auto minval_seq(T &&a, T &&b) {
  using std::min;
  return std::reduce(std::execution::seq, std::forward<T>(a),
                     std::forward<T>(b),
                     std::numeric_limits<typename T::value_type>::max(),
                     [](auto a, auto b) { return min(a, b); });
}

template <typename T> auto minval_vec(T &&a, T &&b) {
  using std::min;
  return std::reduce(std::execution::unseq, std::forward<T>(a),
                     std::forward<T>(b),
                     std::numeric_limits<typename std::decay_t<T>::value_type>::max(),
                     [](auto a, auto b) { return min(a, b); });
}

template <typename T> auto minval_par(T &&a, T &&b) {
  using std::min;
  return std::reduce(std::execution::par_unseq, std::forward<T>(a),
                     std::forward<T>(b),
                     std::numeric_limits<typename T::value_type>::max(),
                     [](auto a, auto b) { return min(a, b); });
}

template <typename T> auto minval_ser(T &&a, T &&b) {
  using std::min;
  return std::accumulate(
      std::forward<T>(a), std::forward<T>(b),
      std::numeric_limits<typename std::decay_t<T>::value_type>::max(),
      [](auto a, auto b) { return min(a, b); });
}
} // namespace detail

template <typename T> auto minval(T a, T b) {
  using std::min;
  if constexpr (Gaudi::cpp17::is_detected_v<reduce_t, T>) {
    return detail::minval_vec(a, b);
  } else {
    return detail::minval_ser(a, b);
  }
}

static void sequential(benchmark::State &state) {
  std::vector<float> store(256, 0);
  for (std::size_t i = 0; i < 256; ++i) {
    store[i] = random_gen();
  }
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        detail::minval_seq(std::begin(store), std::end(store)));
  }
}

static void serial(benchmark::State &state) {
  std::vector<float> store(256, 0);
  for (std::size_t i = 0; i < 256; ++i) {
    store[i] = random_gen();
  }
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        detail::minval_ser(std::begin(store), std::end(store)));
  }
}

static void vectorized(benchmark::State &state) {
  std::vector<float> store(256, 0);
  for (std::size_t i = 0; i < 256; ++i) {
    store[i] = random_gen();
  }
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        detail::minval_vec(std::begin(store), std::end(store)));
  }
}

static void defaulted(benchmark::State &state) {
  std::vector<float> store(256, 0);
  for (std::size_t i = 0; i < 256; ++i) {
    store[i] = random_gen();
  }
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        minval(std::begin(store), std::end(store)));
  }
}

static void parallel(benchmark::State &state) {
  std::vector<float> store(256, 0);
  for (std::size_t i = 0; i < 256; ++i) {
    store[i] = random_gen();
  }
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        detail::minval_par(std::begin(store), std::end(store)));
  }
}

auto compute_min = [](const std::vector<double> &v) -> double {
  return *(std::min_element(std::begin(v), std::end(v)));
};

BENCHMARK(parallel)
    ->ComputeStatistics("min", compute_min)
    ->UseRealTime()
    ->ThreadRange(1, 4);
BENCHMARK(sequential)
    ->ComputeStatistics("min", compute_min)
    ->UseRealTime()
    ->ThreadRange(1, 4);
BENCHMARK(vectorized)
    ->ComputeStatistics("min", compute_min)
    ->UseRealTime()
    ->ThreadRange(1, 4);
BENCHMARK(serial)
    ->ComputeStatistics("min", compute_min)
    ->UseRealTime()
    ->ThreadRange(1, 4);
BENCHMARK(defaulted)
    ->ComputeStatistics("min", compute_min)
    ->UseRealTime()
    ->ThreadRange(1, 4);

BENCHMARK_MAIN();