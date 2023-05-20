#include "cpems/cpeds-math.h"

#include "gtest/gtest.h"

namespace cpems {
namespace {

template <typename T> T* getArray(T a0, long N) {
  T* A = new T[N];
  for (long i = 0; i < N; i++) {
    A[i] = a0 + i;
  }
  return A;
}

int sumLong(int a0, int N) {
  long* A = getArray<long>(a0, N);
  return cpeds_sum(A, N, true);
}

TEST(sumTest, canSumArray) { EXPECT_EQ(sumLong(1, 5), 15); }
// TEST(sumTest, canSumArray2) { EXPECT_EQ(sumLong(1, 5), 16); }

} // namespace
} // namespace cpems
