
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
	#include <intrin.h>
#else
	#include <x86intrin.h>
#endif

const char OUTPUT_FILE[] = "./results.bin";

typedef float float_t;

const float_t NX = 600.f;
const float_t NY = 600.f;
const float_t NT = 2500.f;
const float_t XA = 0.f;
const float_t XB = 4.0f;
const float_t YA = 0.0f;
const float_t YB = 4.0f;
const float_t tau = 0.01f;
const float_t f0 = 1.0f;
const float_t t0 = 1.5f;
const float_t gammaSqr = 16.0f;
const float_t SX = 1.;
const float_t SY = NY / 2;
const float_t twoPiF0 = M_PI * 2 * f0;

inline float_t max(float_t a, float_t b) {
  return a > b ? a : b;
}

inline float_t min(float_t a, float_t b) {
  return a < b ? a : b;
}

void find_min_max(float_t** U, int ny, int nx) {
  float_t mn = 10000;
  float_t mx = -10000;
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      mn = min(mn, U[i][j]);
      mx = max(mx, U[i][j]);
    }
  }
  printf("min=%.9f, max=%.9f\n", mn, mx);
}

float_t sqr(float_t a) {
  return a * a;
}

float_t source_function_coord(int i, int j, int n) {
  if (j != SX || i != SY) {
    return 0;
  }
  float_t arg = twoPiF0 * (n * tau - t0);
  return expf(-sqr(arg) / gammaSqr) * sinf(arg);
}

float_t** create(int ny, int nx) {
  float_t** arr = (float_t**)malloc(ny * sizeof(float_t*));
  for (int i = 0; i < ny; i++) {
    arr[i] = (float_t*)malloc(nx * sizeof(float_t));
    for (int j = 0; j < nx; ++j) {
      arr[i][j] = 0;
    }
  }
  return arr;
}

float_t** calc_P() {
  float_t** arr = (float_t**)malloc(NY * sizeof(float_t*));
  for (int i = 0; i < NY; i++) {
    arr[i] = (float_t*)malloc(NX * sizeof(float_t));
    for (int j = 0; j < NX; ++j) {
      arr[i][j] = j < NX / 2 ? 0.1f * 0.1f : 0.2f * 0.2f;
      ;
    }
  }
  return arr;  // задать матрицей
}

int main(void) {
  const float_t HXdeg2mul2 = 2 * sqr((XB - XA) / (NX - 1));
  const float_t HYdeg2mul2 = 2 * sqr((YB - YA) / (NY - 1));
  const float_t tauSqr = sqr(tau);

  float_t** U = create(NY, NX);
  float_t** uMinusOne = create(NY, NX);
  float_t** uPlusOne = create(NY, NX);
  const float_t** P = calc_P();
  unsigned long long before = __rdtsc();
  printf("before:%lld\n", before);

  for (int n = 0; n < NT; n++) {
    for (int i = 1; i < NX - 2; i++) {
      for (int j = 1; j < NY - 2; j++) {
        float_t fnij =
            source_function_coord(i, j, n);  //вынести массив одномерный

        float_t U1 = U[i][j + 1] - U[i][j];
        float_t P1 = P[i - 1][j] + P[i][j];

        float_t U2 = U[i][j - 1] - U[i][j];
        float_t P2 = P[i - 1][j - 1] + P[i][j - 1];

        float_t first = (U1 * P1 + U2 * P2) / HXdeg2mul2;

        float_t U3 = U[i + 1][j] - U[i][j];
        float_t P3 = P[i][j - 1] + P[i][j];

        float_t U4 = U[i - 1][j] - U[i][j];
        float_t P4 = P[i - 1][j - 1] + P[i - 1][j];

        float_t second = (U3 * P3 + U4 * P4) / HYdeg2mul2;

        uPlusOne[i][j] =
            2 * U[i][j] - uMinusOne[i][j] + tauSqr * (fnij + first + second);
      }
    }

    float_t** t = uMinusOne;
    uMinusOne = U;
    U = uPlusOne;
    uPlusOne = t;
    printf("\rstep:%d\\%.0f", n, NT);
    fflush(0);
    find_min_max(U, NY, NX);
  }
  printf("\n");

  unsigned long long after = __rdtsc();

  printf("nanoseconds :%f\n", (after - before) / 2.4f);

  FILE* output = NULL;

  output = fopen(OUTPUT_FILE, "wb");
  if (output == NULL) {
    printf("Error opening file");
    return (EXIT_FAILURE);
  }

  for (int i = 0; i < NY; ++i) {
    size_t i1 = fwrite(U[i], sizeof(float_t), NX, output);
  }

  return 0;
}
