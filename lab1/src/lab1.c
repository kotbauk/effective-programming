
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const char OUTPUT_FILE[] = "results.bin";

typedef float float_t;
float_t tau = 0.01;
const float_t XA = 0.f;
const float_t XB = 4.0f;
const float_t YA = 0.0f;
const float_t YB = 4.0f;
const float_t f0 = 1.0f;
const float_t t0 = 1.5f;
const float_t gamma = 0.25;
const uint32_t SX = 1;
const float_t twoPiF0 = M_PI * 2 * f0;
const float_t t0twoPiF0 = t0 * twoPiF0;
const float_t t0twoPi = t0 * 2 * M_PI;

float_t* create(int ny, int nx);
float_t* calc_P(int ny, int nx);
__attribute__((always_inline)) inline float_t sqr(const float_t a);
int main(int argc, const char* argv[]) {
  const uint32_t NX = atoi(argv[1]);
  const uint32_t NY = atoi(argv[2]);
  const uint32_t NT = atoi(argv[3]);
  const uint32_t SY = NY / 2;
  const uint64_t xBorder = NX - 2;
  const uint64_t yBorder = NY - 2;
  const float_t HXdeg2mul2 = 2 * sqr((XB - XA) / (NX - 1));
  const float_t HYdeg2mul2 = 2 * sqr((YB - YA) / (NY - 1));

#ifdef DEBUG
  FILE* max_u_file = fopen("max_u_file", "wb");
#endif

  if (NX > 1000 || NY > 1000) {
    tau = 0.001;
  }

  const float_t tauSqr = sqr(tau);

  float_t* U = create(NY, NX);
  float_t* uMinusOne = create(NY, NX);
  float_t* uPlusOne = create(NY, NX);
  float_t* P = calc_P(NY, NX);

#ifdef DEBUG
  float_t temp_max = -0.1;
#endif
  const float_t tautwoPiF0 = tau * t0twoPiF0;
  float_t arg = -tautwoPiF0 - t0twoPi;
  for (int n = 0; n < NT; n++) {
    for (int i = 1; i < yBorder; i++) {
      for (int j = 1; j < xBorder; j++) {
        float_t U1 = U[j + 1 + i * (NX)] - U[j + i * (NX)];
        float_t P1 = P[j + (i - 1) * (NX)] + P[j + i * (NX)];

        float_t U2 = U[j - 1 + i * (NX)] - U[j + i * (NX)];
        float_t P2 = P[(j - 1) + (i - 1) * NX] + P[(j - 1) + i * (NX)];

        float_t first = (U1 * P1 + U2 * P2) / HXdeg2mul2;

        float_t U3 = U[j + (i + 1) * (NX)] - U[j + i * (NX)];
        float_t P3 = P[(j - 1) + i * (NX)] + P[j + i * (NX)];

        float_t U4 = U[j + (i - 1) * (NX)] - U[j + i * (NX)];
        float_t P4 = P[j - 1 + (i - 1) * (NX)] + P[j + (i - 1) * (NX)];

        float_t second = (U3 * P3 + U4 * P4) / HYdeg2mul2;

        uPlusOne[j + i * (NX)] = 2 * U[j + i * (NX)] - uMinusOne[j + i * (NX)] +
                                 tauSqr * (first + second);
#ifdef DEBUG
        temp_max = fmax(temp_max, fabs(uPlusOne[j + i * (NX)]));
#endif
      }
    }

    arg += tautwoPiF0;
    uPlusOne[SX + SY * (NX)] += tauSqr * expf(-sqr(arg * gamma)) * sinf(arg);

    float_t* t = uMinusOne;
    uMinusOne = U;
    U = uPlusOne;
    uPlusOne = t;
#ifdef DEBUG
    fprintf(max_u_file, "\nStep:%d\\%u Max: %f", n, NT, temp_max);
#endif
  }

#ifdef DEBUG
  FILE* output = fopen(OUTPUT_FILE, "wb");
  if (output == NULL) {
    printf("Error opening file");
    return (EXIT_FAILURE);
  }

  for (int i = 0; i < NY; ++i) {
    fwrite(U + i * NX, sizeof(float_t), NX, output);
  }

  fclose(output);
#endif
#ifdef DEBUG
  fclose(max_u_file);
#endif
  free(U);
  free(P);
  return 0;
}

__attribute__((always_inline)) inline float_t sqr(const float_t a) {
  return a * a;
}

float_t* create(int ny, int nx) {
  float_t* matrix = (float_t*)calloc(sizeof(float_t), nx * ny);
  return matrix;
}

float_t* calc_P(int ny, int nx) {
  float_t* matrix = (float_t*)malloc(sizeof(float_t) * nx * ny);
  const uint64_t border = nx / 2;
  for (int i = 0; i < ny; i++) {
    for (int j = 0; j < nx / 2; ++j) {
      matrix[j + (i * ny)] = 0.01;
      matrix[border + j + (i * ny)] = 0.04;
    }
  }
  return matrix;
}
