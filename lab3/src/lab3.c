
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const char OUTPUT_FILE[] = "results.bin";

#define TYPE_ALIGNMENT 32
#define SIZE_OF_VECTOR sizeof(__m256) / sizeof(float)
#define MY_MASK
typedef float float_t;
float_t tau = 0.01;
const float_t XA = 0.f;
const float_t XB = 4.0f;
const float_t YA = 0.0f;
const float_t YB = 4.0f;
const float_t f0 = 1.0f;
const float_t t0 = 1.5f;
const float_t gamma = 0.25f;
const uint32_t SX = 1;
const float_t twoPiF0 = M_PI * 2 * f0;
const float_t t0twoPiF0 = t0 * twoPiF0;
const float_t t0twoPi = t0 * 2 * M_PI;
const uint32_t filler = 0xFFFFFFFF;
float_t* create(int ny, int nx);
float_t* calc_P(int ny, int nx);
inline void swap(float_t* a, float_t* b) {
  float_t tmp = *a;
  *a = *b;
  *b = tmp;
}
float_t sqr(const float_t a);

/*
 *Функция "сдвигающая" вектор влево на один элемент в матрице
 */
/* (a, b, c, d, e, f, g, h)  (a2, b2, c2, d2, e2, f2, g2, h2)->
 * (h, a2, b2, c2, d2, e2, f2, g2) */
inline __m256 get_left_vector(__m256 a, __m256 b) {
  float temp_vec[8];
  _mm256_store_ps(temp_vec, _mm256_permute_ps(_mm256_blend_ps(a, b, 0b01111111),
                                              0b10010011));
  swap(temp_vec, temp_vec + 4);
  return _mm256_load_ps(temp_vec);
}
/*
 *Функция "сдвигающая" вектор вправо на один элемент в матрице
 */
/* (a, b, c, d, e, f, g, h)  (a2, b2, c2, d2, e2, f2, g2, h2)->
 * (b, c, d, e, f, g, a2) */
inline __m256 get_right_vector(__m256 a, __m256 b) {
  float temp_vec[8];
  _mm256_store_ps(temp_vec, _mm256_permute_ps(_mm256_blend_ps(a, b, 0b00000001),
                                              0b00111001));
  swap(temp_vec + 3, temp_vec + 7);
  return _mm256_load_ps(temp_vec);
}

int main(int argc, const char* argv[]) {
  const uint32_t NX = atoi(argv[1]);
  const uint32_t NY = atoi(argv[2]);
  const uint32_t NT = atoi(argv[3]);
  const uint32_t SY = NY / 2;
  const uint64_t y_bord = NY - 1;
  const float_t HXdeg2mul2 = 1 / (2 * sqr((XB - XA) / (NX - 1)));
  const float_t HYdeg2mul2 = 1 / (2 * sqr((YB - YA) / (NY - 1)));
  const __m256 vec_HXdeg2mul2 = _mm256_set1_ps(HXdeg2mul2);
  const __m256 vec_HYdeg2mul2 = _mm256_set1_ps(HYdeg2mul2);
  const __m256 vec_two = _mm256_set1_ps(2.0f);

#ifdef DEBUG
  FILE* max_u_file = fopen("max_u_file", "wb");
#endif

  if (NX > 1000 || NY > 1000) {
    tau = 0.001;
  }

  const float_t tauSqr = sqr(tau);
  const __m256 vec_tauSqr = _mm256_set1_ps(tauSqr);

  float_t* U = create(NY, NX);
  float_t* uMinusOne = create(NY, NX);
  float_t* uPlusOne = create(NY, NX);
  float_t* P = calc_P(NY, NX);

#ifdef DEBUG
  float_t unpacked_max_vec[SIZE_OF_VECTOR];
  float_t temp_max = -0.1;
  __m256 vec_temp_max = _mm256_setzero_ps();
  const __m256 abs_mask =
      _mm256_castsi256_ps(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
#endif
  const float_t tautwoPiF0 = tau * t0twoPiF0;
  float_t arg = -tautwoPiF0 - t0twoPi;
  const float_t unborder_filler = *(float_t*)&filler;
  for (size_t n = 0; n < NT; n++) {
    for (size_t i = 1; i < y_bord; i++) {
      size_t cur_line = i * NX;

      __m256* vec_U_cur = (__m256*)(U + cur_line);
      __m256* vec_uPlusOne_cur = (__m256*)(uPlusOne + cur_line);
      __m256* vec_uMinusOne_cur = (__m256*)(uMinusOne + cur_line);
      __m256* vec_U_line_upper = (__m256*)(U + cur_line - NX);
      __m256* vec_U_line_lower = (__m256*)(U + cur_line + NX);

      __m256 vec_U_prev = _mm256_setzero_ps();
      __m256 vec_U_mid = vec_U_cur[0];
      __m256 vec_U_next = vec_U_cur[1];

      __m256* vec_P_upper = (__m256*)(P + cur_line - NX);
      __m256* vec_P_lower = (__m256*)(P + cur_line);
      __m256 vec_P_upper_prev = _mm256_setzero_ps();
      __m256 vec_P_lower_prev = _mm256_setzero_ps();

      for (size_t j = 0; j < NX; j += SIZE_OF_VECTOR) {
        __m256 P_01 = *vec_P_upper;                             // P[i-1][j]
        __m256 P_11 = *vec_P_lower;                             // P[i][j]
        __m256 P_00 = get_left_vector(vec_P_upper_prev, P_01);  // P[i-1][j-1]
        __m256 P_10 = get_left_vector(vec_P_lower_prev, P_11);  // P[i][j]
        float_t left_border_elem = 0.0f;
        float_t right_border_elem = 0.0f;

        if (j != 0) {
          left_border_elem = unborder_filler;
        }

        if (j + SIZE_OF_VECTOR < NX) {
          right_border_elem = unborder_filler;
        }
        __m256 border_mask =
            _mm256_set_ps(left_border_elem, unborder_filler, unborder_filler,
                          unborder_filler, unborder_filler, unborder_filler,
                          unborder_filler, right_border_elem);

        __m256 vec_U_left_val =
            _mm256_and_ps(get_left_vector(vec_U_prev, vec_U_mid), border_mask);
        __m256 vec_U_right_val =
            _mm256_and_ps(get_right_vector(vec_U_mid, vec_U_next), border_mask);

        __m256 vec_U1 = _mm256_sub_ps(*vec_U_line_upper, vec_U_mid);
        __m256 vec_P1 = _mm256_add_ps(P_00, P_01);

        __m256 vec_U2 = _mm256_sub_ps(*vec_U_line_lower, vec_U_mid);
        __m256 vec_P2 = _mm256_add_ps(P_10, P_11);

        __m256 vec_first =
            _mm256_mul_ps((_mm256_add_ps(_mm256_mul_ps(vec_U1, vec_P1),
                                         _mm256_mul_ps(vec_U2, vec_P2))),
                          vec_HXdeg2mul2);

        __m256 vec_U3 = _mm256_sub_ps(vec_U_right_val, vec_U_mid);
        __m256 vec_P3 = _mm256_add_ps(P_01, P_11);

        __m256 vec_U4 = _mm256_sub_ps(vec_U_left_val, vec_U_mid);
        __m256 vec_P4 = _mm256_add_ps(P_00, P_10);

        __m256 vec_second =
            _mm256_mul_ps((_mm256_add_ps(_mm256_mul_ps(vec_U3, vec_P3),
                                         _mm256_mul_ps(vec_U4, vec_P4))),
                          vec_HYdeg2mul2);

        *vec_uPlusOne_cur = _mm256_fmadd_ps(
            vec_tauSqr, _mm256_add_ps(vec_first, vec_second),
            _mm256_fmsub_ps(vec_two, vec_U_mid, *vec_uMinusOne_cur));

#ifdef DEBUG
        vec_temp_max = _mm256_max_ps(
            vec_temp_max, _mm256_and_ps(*vec_uPlusOne_cur, abs_mask));
#endif
        ++vec_U_cur;
        ++vec_uPlusOne_cur;
        ++vec_uMinusOne_cur;
        ++vec_U_line_lower;
        ++vec_U_line_upper;

        vec_U_prev = vec_U_mid;
        vec_U_mid = vec_U_next;
        vec_U_next = vec_U_cur[1];

        ++vec_P_upper;
        ++vec_P_lower;
        vec_P_upper_prev = P_01;
        vec_P_lower_prev = P_11;
      }
    }

    arg += tautwoPiF0;
    uPlusOne[SX + SY * (NX)] += tauSqr * expf(-sqr(arg * gamma)) * sinf(arg);

    float_t* t = uMinusOne;
    uMinusOne = U;
    U = uPlusOne;
    uPlusOne = t;
#ifdef DEBUG
    _mm256_store_ps(unpacked_max_vec, vec_temp_max);
    for (size_t i = 0; i < SIZE_OF_VECTOR; ++i) {
      if (temp_max < unpacked_max_vec[i]) {
        temp_max = unpacked_max_vec[i];
      }
    }
    fprintf(max_u_file, "\nStep:%zu\\%u Max: %f", n, NT, temp_max);
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

inline float_t sqr(const float_t a) {
  return a * a;
}

float_t* create(int ny, int nx) {
  size_t my_size = ny * nx * sizeof(float_t);
  float_t* matrix = (float_t*)_mm_malloc(my_size, sizeof(__m256));
  memset(matrix, 0, my_size);
  return matrix;
}

float_t* calc_P(int ny, int nx) {
  size_t my_size = ny * nx * sizeof(float_t);
  float_t* matrix = (float_t*)_mm_malloc(my_size, sizeof(__m256));
  const uint64_t border = nx / 2;
  for (int i = 0; i < ny; i++) {
    for (int j = 0; j < nx / 2; ++j) {
      matrix[j + (i * ny)] = 0.01;
      matrix[border + j + (i * ny)] = 0.04;
    }
  }
  return matrix;
}
