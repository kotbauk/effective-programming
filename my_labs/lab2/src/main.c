#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <xmmintrin.h>
#include "inline_function.h"
#include "my_matrix.h"
#define DELTA_OUTPUT_FILE "delta_file.out"
#define RESULTS_FILE "results_matrix.bin"
#define BEGIN_X 0.0f
#define END_X 4.0f
#define BEGIN_Y 0.0f
#define END_Y 4.0f
#define VEC_TYPE_ALIGN 32
const uint64_t filler = 0xFFFFFFFFFFFFFFFF;
my_float_t unborder_filler;
typedef __m256d my_vec_t;
#define SCAL_IN_VEC sizeof(my_vec_t) / sizeof(my_float_t)
typedef struct modeling_area {
  my_float_t begin_x;
  my_float_t end_x;
  my_float_t begin_y;
  my_float_t end_y;
  my_float_t hx;
  my_float_t hy;
  size_t number_x_points;
  size_t number_y_points;
  size_t max_j;
  size_t max_i;
  size_t lenght_of_x_segment;
  size_t lenght_of_y_segment;

} modeling_area;
modeling_area init_modeling_area(const my_float_t begin_x,
                                 const my_float_t end_x,
                                 const my_float_t begin_y,
                                 const my_float_t end_y,
                                 const size_t number_x_points,
                                 const size_t number_y_points);
my_matrix_t* compute_po(my_matrix_t* po_matrix,
                        my_matrix_t* temp_matrix,
                        const modeling_area mod_area);
my_matrix_t* compute_results(my_matrix_t* results_matrix,
                             my_matrix_t* temp_matrix,
                             const my_matrix_t* po_matrix,
                             const size_t number_of_steps,
                             const modeling_area mod_area);
/* (a, b, c, d); (e, f, g, h) -> (d, e, f, g) */
my_vec_t get_left_vector(my_vec_t a, my_vec_t b) {
  return _mm256_permute4x64_pd(_mm256_blend_pd((a), (b), 0x07), 0b10010011);
}

/* (a, b, c, d); (e, f, g, h) -> (b, c, d, e) */
my_vec_t get_right_vector(my_vec_t a, my_vec_t b) {
  return _mm256_permute4x64_pd(_mm256_blend_pd((a), (b), 0x01), 0b111001);
}
int main(int argc, const char* argv[]) {
  const size_t number_x_points = atoi(argv[1]);
  const size_t number_y_points = atoi(argv[2]);
  const size_t steps_numbers = atoi(argv[3]);
  modeling_area mod_area = init_modeling_area(BEGIN_X, END_X, BEGIN_Y, END_Y,
                                              number_x_points, number_y_points);
  const size_t matrix_size =
      (mod_area.number_x_points) * (mod_area.number_y_points);

  my_float_t* matrixes_buffers =
      _mm_malloc(3 * (matrix_size) * sizeof(my_float_t), VEC_TYPE_ALIGN);
  memset(matrixes_buffers, 0, matrix_size);
  my_matrix_t* result_matrix = matrix_init(
      matrixes_buffers, mod_area.number_x_points, mod_area.number_y_points);
  my_matrix_t* temp_matrix =
      matrix_init(matrixes_buffers + matrix_size, mod_area.number_x_points,
                  mod_area.number_y_points);
  my_matrix_t* po_matrix =
      matrix_init(matrixes_buffers + 2 * matrix_size, mod_area.number_x_points,
                  mod_area.number_y_points);

  po_matrix = compute_po(po_matrix, temp_matrix, mod_area);
  FILE* results_file = fopen(RESULTS_FILE, "wb");
  print_matrix_bin_format(compute_results(result_matrix, temp_matrix, po_matrix,
                                          steps_numbers, mod_area),
                          results_file);

  fclose(results_file);
  _mm_free(matrixes_buffers);
  matrix_free(result_matrix);
  matrix_free(po_matrix);
  matrix_free(temp_matrix);
  return (EXIT_SUCCESS);
}

void swap_my_matrix_t(my_matrix_t** a, my_matrix_t** b) {
  my_matrix_t** tmp_ptr = a;
  *a = *b;
  *b = *tmp_ptr;
}
modeling_area init_modeling_area(const my_float_t begin_x,
                                 const my_float_t end_x,
                                 const my_float_t begin_y,
                                 const my_float_t end_y,
                                 const size_t number_x_points,
                                 const size_t number_y_points) {
  modeling_area mod_area = {.begin_x = begin_x,
                            .end_x = end_x,
                            .begin_y = begin_y,
                            .end_y = end_y,
                            .lenght_of_x_segment = end_x - begin_x,
                            .lenght_of_y_segment = end_y - begin_y,
                            .number_x_points = number_x_points,
                            .number_y_points = number_y_points,
                            .max_j = number_x_points - 2,
                            .max_i = number_y_points - 2,
                            .hx = (end_x - begin_x) / (number_x_points - 1),
                            .hy = (end_y - begin_y) / (number_y_points - 1)};
  return mod_area;
}
my_float_t compute_po_in_point(my_float_t xj,
                               my_float_t yi,
                               my_float_t xs1,
                               my_float_t xs2,
                               my_float_t ys1,
                               my_float_t ys2,
                               my_float_t sqr_r) {
  if (pow((xj - xs1), 2) + pow((yi - ys1), 2) < sqr_r) {
    return 0.1;
  } else if (pow((xj - xs2), 2) + pow((yi - ys2), 2) < sqr_r) {
    return -0.1;
  }
  return 0.0;
}

my_matrix_t* compute_po(my_matrix_t* po_matrix,
                        my_matrix_t* temp_matrix,
                        const struct modeling_area mod_area) {
  const my_float_t ys1 =
      mod_area.begin_y + ((2. / 3.) * mod_area.lenght_of_y_segment);
  const my_float_t ys2 = mod_area.begin_y + (mod_area.lenght_of_y_segment) / 3.;
  const my_float_t xs1 = mod_area.begin_x + (mod_area.lenght_of_x_segment) / 3.;
  const my_float_t xs2 =
      mod_area.begin_x + ((2. / 3.) * mod_area.lenght_of_x_segment);
  const my_float_t sqr_r =
      0.01 * fmin(mod_area.lenght_of_x_segment, mod_area.lenght_of_y_segment) *
      fmin(mod_area.lenght_of_x_segment, mod_area.lenght_of_y_segment);
  for (size_t i = 0; i < temp_matrix->y_dim_size; ++i) {
    for (size_t j = 0; j < temp_matrix->x_dim_size; ++j) {
      set_matrix_value(po_matrix, i, j,
                       compute_po_in_point(mod_area.begin_x + j * mod_area.hx,
                                           mod_area.begin_y + i * mod_area.hy,
                                           xs1, xs2, ys1, ys2, sqr_r));
    }
  }

  return po_matrix;
}

my_matrix_t* compute_results(my_matrix_t* results_matrix,
                             my_matrix_t* temp_matrix,
                             const my_matrix_t* po_matrix,
                             const size_t number_of_steps,
                             const modeling_area mod_area) {
#ifdef DEBUG_MODE
  FILE* delta_output_file = fopen(DELTA_OUTPUT_FILE, "w");
#endif
  const my_float_t sqr_hx = mod_area.hx * mod_area.hx;
  const my_float_t sqr_hy = mod_area.hy * mod_area.hy;
  const my_float_t sum_of_sqr_h = 1 / sqr_hx + 1 / sqr_hy;
  const my_vec_t vec_sqr_hx = _mm256_set1_pd(sqr_hx);
  const my_vec_t vec_sqr_hy = _mm256_set1_pd(sqr_hy);
  const my_vec_t vec_sum_of_sqr_h = _mm256_set1_pd(1 / sqr_hx + 1 / sqr_hy);
  const my_vec_t vec_first_multiplier = _mm256_set1_pd(0.2 / sum_of_sqr_h);
  const my_vec_t vec_second_multiplier =
      _mm256_set1_pd(2.5 / sqr_hx - 0.5 / sqr_hy);
  const my_vec_t vec_third_multiplier =
      _mm256_set1_pd(2.5 / sqr_hy - 0.5 / sqr_hx);
  const my_vec_t vec_fourth_multiplier = _mm256_set1_pd(0.25 * (sum_of_sqr_h));
  const my_vec_t vec_two = _mm256_set1_pd(2.f);
  const my_vec_t vec_quarter = _mm256_set1_pd(0.25f);
#ifdef DEBUG_MODE
  my_float_t max_delta = .0;
  my_float_t temp_delta = .0;
#endif
  const my_float_t border_filler = *(my_float_t*)&filler;
  for (size_t n = 0; n < number_of_steps; ++n) {
    for (size_t i = 1; i < mod_area.max_i + 1; ++i) {
      my_vec_t* temp_matrix_cur = (my_vec_t*)(get_matrix_buffer(temp_matrix) +
                                              i * mod_area.number_x_points);
      my_vec_t* results_matrix_up =
          (my_vec_t*)(get_matrix_buffer(results_matrix) +
                      (i - 1) * mod_area.number_x_points);
      my_vec_t* results_matrix_cur =
          (my_vec_t*)(get_matrix_buffer(results_matrix) +
                      i * mod_area.number_x_points);
      my_vec_t* results_matrix_low =
          (my_vec_t*)(get_matrix_buffer(results_matrix) +
                      (i + 1) * mod_area.number_x_points);

      my_vec_t vec_res_prev = _mm256_setzero_pd();
      my_vec_t vec_res_cur = results_matrix_cur[0];
      my_vec_t vec_res_next = results_matrix_cur[1];

      my_vec_t vec_res_up_prev = _mm256_setzero_pd();
      my_vec_t vec_res_up_cur = results_matrix_up[0];
      my_vec_t vec_res_up_next = results_matrix_up[1];

      my_vec_t vec_res_low_prev = _mm256_setzero_pd();
      my_vec_t vec_res_low_cur = results_matrix_low[0];
      my_vec_t vec_res_low_next = results_matrix_low[1];

      my_vec_t* po_low = (my_vec_t*)(get_matrix_buffer(po_matrix) +
                                     (i + 1) * mod_area.number_x_points);
      my_vec_t* po_cur = (my_vec_t*)(get_matrix_buffer(po_matrix) +
                                     i * mod_area.number_x_points);
      my_vec_t* po_up = (my_vec_t*)(get_matrix_buffer(po_matrix) +
                                    (i - 1) * mod_area.number_x_points);

      my_vec_t vec_po_prev = _mm256_setzero_pd();
      my_vec_t vec_po_cur = po_cur[0];
      my_vec_t vec_po_next = po_cur[1];

      my_vec_t vec_po_up = po_up[0];
      my_vec_t vec_po_low = po_low[0];

      for (size_t j = 0; j < mod_area.number_x_points; j += SCAL_IN_VEC) {
        my_float_t left_border_flag = 0.0;
        my_float_t right_border_flag = 0.0;
        /*if (j != 0) {*/
          /*left_border_flag = border_filler;*/
        /*}*/
        /*if (j < mod_area.number_x_points) {*/
          /*right_border_flag = border_filler;*/
        /*}*/
        /*__m256d border_mask = _mm256_set_pd(right_border_flag, border_filler,*/
                                            /*border_filler, left_border_flag);*/
        /*my_vec_t vec_po_left = _mm256_and_pd(*/
            /*get_left_vector(vec_po_prev, vec_po_cur), border_mask);*/
        /*my_vec_t vec_po_right = _mm256_and_pd(*/
            /*get_right_vector(vec_po_cur, vec_po_next), border_mask);*/

        /*my_vec_t vec_res_left = _mm256_and_pd(*/
            /*get_left_vector(vec_res_prev, vec_res_cur), border_mask);*/
        /*my_vec_t vec_res_right = _mm256_and_pd(*/
            /*get_right_vector(vec_res_cur, vec_res_next), border_mask);*/

        /*my_vec_t vec_res_up_left = _mm256_and_pd(*/
            /*get_left_vector(vec_res_up_prev, vec_res_up_cur), border_mask);*/
        /*my_vec_t vec_res_up_right = _mm256_and_pd(*/
            /*get_right_vector(vec_res_up_cur, vec_res_up_next), border_mask);*/

        /*my_vec_t vec_res_low_left = _mm256_and_pd(*/
            /*get_left_vector(vec_res_low_prev, vec_res_low_cur), border_mask);*/
        /*my_vec_t vec_res_low_right = _mm256_and_pd(*/
            /*get_right_vector(vec_res_low_cur, vec_res_low_next), border_mask);*/
        my_vec_t vec_po_left = 
            get_left_vector(vec_po_prev, vec_po_cur);
        my_vec_t vec_po_right = 
            get_right_vector(vec_po_cur, vec_po_next);

        my_vec_t vec_res_left = 
            get_left_vector(vec_res_prev, vec_res_cur);
        my_vec_t vec_res_right =
            get_right_vector(vec_res_cur, vec_res_next);

        my_vec_t vec_res_up_left =
            get_left_vector(vec_res_up_prev, vec_res_up_cur);
        my_vec_t vec_res_up_right = 
            get_right_vector(vec_res_up_cur, vec_res_up_next);

        my_vec_t vec_res_low_left = 
            get_left_vector(vec_res_low_prev, vec_res_low_cur);
        my_vec_t vec_res_low_right = 
            get_right_vector(vec_res_low_cur, vec_res_low_next);

        my_vec_t first_bracket = _mm256_mul_pd(
            vec_second_multiplier, _mm256_add_pd(vec_res_left, vec_res_right));

        my_vec_t second_bracket =
            _mm256_mul_pd(vec_third_multiplier,
                          _mm256_add_pd(vec_res_up_cur, vec_res_low_cur));

        my_vec_t third_bracket = _mm256_mul_pd(
            vec_fourth_multiplier,
            _mm256_add_pd(_mm256_add_pd(vec_res_up_left, vec_res_up_right),
                          _mm256_add_pd(vec_res_low_left, vec_res_low_right)));

        my_vec_t po_val = _mm256_fmadd_pd(
            vec_two, vec_po_cur,
            _mm256_mul_pd(
                vec_quarter,
                _mm256_add_pd(_mm256_add_pd(*po_low,*po_up),
                              _mm256_add_pd(vec_po_left, vec_po_right))));
        my_vec_t fourth_bracket =
            _mm256_add_pd(_mm256_add_pd(first_bracket, second_bracket),
                          _mm256_add_pd(third_bracket, po_val));
        *temp_matrix_cur = _mm256_mul_pd(vec_first_multiplier, fourth_bracket);

        ++temp_matrix_cur;
        ++results_matrix_up;
        ++results_matrix_cur;
        ++results_matrix_low;
        ++po_low;
        ++po_cur;
        ++po_up;

        vec_res_prev = vec_res_cur;
        vec_res_cur = vec_res_next;
        vec_res_next = results_matrix_cur[1];

        vec_res_up_prev = vec_res_up_cur;
        vec_res_up_cur = vec_res_up_next;
        vec_res_up_next = results_matrix_up[1];

        vec_res_low_prev = vec_res_low_cur;
        vec_res_low_cur = vec_res_low_next;
        vec_res_low_next = results_matrix_low[1];

        vec_po_prev = vec_po_cur;
        vec_po_cur = vec_po_next;
        vec_po_next = po_cur[1];

        vec_po_low = *po_low;
        vec_po_up = *po_up;
#ifdef DEBUG_MODE
        temp_delta = fabs(get_matrix_value(temp_matrix, i, j) -
                          get_matrix_value(results_matrix, i, j));
        max_delta = fmax(max_delta, temp_delta);
#endif
      }
    }
    my_matrix_t* temp_ptr = results_matrix;
    results_matrix = temp_matrix;
    temp_matrix = temp_ptr;
#ifdef DEBUG_MODE
    fprintf(delta_output_file, "Step:%zu, delta:%.7f\n", n, max_delta);
    max_delta = -0.1f;
#endif
  }
#ifdef DEBUG_MODE
  fclose(delta_output_file);
#endif
  return results_matrix;
}
