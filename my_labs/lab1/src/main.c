#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "inline_function.h"
#include "my_matrix.h"
#define DELTA_OUTPUT_FILE "delta_file.out"
#define RESULTS_FILE "results_matrix.bin"
#define BEGIN_X 0.0f
#define END_X 4.0f
#define BEGIN_Y 0.0f
#define END_Y 4.0f
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
int main(int argc, const char* argv[]) {
  const size_t number_x_points = atoi(argv[1]);
  const size_t number_y_points = atoi(argv[2]);
  const size_t steps_numbers = atoi(argv[3]);
  modeling_area mod_area = init_modeling_area(BEGIN_X, END_X, BEGIN_Y, END_Y,
                                              number_x_points, number_y_points);
  const size_t matrix_size =
      (mod_area.number_x_points) * (mod_area.number_y_points);

  my_float_t* matrixes_buffers = calloc(4 * (matrix_size), sizeof(my_float_t));
  my_matrix_t* result_matrix = matrix_init(
      matrixes_buffers, mod_area.number_x_points, mod_area.number_y_points);
  my_matrix_t* po_matrix =
      matrix_init(matrixes_buffers + matrix_size, mod_area.number_x_points,
                  mod_area.number_y_points);
  my_matrix_t* temp_matrix =
      matrix_init(matrixes_buffers + 2 * matrix_size, mod_area.number_x_points,
                  mod_area.number_y_points);

  po_matrix = compute_po(po_matrix, temp_matrix, mod_area);
  FILE* results_file = fopen(RESULTS_FILE, "wb");
  print_matrix_bin_format(compute_results(result_matrix, temp_matrix, po_matrix,
                                          steps_numbers, mod_area),
                          results_file);

  fclose(results_file);
  free(matrixes_buffers);
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
      set_matrix_value(temp_matrix, i, j,
                       compute_po_in_point(mod_area.begin_x + j * mod_area.hx,
                                           mod_area.begin_y + i * mod_area.hy,
                                           xs1, xs2, ys1, ys2, sqr_r));
    }
  }

  for (size_t i = 0; i < po_matrix->y_dim_size; ++i) {
    for (size_t j = 0; j < po_matrix->x_dim_size; ++j) {
      set_matrix_value(po_matrix, i, j,
                       2. * get_matrix_value(temp_matrix, i, j) +
                           0.25 * (get_matrix_value(temp_matrix, i - 1, j) +
                                   get_matrix_value(temp_matrix, i + 1, j) +
                                   get_matrix_value(temp_matrix, i, j - 1) +
                                   get_matrix_value(temp_matrix, i, j + 1)));
    }
  }

  return po_matrix;
}
inline my_float_t __attribute__((__always_inline__))
compute_function_in_point(const my_matrix_t* results_matrix,
                          const my_matrix_t* po_matrix,
                          const uint32_t i,
                          const uint32_t j,
                          const my_float_t first_multiplier,
                          const my_float_t second_multiplier,
                          const my_float_t third_multiplier,
                          const my_float_t fourth_multiplier) {
  return first_multiplier *
         (second_multiplier * (get_matrix_value(results_matrix, i, j - 1) +
                               get_matrix_value(results_matrix, i, j + 1)) +
          third_multiplier * (get_matrix_value(results_matrix, i - 1, j) +
                              get_matrix_value(results_matrix, i + 1, j)) +
          fourth_multiplier * (get_matrix_value(results_matrix, i - 1, j - 1) +
                               get_matrix_value(results_matrix, i - 1, j + 1) +
                               get_matrix_value(results_matrix, i + 1, j - 1) +
                               get_matrix_value(results_matrix, i + 1, j + 1)) +
          get_matrix_value(po_matrix, i, j));
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
  const my_float_t first_multiplier = 0.2 / sum_of_sqr_h;
  const my_float_t second_multiplier = (2.5 / sqr_hx - 0.5 / sqr_hy);
  const my_float_t third_multiplier = (2.5 / sqr_hy - 0.5 / sqr_hx);
  const my_float_t fourth_multiplier = 0.25 * (sum_of_sqr_h);
#ifdef DEBUG_MODE
  my_float_t max_delta = .0;
  my_float_t temp_delta = .0;
#endif
  for (size_t n = 1; n < number_of_steps; ++n) {
    for (size_t i = 1; i < mod_area.max_i; ++i) {
      for (size_t j = 1; j < mod_area.max_j; ++j) {
        set_matrix_value(
            temp_matrix, i, j,
            compute_function_in_point(results_matrix, po_matrix, i, j,
                                      first_multiplier, second_multiplier,
                                      third_multiplier, fourth_multiplier));
#ifdef DEBUG_MODE
        temp_delta = fabs(get_matrix_value(temp_matrix, i, j) -
                          get_matrix_value(results_matrix, i, j));
        max_delta = fmax(max_delta, temp_delta);
#endif
      }
    }
    swap_my_matrix_t(&results_matrix, &temp_matrix);
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
