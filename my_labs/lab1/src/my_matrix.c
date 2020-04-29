#include "my_matrix.h"

#include <stdio.h>

my_matrix_t* matrix_init(my_float_t* buf,
                         const size_t x_dim_size,
                         const size_t y_dim_size) {
  my_matrix_t* matrix = malloc(sizeof(my_matrix_t));
  if (matrix == NULL) {
    fprintf(stderr, "Can't allocate memory for a matrix");
    return matrix;
  }
  matrix->matrix_buf = buf;
  matrix->x_dim_size = x_dim_size;
  matrix->y_dim_size = y_dim_size;
  return matrix;
}

void matrix_free(my_matrix_t* matrix) {
  free(matrix);
}

my_matrix_t* set_matrix_buf(my_matrix_t* matrix, my_float_t* new_buf) {
  matrix->matrix_buf = new_buf;
  return matrix;
}

my_float_t* get_matrix_buffer(const my_matrix_t* matrix) {
	return matrix->matrix_buf;
}


void print_matrix_bin_format(const my_matrix_t* matrix, FILE* output_file) {
  for (size_t i = 0; i < matrix->y_dim_size; ++i) {
    fwrite(get_matrix_buffer(matrix) + i * (matrix->x_dim_size),
           sizeof(my_float_t), matrix->x_dim_size, output_file);
  }
}
/*[>void print_matrix(const my_matrix_t* matrix, FILE* output_file) {<]*/
  /*[>for (size_t i = 0; i < matrix->y_dim_size; ++i) {<]*/
    /*[>for (size_t j = 0; j < matrix->x_dim_size; ++j) {<]*/
      /*[>fprintf(output_file, "%f|", get_matrix_value(matrix, i, j));<]*/
    /*[>}<]*/
    /*[>fprintf(output_file, "\n");<]*/
  /*[>}<]*/
/*}*/
