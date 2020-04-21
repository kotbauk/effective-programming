#ifndef MY_MATRIX_H
#define MY_MATRIX_H
#include <stdlib.h>
#include <stdio.h>

typedef double my_float_t;
typedef struct my_matrix_t {
  my_float_t* matrix_buf;
  size_t x_dim_size;
  size_t y_dim_size;
} my_matrix_t;

my_matrix_t* matrix_init(my_float_t* buf,
                         const size_t x_dim_size,
                         const size_t y_dim_size);
void matrix_free(my_matrix_t* matrix);
my_matrix_t* set_matrix_buf(my_matrix_t* matrix, my_float_t* new_buf);
void set_matrix_value(my_matrix_t* matrix,
                      const size_t i,
                      const size_t j,
                      const my_float_t val);
my_float_t get_matrix_value(const my_matrix_t* matrix,
                            const size_t i,
                            const size_t j);
my_float_t* get_matrix_buffer(const my_matrix_t* matrix);


void print_matrix_bin_format(const my_matrix_t* matrix, FILE* output_file);
void print_matrix_txt_format(const my_matrix_t* matrix, FILE* output_file);
#endif
