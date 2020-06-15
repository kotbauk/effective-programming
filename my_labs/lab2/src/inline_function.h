#ifndef INLINE_FUNCTION_H
#define INLINE_FUNCTION_H
#include <xmmintrin.h>
#include "my_matrix.h"
static inline my_float_t* __attribute__((__always_inline__))
at_matrix_value(const my_matrix_t* matrix, const size_t i, const size_t j) {
  return matrix->matrix_buf + (j + i * matrix->x_dim_size);
}

inline void __attribute__((__always_inline__))
set_matrix_value(my_matrix_t* matrix,
                 const size_t i,
                 const size_t j,
                 const my_float_t val) {
  *at_matrix_value(matrix, i, j) = val;
}
inline my_float_t __attribute__((__always_inline__))
get_matrix_value(const my_matrix_t* matrix, const size_t i, const size_t j) {
  return *at_matrix_value(matrix, i, j);
}
#endif
