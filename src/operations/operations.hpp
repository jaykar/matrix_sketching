#ifndef __BOOST_WRAPPER_H__
#define __BOOST_WRAPPER_H__
#include <math.h>

template <typename SKMatrix>
SKMatrix gaussian_projection(SKMatrix& input, int s);

template <typename SKMatrix>
SKMatrix solve_x_sketch(SKMatrix& A, SKMatrix& B, int s);

template <typename SKMatrix>
SKMatrix count_sketch(const SKMatrix& A, int s);

#endif