#include <iostream>
#include "../sk_boost.hpp"
#include "../../interface/SKMatrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>

namespace bnu = boost::numeric::ublas;

template <typename SKMatrix>
void test(SKMatrix value){
}

int main(){
  bnu::matrix<float> m1(5,5, 1.0);
  sk_boost<float> mat1(m1);
  test<sk_boost<float>>(mat1);
  // bnu::matrix<float> m2(5,5, 2.0);
  // sk_boost<float> mat2(m2);

  // std::cout << "mat 1" << std::endl;
  // std::cout << mat1.data() << std::endl;

  // std::cout << "mat 2" << std::endl;
  // std::cout << mat2.data() << std::endl;

  // std::cout << "mat 1 subtract mat 2" << std::endl;
  // std::cout << mat1.subtract(mat2).data() << std::endl;

  // std::cout << "mat 1 concat mat 2" << std::endl;
  // std::cout << mat1.concat(mat2).data() << std::endl;

  // std::cout << "mat 1 mult mat 2" << std::endl;
  // std::cout << mat1.mult(mat2).data() << std::endl;

  // auto buckets = mat1.bucket(3);
  // std::cout << buckets[0].size() << std::endl;
  // std::cout << buckets[1].size() << std::endl;
  // std::cout << buckets[2].size() << std::endl;

  // std::cout << "rand_n(2, 2)" << std::endl;
  // std::cout << mat1.rand_n(2, 2).data() << std::endl;

  // std::cout << "mat1 elem_div 1.5" << std::endl;
  // std::cout << mat1.elem_div(1.5).data() << std::endl;

  // std::cout << "mat1 get_col 1" << std::endl;
  // std::cout << mat1.get_col(2).data() << std::endl;

  // std::cout << "mat1 get_col 2, 4" << std::endl;
  // std::cout << mat1.get_cols(2,4).data() << std::endl;

  // std::cout << "mat1 accumulate" << std::endl;
  // std::cout << mat1.accumulate() << std::endl;
  // bnu::matrix<double> B(3, 3);
  // B(0,0) = 12;
  // B(0,1) = -51;
  // B(0,2) = 4;

  // B(1,0) = 6;
  // B(1,1) = 167;
  // B(1,2) = -68;

  // B(2,0) = -4;
  // B(2,1) = 24;
  // B(2,2) = -41;
  // std::cout << B << std::endl;

  // bnu::matrix<double> q(3, 3);
  // bnu::matrix<double> r(3, 3);

  // sk_boost<double> matB(B);
  // sk_boost<double> Q(q);
  // sk_boost<double> R(r);
  // std::cout << std::endl;
  // matB.qr_decompose(Q, R);

  // std::cout << Q.data() << std::endl;
  // std::cout << R.data() << std::endl;

  // bnu::matrix<float> D(2,1);
  // D(0,0) = 5;
  // D(1,0) = 7;
  // sk_boost<float> matD(D);

  // bnu::matrix<float> E(2,2);
  // E(0,0) = 3;
  // E(0,1) = 4;

  // E(1,0) = 2;
  // E(1,1) = -1;
  // sk_boost<float> matE(E);

  // std::cout << "Ax = B" << std::endl;

  // std::cout << matE.solve_x(matD).data() << std::endl;

  return 0;
}