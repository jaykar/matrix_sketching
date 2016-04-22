#include "../sk_boost.hpp"
#include <boost/numeric/ublas/matrix.hpp>

namespace bnu = boost::numeric::ublas;

int main(){
  bnu::matrix<double> m(1,2);
  sk_boost<double> bw(m);
  return 0;
}