#ifndef __OPERATIONS_H__
#define __OPERATIONS_H__

#include <math.h>
#include <vector>

namespace sk {
    namespace ops {
        /**
         * Funciton that returns a gaussian projection of given
         * @param input SKMatrix to take gaussian projection of
         * @param n a constant number of dimension to reduce input to
         * @return gaussian projection of a given graph in n dimensions
         */
        template <typename SKMatrix>
        SKMatrix gaussian_projection(const SKMatrix& input, const int n){
            int n_col = input.num_cols();
            SKMatrix a = SKMatrix();
            a.rand_n(n_col, n);
            a.elem_div(sqrt(n*1.0));
            return input.mult(a);
        }

        /**
         * Takes a linear regression on Ax=B and returns the solution
         * using gaussian projection to sketch
         * @param A matrix A in Ax=B
         * @param B matrix B in Ax=B
         * @param n a constant number of dimension to reduce A to
         * @return approximated matrix from linear regression
         */
        template <typename SKMatrix>
        SKMatrix& lin_regress(const SKMatrix& A, const SKMatrix& B, const int n){
            int n_col = A.num_cols();

            auto concat_mat = A.concat(B);
            concat_mat.transpose();

            auto sketch = gaussian_projection<SKMatrix>(concat_mat, n);
            sketch.transpose();

            auto A_sketch = sketch.get_cols(0, n_col -1);
            auto B_sketch = sketch.get_col(n_col);
            auto X_tilde  = A_sketch.solve_x(B_sketch);

            return X_tilde;
        }

        /**
         * Count sketch algorithm in map-reduce
         * @param A matrix to perform count sketch on
         * @param n number of buckets/hash sets
         * @return count sketch matrix
         */
        template <typename SKMatrix>
        SKMatrix count_sketch(const SKMatrix& A, const int n) {
            std::vector<std::vector<int> > buckets = A.bucket(n);
            std::vector<bool> flipped = A.flip_signs();

            SKMatrix sum(A.num_rows(), 0);

            for(auto hash_set : buckets){
                SKMatrix set(A.num_rows(), 1);
                for(int index : hash_set) {
                    SKMatrix col = A.get_col(index);
                    if(flipped[index]) {
                        col.data() *= -1;
                    }
                    set.data() += col.data();
                }
                sum = sum.concat(set);
            }
            return sum;
        }

        template <typename SK>
        std::vector<SK> k_svd(const SK&A, int k, int s){
            auto ans = std::vector<SK>(3);
            //auto sketch = count_sketch(A, s);
            auto sketch = gaussian_projection(A, s);
            std::cout << sketch.size() << std::endl;
            SK Q;
            SK R;
            sketch.qr_decompose(Q,R);

            std::cout << Q.size() << std::endl;
            std::cout << R.size() << std::endl;
            auto Q_temp = Q;
            Q_temp.transpose();
            Q_temp = Q_temp.mult(A);
            auto svd_vec = Q_temp.svds(k);
            auto U_tilde = Q.mult(svd_vec[0]);
            ans[0] = U_tilde;
            ans[1] = svd_vec[1];
            ans[2] = svd_vec[2];
            return ans;
        }
    }
}
#endif