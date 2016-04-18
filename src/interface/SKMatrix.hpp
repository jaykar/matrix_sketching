#import <vector>
#include <random>
#include <cstdarg>
#include <time.h>

template <typename T, typename C>
class SKMatrix {
    public:
        virtual int size(void) const = 0;

        virtual int num_rows(void) const = 0;
        virtual int num_cols(void) const = 0;

        virtual std::vector<int>& dimensions(void) const = 0;
        virtual C& data() const = 0;

        virtual SKMatrix<T>& mult(SKMatrix<T>& rhs) const = 0;

        // Gaussian projection
        virtual T& rand_n(int row, int col, int mean, int std) const = 0;
        virtual T& elem_div(const T a) const = 0;

        virtual T& flip_signs(const std::vector<int> cols) = 0;

        std::vector<std::vector<int> >& bucket(const int num_buckets) const {
            if(num_buckets > this.cols){
                std::cout << "Number of buckets must be less than or equal to the number of columns";
                std::cout << '\n';
                throw;
            } else {
                std::vector<std::vector<int> > buckets(num_buckets);
                srand( time(NULL) );
                for(int i = 0; i < this.cols; i++){
                    int bucket = rand() % num_buckets;
                    buckets[bucket].push_back(i);
                }
            }

            return buckets;
        }

        virtual T& count_sketch(void) const = 0;

        /* Regression */
        virtual T& concat(const T& col) const  = 0;
        virtual T& solve_x(const T& A, const T& B) const = 0;

        /* TODO: K-SVD */
        virtual T& overridce_col(const int col, const SKMatrix& B) const = 0;
        virtual void qr_decompose(T& a, T& b) const = 0;
};

