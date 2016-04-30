#include "SKMatrix.hpp"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <helper_cuda.h>
#include <sstream>
#include <chrono>

namespace sketchy {

    class cuda : SKMatrix<cuda, float *>  {
        private:
            int rows;
            int cols;
            float *gpu_data;
            mutable cublasStatus_t stat;
            mutable cudaError_t err;
            cublasHandle_t handle;

            /* copies matrix from GPU memory to CPU memory */
            void getMatrix() const {
                stat = cublasGetMatrix(rows, cols, sizeof(float), (void *) gpu_data, rows, (void *) matrix_data, rows);
                std::cout << "getting" << std::endl;
                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "upload failed: " << _cudaGetErrorEnum(stat) << std::endl;
                    exit(1);
                }
            }

            /* copies matrix from CPU memory to GPU memory */
            void setMatrix() {
                stat = cublasSetMatrix(rows, cols, sizeof(float), matrix_data, rows, gpu_data, rows);
                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "data download failed: " << _cudaGetErrorEnum(stat) << stat << std::endl;
                    cudaFree(gpu_data);
                    cublasDestroy(handle);
                    exit(1);
                }
            }

        public:
            cuda() {
                matrix_data = NULL;
                gpu_data = NULL;
                rows = 0;
                cols = 0;
            }

            cuda(const int row, const int col) {
                std::cout << "parametric constructor" << std::endl;

                matrix_data = (float *) malloc(row * col * sizeof(float));
                memset(matrix_data, 0, row * col * sizeof(float));
                
                stat = cublasCreate(&handle);
                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "CUBLAS initialization failed" << std::endl;
                    exit(1);
                }

                err = cudaMalloc((void **) &gpu_data, row * col * sizeof(float));

                if(err != cudaSuccess) {
                    std::cout << "gpu memory allocation failed" << std::endl;
                    cublasDestroy(handle);
                    exit(1);
                }

                rows = row;
                cols = col;

                setMatrix();

            }

            void identity(const int n) {
                *this = cuda(n, n);
                int i;
                for(i = 0; i < n; i++)
                    matrix_data[i*cols + i] = 1.0;
            }

            void set(int r, int c, float a) {

                this->matrix_data[r*this->rows + c] = a;

            }

            ~cuda() {
                
                cudaFree(gpu_data);
                stat = cublasDestroy(handle);
                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "CUBLAS destruction failed" << std::endl;
                    exit(1);
                }
                free(matrix_data);
            }

            /*
            cuda(const cuda& c) {
                std::cout << "copy constructor" << std::endl;
                matrix_data = (float *) malloc(c.size() * sizeof(float));

                err = cudaMalloc((void **) &gpu_data, c.size() * sizeof(float));

                if(err != cudaSuccess) {
                    std::cout << "gpu memory allocation failed " << err << std::endl;
                    cublasDestroy(handle);
                    exit(1);
                }

                stat = cublasCreate(&handle);
                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "CUBLAS initialization failed" << std::endl;
                    exit(1);
                }
                c.getMatrix();

                stat = cublasScopy(c.handle, c.size(), c.matrix_data, 1, matrix_data, 1);

                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "cublasScopy failed: " << _cudaGetErrorEnum(stat) << std::endl;
                    cudaFree(gpu_data);
                    cublasDestroy(handle);
                    free(matrix_data);
                    exit(1);
                }

                rows = c.rows;
                cols = c.cols;

                setMatrix();

            }

            cuda& operator=(const cuda& rhs) {
                std::cout << "assignment" <<std::endl;
                if(this->matrix_data)
                    free(this->matrix_data);
                if(this->gpu_data)
                    cudaFree(this->gpu_data);
                stat = cublasCreate(&handle);
                if(stat != CUBLAS_STATUS_SUCCESS) {
                    std::cout << "CUBLAS initialization failed" << std::endl;
                    exit(1);
                }
                this->rows = rhs.rows;
                this->cols = rhs.cols;
                this->matrix_data = (float *) malloc(rhs.size() * sizeof(float));
                cublasScopy(rhs.handle, rhs.size(), rhs.matrix_data, 1, matrix_data, 1);

                err = cudaMalloc((void **) &gpu_data, rhs.rows * rhs.cols * sizeof(float));

                if(err != cudaSuccess) {
                    std::cout << "gpu memory allocation failed" << std::endl;
                    cublasDestroy(handle);
                    exit(1);
                }

                setMatrix();
                return *this;
            }
            */

            /*
            void clear() {
                if(this->matrix_data)
                    memset(this->matrix_data, 0, this->rows * this->cols * sizeof(float));
            }
            */

            int size() const {
                return this->rows * this->cols;
            }

            int num_rows(void) const {
                return this->rows;
            }

            int num_cols(void) const {
                return this->cols;
            }

            float *data() const {
                getMatrix();
                return this->matrix_data;
            }

            cuda rand_n(const int row, const int col) {

                if(this->matrix_data)
                    free(this->matrix_data);
                if(this->gpu_data)
                    cudaFree(gpu_data);

                this->rows = row;
                this->cols = col;

                this->matrix_data = (float *) malloc(row * col * sizeof(float));
                err = cudaMalloc((void **) &this->gpu_data, row * col * sizeof(float));
                if(err != cudaSuccess) {
                    std::cout << "gpu memory allocation failed" << std::endl;
                    cublasDestroy(handle);
                    exit(1);
                }
                stat = cublasCreate(&handle);
                int i;
                unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                std::mt19937 gen(seed);
                std::normal_distribution<float> n;

                for(i = 0; i < this->size(); i++) {
                    this->matrix_data[i] = n(gen);
                }

                setMatrix();

                return *this;
            }

            /*
            cuda mult(const cuda& rhs) const {
                if(this->cols != rhs.rows)
                    throw std::invalid_argument("mismatched dimensions");
                cuda product(this->rows, rhs.cols);
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->rows,
                        rhs.cols, this->cols, 1, this->matrix_data, rhs.rows,
                        rhs.matrix_data, rhs.cols, 0, product.matrix_data, rhs.cols);
                //std::cout << "p: " << product << std::endl;
                return product;
            }

            cuda elem_div(const double a) const {
                if(a == 0.0)
                    throw std::invalid_argument("can't divide by zero");
                int i;
                cuda temp;
                for(i = 0; i < this->size(); i++)
                    temp.matrix_data[i] = this->matrix_data[i] / a;
                return temp;
            }

            //TODO
            cuda concat(const cuda& col) const {
                if(this->rows != col.rows)
                    throw std::invalid_argument("mismatched dimensions in concat");

                cuda ret(this->rows, this->cols + col.cols);
                float *current = ret.matrix_data;
                int i;
                for(i = 0; i < this->rows; i++) {
                    cblas_scopy(this->cols, current, 1, this->matrix_data, 1);
                    current += this->cols;
                    cblas_scopy(col.cols, current, 1, col.matrix_data, 1);
                    current += col.cols;
                }
                return ret;
            }

            cuda solve_x(const cuda& B) const {
                if(this->rows != this->cols)
                    throw std::invalid_argument("A must be square in Ax = B");
                if(this->cols != B.rows)
                    throw std::invalid_argument("B.rows must equal A.cols in Ax = B");
                if(B.cols != 1)
                    throw std::invalid_argument("B must be nx1 vector in Ax = B");
                cuda X(B.rows, 1);
                MKL_INT lda = this->rows; // rows in A
                MKL_INT n = this->cols; // cols in A
                MKL_INT nrhs = 1; // cols in B
                MKL_INT ldb = nrhs;
                MKL_INT ipiv[n];

                cuda copy(B);
                LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, nrhs, this->matrix_data, lda, ipiv, copy.matrix_data, ldb);

                return copy;
            }

            //TODO
            cuda get_col(const int n) const {
                return get_cols(n, n);
            }

            //TODO
            cuda get_cols(const int start, const int end) const {
                return cuda();
            }

            void transpose() {
                mkl_simatcopy('r', 't', this->rows, this->cols, 1.0, this->matrix_data,
                        this->cols, this->rows);
                int i = this->rows;
                this->rows = this->cols;
                this->cols = i;
            }

            cuda subtract(const cuda& rhs) const {
                if(this->dimensions() != rhs.dimensions())
                    throw std::invalid_argument("mismatched dimensions");
                cuda temp(*this);
                cblas_saxpy(temp.size(), -1, rhs.matrix_data, 1, temp.matrix_data, 1);
                return temp;
            }

            //TODO
            float accumulate() const {
                return 0.0;
            }

            //TODO
            void qr_decompose(cuda& Q, cuda& R) const {
            }

            //TODO
            void svd(cuda& U, cuda& S, cuda& V, const int k) const {
                // return std::vector<cuda>();
            }

            std::vector<int> dimensions() const {
                std::vector<int> ret = {this->rows, this->cols};
                return ret;
            }

            //gpu parallelizable
            cuda add(const cuda& rhs) const {
                if(this->dimensions() != rhs.dimensions())
                    throw std::invalid_argument("mismatched dimensions");
                cuda temp(*this);
                std::cout << "temp: " << temp << std::endl;
                cblas_saxpy(temp.size(), 1, rhs.matrix_data, 1, temp.matrix_data, 1);
                return temp;
            }

            cuda operator+(const cuda& rhs) const {
                return this->add(rhs);
            }

            cuda& operator+=(const cuda& rhs) {
                *this = *this + rhs;
                return *this;
            }

            cuda operator-(const cuda& rhs) const {
                return this->subtract(rhs);
            }

            cuda& operator-=(const cuda& rhs) {
                *this = *this - rhs;
                return *this;
            }

            cuda operator*(const cuda& rhs) const {
                return this->mult(rhs);
            }

            cuda& operator*=(const cuda& rhs) {
                *this = this->mult(rhs);
                return *this;
            }

            cuda operator/(const float a) {
                return this->elem_div(a);
            }

            cuda& operator/=(const float a) {
                *this = this->elem_div(a);
                return *this;
            }
            */

            /*
            // Count Sketch

            // Regression
             */


            /*
            // TODO: K-SVD
            cuda& override_col(const int col, const cuda& B) const;

            friend bool operator==(const cuda& lhs, const cuda& rhs);
             */
            friend std::ostream& operator<<(std::ostream&os, const cuda& im);
    };

    /*
    bool operator==(const cuda& lhs, const cuda& rhs) {
        if(rhs.rows != lhs.rows || rhs.cols != lhs.cols)
            return false;
        int i;
        for(i = 0; i < lhs.size(); i++)
            if(lhs.matrix_data[i] != rhs.matrix_data[i])
                return false;
        return true;
    }
    */

    std::ostream& operator<<(std::ostream&os, const cuda& im) {
        int i;
        std::ostringstream out;
        im.getMatrix();
        for(i = 0; i < im.size(); i++) {
            if(i > 0 && i % im.cols == 0)
                out << '\n';
            out << std::to_string(im.matrix_data[i]);
            out << '\t';
        }
        os << out.str();
        return os;
    }

}
