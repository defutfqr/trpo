#ifndef CBLAS_H
#define CBLAS_H
#include <stddef.h>

#define CBLAS_INDEX size_t

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

#ifdef __cplusplus
extern "C" {
#endif

void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *B, const int ldb, const float beta,
                 float *C, const int ldc);
void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc);
void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const float alpha, const float *A, const int lda,
                  const float *B, const int ldb, const float beta,
                  float *C, const int ldc);
void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);
void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb);

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc);
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double beta,
                  double *C, const int ldc);
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);

void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);

void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc);
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc);
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc);
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb);

void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const void *A, const int lda,
                 const float beta, void *C, const int ldc);
void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const float beta,
                  void *C, const int ldc);

void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc);
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const void *A, const int lda,
                 const double beta, void *C, const int ldc);
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const double beta,
                  void *C, const int ldc);

void cblas_xerbla(int p, const char *rout, const char *form, ...);

#ifdef __cplusplus
}
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define FLOAT_TOL 1e-5
#define DOUBLE_TOL 1e-10

#define ASSERT_NEAR(a, b, tol, test_name) \
    do { \
        tests_run++; \
        if (fabs((a) - (b)) < (tol)) { \
            tests_passed++; \
            printf("  [PASS] %s\n", test_name); \
        } else { \
            tests_failed++; \
            printf("  [FAIL] %s: expected %f, got %f\n", test_name, (double)(b), (double)(a)); \
        } \
    } while(0)

#define ASSERT_COMPLEX_NEAR(a, b, tol, test_name) \
    do { \
        tests_run++; \
        double diff = cabs((a) - (b)); \
        if (diff < (tol)) { \
            tests_passed++; \
            printf("  [PASS] %s\n", test_name); \
        } else { \
            tests_failed++; \
            printf("  [FAIL] %s: expected %f+%fi, got %f+%fi\n", \
                   test_name, creal(b), cimag(b), creal(a), cimag(a)); \
        } \
    } while(0)

void test_sgemm_basic(void) {
    printf("Testing cblas_sgemm - basic multiplication\n");
    
    float A[6] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    float B[6] = {7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 3, alpha, A, 3, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 58.0f, FLOAT_TOL, "C[0,0]");
    ASSERT_NEAR(C[1], 64.0f, FLOAT_TOL, "C[0,1]");
    ASSERT_NEAR(C[2], 139.0f, FLOAT_TOL, "C[1,0]");
    ASSERT_NEAR(C[3], 154.0f, FLOAT_TOL, "C[1,1]");
}

void test_sgemm_transpose(void) {
    printf("Testing cblas_sgemm - with transpose\n");
    
    float A[6] = {1.0f, 3.0f, 5.0f, 2.0f, 4.0f, 6.0f};
    float B[6] = {7.0f, 9.0f, 11.0f, 8.0f, 10.0f, 12.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans,
                2, 2, 3, alpha, A, 3, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 58.0f, FLOAT_TOL, "C[0,0] with transpose");
    ASSERT_NEAR(C[1], 64.0f, FLOAT_TOL, "C[0,1] with transpose");
}

void test_sgemm_with_beta(void) {
    printf("Testing cblas_sgemm - with non-zero beta\n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    float alpha = 2.0f, beta = 3.0f;
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 2.0f*19.0f + 3.0f*1.0f, FLOAT_TOL, "C[0,0] with beta");
}

void test_ssymm_left(void) {
    printf("Testing cblas_ssymm - left multiplication\n");
    
    float A[4] = {1.0f, 2.0f, 2.0f, 3.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_ssymm(CblasRowMajor, CblasLeft, CblasLower,
                2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 1.0f*1.0f + 2.0f*3.0f, FLOAT_TOL, "C[0,0] symm left");
    ASSERT_NEAR(C[1], 1.0f*2.0f + 2.0f*4.0f, FLOAT_TOL, "C[0,1] symm left");
}

void test_ssymm_right(void) {
    printf("Testing cblas_ssymm - right multiplication\n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {1.0f, 2.0f, 2.0f, 3.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_ssymm(CblasRowMajor, CblasRight, CblasLower,
                2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 1.0f*1.0f + 2.0f*2.0f, FLOAT_TOL, "C[0,0] symm right");
}

void test_ssyrk_upper(void) {
    printf("Testing cblas_ssyrk - upper triangular\n");
    
    float A[6] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_ssyrk(CblasRowMajor, CblasUpper, CblasNoTrans,
                2, 3, alpha, A, 3, beta, C, 2);
    
    ASSERT_NEAR(C[0], 1.0f*1.0f + 2.0f*2.0f + 3.0f*3.0f, FLOAT_TOL, "C[0,0] syrk");
    ASSERT_NEAR(C[1], 1.0f*4.0f + 2.0f*5.0f + 3.0f*6.0f, FLOAT_TOL, "C[0,1] syrk");
    ASSERT_NEAR(C[3], 4.0f*4.0f + 5.0f*5.0f + 6.0f*6.0f, FLOAT_TOL, "C[1,1] syrk");
}

void test_ssyr2k_basic(void) {
    printf("Testing cblas_ssyr2k - basic\n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_ssyr2k(CblasRowMajor, CblasUpper, CblasNoTrans,
                 2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 2.0f*(1.0f*5.0f + 2.0f*6.0f), FLOAT_TOL, "C[0,0] syr2k");
}

void test_strmm_left(void) {
    printf("Testing cblas_strmm - left multiplication\n");
    
    float A[4] = {1.0f, 0.0f, 2.0f, 3.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float alpha = 2.0f;
    
    cblas_strmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                2, 2, alpha, A, 2, B, 2);
    
    ASSERT_NEAR(B[0], 2.0f*1.0f*1.0f, FLOAT_TOL, "B[0,0] trmm left");
    ASSERT_NEAR(B[1], 2.0f*(2.0f*1.0f + 3.0f*3.0f), FLOAT_TOL, "B[1,0] trmm left");
}

void test_strsm_left(void) {
    printf("Testing cblas_strsm - left solve\n");
    
    float A[4] = {2.0f, 0.0f, 4.0f, 3.0f};
    float B[4] = {2.0f, 4.0f, 14.0f, 12.0f};
    float alpha = 1.0f;
    
    cblas_strsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                2, 2, alpha, A, 2, B, 2);
    
    ASSERT_NEAR(B[0], 1.0f, FLOAT_TOL, "B[0,0] trsm left");
    ASSERT_NEAR(B[2], 2.0f, FLOAT_TOL, "B[0,1] trsm left");
}

void test_dgemm_basic(void) {
    printf("Testing cblas_dgemm - basic multiplication\n");
    
    double A[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double B[6] = {7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};
    double alpha = 1.0, beta = 0.0;
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 3, alpha, A, 3, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 58.0, DOUBLE_TOL, "C[0,0] dgemm");
    ASSERT_NEAR(C[1], 64.0, DOUBLE_TOL, "C[0,1] dgemm");
    ASSERT_NEAR(C[2], 139.0, DOUBLE_TOL, "C[1,0] dgemm");
    ASSERT_NEAR(C[3], 154.0, DOUBLE_TOL, "C[1,1] dgemm");
}

void test_dsymm_basic(void) {
    printf("Testing cblas_dsymm - basic\n");
    
    double A[4] = {1.0, 2.0, 2.0, 3.0};
    double B[4] = {1.0, 2.0, 3.0, 4.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};
    double alpha = 1.0, beta = 0.0;
    
    cblas_dsymm(CblasRowMajor, CblasLeft, CblasLower,
                2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 7.0, DOUBLE_TOL, "C[0,0] dsymm");
}

void test_dsyrk_basic(void) {
    printf("Testing cblas_dsyrk - basic\n");
    
    double A[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};
    double alpha = 1.0, beta = 0.0;
    
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans,
                2, 3, alpha, A, 3, beta, C, 2);
    
    ASSERT_NEAR(C[0], 14.0, DOUBLE_TOL, "C[0,0] dsyrk");
    ASSERT_NEAR(C[3], 77.0, DOUBLE_TOL, "C[1,1] dsyrk");
}

void test_dtrmm_basic(void) {
    printf("Testing cblas_dtrmm - basic\n");
    
    double A[4] = {1.0, 0.0, 2.0, 3.0};
    double B[4] = {1.0, 2.0, 3.0, 4.0};
    double alpha = 1.0;
    
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                2, 2, alpha, A, 2, B, 2);
    
    ASSERT_NEAR(B[0], 1.0, DOUBLE_TOL, "B[0,0] dtrmm");
}

void test_dtrsm_basic(void) {
    printf("Testing cblas_dtrsm - basic\n");
    
    double A[4] = {2.0, 0.0, 4.0, 3.0};
    double B[4] = {2.0, 4.0, 14.0, 12.0};
    double alpha = 1.0;
    
    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                2, 2, alpha, A, 2, B, 2);
    
    ASSERT_NEAR(B[0], 1.0, DOUBLE_TOL, "B[0,0] dtrsm");
}

void test_cgemm_basic(void) {
    printf("Testing cblas_cgemm - basic multiplication\n");
    
    float _Complex A[4] = {1.0f+0.0f*I, 2.0f+0.0f*I, 3.0f+0.0f*I, 4.0f+0.0f*I};
    float _Complex B[4] = {5.0f+0.0f*I, 6.0f+0.0f*I, 7.0f+0.0f*I, 8.0f+0.0f*I};
    float _Complex C[4] = {0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I};
    float _Complex alpha = 1.0f+0.0f*I, beta = 0.0f+0.0f*I;
    
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, &alpha, A, 2, B, 2, &beta, C, 2);
    
    ASSERT_COMPLEX_NEAR(C[0], 19.0f+0.0f*I, FLOAT_TOL, "C[0,0] cgemm");
    ASSERT_COMPLEX_NEAR(C[1], 22.0f+0.0f*I, FLOAT_TOL, "C[0,1] cgemm");
}

void test_csymm_basic(void) {
    printf("Testing cblas_csymm - basic\n");
    
    float _Complex A[4] = {1.0f+0.0f*I, 2.0f+0.0f*I, 2.0f+0.0f*I, 3.0f+0.0f*I};
    float _Complex B[4] = {1.0f+0.0f*I, 2.0f+0.0f*I, 3.0f+0.0f*I, 4.0f+0.0f*I};
    float _Complex C[4] = {0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I};
    float _Complex alpha = 1.0f+0.0f*I, beta = 0.0f+0.0f*I;
    
    cblas_csymm(CblasRowMajor, CblasLeft, CblasLower,
                2, 2, &alpha, A, 2, B, 2, &beta, C, 2);
    
    ASSERT_COMPLEX_NEAR(C[0], 7.0f+0.0f*I, FLOAT_TOL, "C[0,0] csymm");
}

void test_chemm_basic(void) {
    printf("Testing cblas_chemm - basic Hermitian\n");
    
    float _Complex A[4] = {1.0f+0.0f*I, 2.0f-3.0f*I, 2.0f+3.0f*I, 4.0f+0.0f*I};
    float _Complex B[4] = {1.0f+0.0f*I, 2.0f+0.0f*I, 3.0f+0.0f*I, 4.0f+0.0f*I};
    float _Complex C[4] = {0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I};
    float _Complex alpha = 1.0f+0.0f*I, beta = 0.0f+0.0f*I;
    
    cblas_chemm(CblasRowMajor, CblasLeft, CblasLower,
                2, 2, &alpha, A, 2, B, 2, &beta, C, 2);
    
    ASSERT_COMPLEX_NEAR(C[0], 1.0f*(1.0f+0.0f*I) + (2.0f-3.0f*I)*(3.0f+0.0f*I), 
                        FLOAT_TOL, "C[0,0] chemm");
}

void test_cherk_basic(void) {
    printf("Testing cblas_cherk - basic Hermitian rank-k\n");
    
    float _Complex A[4] = {1.0f+1.0f*I, 2.0f+0.0f*I, 3.0f-1.0f*I, 4.0f+0.0f*I};
    float _Complex C[4] = {0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I, 0.0f+0.0f*I};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_cherk(CblasRowMajor, CblasUpper, CblasNoTrans,
                2, 2, alpha, A, 2, beta, C, 2);
    
    ASSERT_NEAR(creal(C[0]), 2.0f, FLOAT_TOL, "C[0,0] real cherk");
}

void test_zgemm_basic(void) {
    printf("Testing cblas_zgemm - basic multiplication\n");
    
    double _Complex A[4] = {1.0+0.0*I, 2.0+0.0*I, 3.0+0.0*I, 4.0+0.0*I};
    double _Complex B[4] = {5.0+0.0*I, 6.0+0.0*I, 7.0+0.0*I, 8.0+0.0*I};
    double _Complex C[4] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
    double _Complex alpha = 1.0+0.0*I, beta = 0.0+0.0*I;
    
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, &alpha, A, 2, B, 2, &beta, C, 2);
    
    ASSERT_COMPLEX_NEAR(C[0], 19.0+0.0*I, DOUBLE_TOL, "C[0,0] zgemm");
    ASSERT_COMPLEX_NEAR(C[1], 22.0+0.0*I, DOUBLE_TOL, "C[0,1] zgemm");
}

void test_zsymm_basic(void) {
    printf("Testing cblas_zsymm - basic\n");
    
    double _Complex A[4] = {1.0+0.0*I, 2.0+0.0*I, 2.0+0.0*I, 3.0+0.0*I};
    double _Complex B[4] = {1.0+0.0*I, 2.0+0.0*I, 3.0+0.0*I, 4.0+0.0*I};
    double _Complex C[4] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
    double _Complex alpha = 1.0+0.0*I, beta = 0.0+0.0*I;
    
    cblas_zsymm(CblasRowMajor, CblasLeft, CblasLower,
                2, 2, &alpha, A, 2, B, 2, &beta, C, 2);
    
    ASSERT_COMPLEX_NEAR(C[0], 7.0+0.0*I, DOUBLE_TOL, "C[0,0] zsymm");
}

void test_zhemm_basic(void) {
    printf("Testing cblas_zhemm - basic Hermitian\n");
    
    double _Complex A[4] = {1.0+0.0*I, 2.0-3.0*I, 2.0+3.0*I, 4.0+0.0*I};
    double _Complex B[4] = {1.0+0.0*I, 2.0+0.0*I, 3.0+0.0*I, 4.0+0.0*I};
    double _Complex C[4] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
    double _Complex alpha = 1.0+0.0*I, beta = 0.0+0.0*I;
    
    cblas_zhemm(CblasRowMajor, CblasLeft, CblasLower,
                2, 2, &alpha, A, 2, B, 2, &beta, C, 2);
    
    ASSERT_COMPLEX_NEAR(C[0], 1.0*(1.0+0.0*I) + (2.0-3.0*I)*(3.0+0.0*I), 
                        DOUBLE_TOL, "C[0,0] zhemm");
}

void test_zherk_basic(void) {
    printf("Testing cblas_zherk - basic Hermitian rank-k\n");
    
    double _Complex A[4] = {1.0+1.0*I, 2.0+0.0*I, 3.0-1.0*I, 4.0+0.0*I};
    double _Complex C[4] = {0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I, 0.0+0.0*I};
    double alpha = 1.0, beta = 0.0;
    
    cblas_zherk(CblasRowMajor, CblasUpper, CblasNoTrans,
                2, 2, alpha, A, 2, beta, C, 2);
    
    ASSERT_NEAR(creal(C[0]), 2.0, DOUBLE_TOL, "C[0,0] real zherk");
}

void test_edge_zero_dimensions(void) {
    printf("Testing edge case - zero dimensions\n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {10.0f, 10.0f, 10.0f, 10.0f};
    float alpha = 1.0f, beta = 1.0f;
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                0, 2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 10.0f, FLOAT_TOL, "C[0,0] with M=0");
}

void test_edge_alpha_zero(void) {
    printf("Testing edge case - alpha=0\n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float alpha = 0.0f, beta = 1.0f;
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 1.0f, FLOAT_TOL, "C[0,0] with alpha=0");
}

void test_edge_beta_zero(void) {
    printf("Testing edge case - beta=0\n");
    
    float A[4] = {1.0f, 0.0f, 0.0f, 1.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {100.0f, 100.0f, 100.0f, 100.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 5.0f, FLOAT_TOL, "C[0,0] with beta=0");
}

void test_colmajor_layout(void) {
    printf("Testing column-major layout\n");
    
    float A[4] = {1.0f, 3.0f, 2.0f, 4.0f};
    float B[4] = {5.0f, 7.0f, 6.0f, 8.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float alpha = 1.0f, beta = 0.0f;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, alpha, A, 2, B, 2, beta, C, 2);
    
    ASSERT_NEAR(C[0], 19.0f, FLOAT_TOL, "C[0,0] col-major");
    ASSERT_NEAR(C[1], 43.0f, FLOAT_TOL, "C[1,0] col-major");
}

#ifdef OPENBLAS_THREAD_SUPPORT
#include <omp.h>

void test_multithreaded_gemm(void) {
    printf("Testing multi-threaded GEMM\n");
    
    int M = 100, N = 100, K = 100;
    float *A = malloc(M * K * sizeof(float));
    float *B = malloc(K * N * sizeof(float));
    float *C = malloc(M * N * sizeof(float));
    
    for (int i = 0; i < M * K; i++) A[i] = 1.0f;
    for (int i = 0; i < K * N; i++) B[i] = 1.0f;
    for (int i = 0; i < M * N; i++) C[i] = 0.0f;
    
    float alpha = 1.0f, beta = 0.0f;
    
    omp_set_num_threads(4);
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, alpha, A, K, B, N, beta, C, N);
    
    int correct = 1;
    for (int i = 0; i < M * N; i++) {
        if (fabs(C[i] - K) > FLOAT_TOL) {
            correct = 0;
            break;
        }
    }
    
    tests_run++;
    if (correct) {
        tests_passed++;
        printf("  [PASS] Multi-threaded GEMM\n");
    } else {
        tests_failed++;
        printf("  [FAIL] Multi-threaded GEMM\n");
    }
    
    free(A);
    free(B);
    free(C);
}
#endif

void run_all_tests(void) {
    printf("\n=== Level 3 BLAS Interface Tests ===\n\n");
    
    printf("--- Single Precision (S) Tests ---\n");
    test_sgemm_basic();
    test_sgemm_transpose();
    test_sgemm_with_beta();
    test_ssymm_left();
    test_ssymm_right();
    test_ssyrk_upper();
    test_ssyr2k_basic();
    test_strmm_left();
    test_strsm_left();
    printf("\n");
    
    printf("--- Double Precision (D) Tests ---\n");
    test_dgemm_basic();
    test_dsymm_basic();
    test_dsyrk_basic();
    test_dtrmm_basic();
    test_dtrsm_basic();
    printf("\n");
    
    printf("--- Single Precision Complex (C) Tests ---\n");
    test_cgemm_basic();
    test_csymm_basic();
    test_chemm_basic();
    test_cherk_basic();
    printf("\n");
    
    printf("--- Double Precision Complex (Z) Tests ---\n");
    test_zgemm_basic();
    test_zsymm_basic();
    test_zhemm_basic();
    test_zherk_basic();
    printf("\n");
    
    printf("--- Edge Cases and Special Tests ---\n");
    test_edge_zero_dimensions();
    test_edge_alpha_zero();
    test_edge_beta_zero();
    test_colmajor_layout();
    printf("\n");
    
#ifdef OPENBLAS_THREAD_SUPPORT
    printf("--- Multi-threaded Tests ---\n");
    test_multithreaded_gemm();
    printf("\n");
#endif
}

void print_summary(void) {
    printf("\n=== Test Summary ===\n");
    printf("Total tests run:    %d\n", tests_run);
    printf("Tests passed:       %d\n", tests_passed);
    printf("Tests failed:       %d\n", tests_failed);
    printf("Coverage:           %.1f%%\n", 
           100.0 * tests_passed / (tests_run > 0 ? tests_run : 1));
    
    if (tests_failed == 0) {
        printf("\n✓ All tests PASSED!\n");
    } else {
        printf("\n✗ Some tests FAILED!\n");
    }
}

int main(int argc, char *argv[]) {
    printf("OpenBLAS Level 3 CBLAS Interface Tests\n");
    printf("======================================\n");
    
#ifdef OPENBLAS_THREAD_SUPPORT
    int num_threads = 1;
    if (argc > 1) {
        num_threads = atoi(argv[1]);
        if (num_threads < 1) num_threads = 1;
        omp_set_num_threads(num_threads);
        printf("Using %d thread(s)\n\n", num_threads);
    }
#endif
    
    run_all_tests();
    print_summary();
    
    return (tests_failed == 0) ? 0 : 1;
}