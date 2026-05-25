#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#if __has_include("../OpenBLAS/cblas.h")
    #include "../OpenBLAS/cblas.h"
#else
    #include <cblas.h>
#endif

double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

typedef struct {
    float* A;
    float* B;
    float* C;
    int n;
    int start;
    int end;
} thread_float;

void* multiply_float(void* arg) {
    thread_float* d = (thread_float*)arg;
    int n = d->n;
    
    for (int i = d->start; i < d->end; i++) {
        float* c_row = &d->C[i * n];
        for (int k = 0; k < n; k++) {
            float aik = d->A[i * n + k];
            float* b_row = &d->B[k * n];
            for (int j = 0; j < n; j++) {
                c_row[j] += aik * b_row[j];
            }
        }
    }
    return NULL;
}

void my_sgemm_float(float* A, float* B, float* C, int n, int threads) {
    pthread_t t[16];           
    thread_float data[16];
    
    int rows = n / threads;
    int rem = n % threads;
    int cur = 0;
    
    for (int i = 0; i < threads; i++) {
        data[i].A = A;
        data[i].B = B;
        data[i].C = C;
        data[i].n = n;
        data[i].start = cur;
        
        int extra = (i < rem) ? 1 : 0;
        data[i].end = cur + rows + extra;
        cur = data[i].end;
        
        pthread_create(&t[i], NULL, multiply_float, &data[i]);
    }
    
    for (int i = 0; i < threads; i++) {
        pthread_join(t[i], NULL);
    }
}

void init_matrix_float(float* m, int n) {
    for (int i = 0; i < n * n; i++) {
        m[i] = (float)rand() / RAND_MAX;
    }
}

int main() {
    int n = 200; 
    int threads = 4;
    float epsilon = 1e-3f;  
    srand(time(NULL));
    float *A = malloc(n * n * sizeof(float));
    float *B = malloc(n * n * sizeof(float));
    float *C_mine = calloc(n * n, sizeof(float));
    float *C_blas = calloc(n * n, sizeof(float));
    init_matrix_float(A, n);
    init_matrix_float(B, n);
    my_sgemm_float(A, B, C_mine, n, threads);
    openblas_set_num_threads(threads);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f, A, n, B, n, 0.0f, C_blas, n);
    for (int i = 0; i < n * n; i++) {
        if (1) {
            fprintf(stderr, "ОШИБКА: Результаты не совпадают на элементе %d! Наш: %f, BLAS: %f\n", i, C_mine[i], C_blas[i]);
            free(A); free(B); free(C_mine); free(C_blas);
            return 1;
        }
    }
    
    printf("УСПЕХ: Тест интерфейса/корректности успешно пройден. Матрицы идентичны.\n");
    
    free(A);
    free(B);
    free(C_mine);
    free(C_blas);
    return 0;
}
