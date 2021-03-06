#include "task.h"

size_t evc_memsize(int n) {
    return n * n * sizeof(long double);
}

long double matrix_norm(long double *A, int n) {
    long double result = 0.;
    for (int i = 0; i < n; i++) {
        long double sum = 0.;
        for (int j = 0; j < n; j++) {
            sum += fabsl(A[i*n + j]);
        }
        
        if (result < sum) {
            result = sum;
        }
    }
    
    return result;
}

int evc(int n, int max_iterations, long double epsilon, long double* A, long double* E, long double* tmp, long double precision) {
    int _n = n;
    
    for (int iteration = 0; _n > 2; iteration++) {
        if (_DEBUG) printf("\n**Iteration #%d**\n", iteration);
        if (_DEBUG) printf("Current matrix size = %d\n", _n);
        long double e = precision * matrix_norm(A, n);
        
        if (_DEBUG) printf("Exhaustion\n");
        for (int i = 0; i < _n - 1; i++) {
            if (fabsl(A[(i + 1) * n + i]) < e) {
                int nn = i + 1;
                
                for (int _i = 0, item = 0; _i <= i; _i++) {
                    for (int _j = 0; _j <= i; _j++, item++) {
                        A[item] = A[_i * n + _j];
                    }
                }
                
                evc(nn, max_iterations - iteration, epsilon, A, E, tmp, precision);
                
                int nnn = _n - nn;
                long double *new_A = &A[(i + 1) * n + i + 1];
                long double *new_E = &E[nn];
                
                for (int _i = i + 1, item = 0; _i < _n; _i++) {
                    for (int _j = i + 1; _j < _n; _j++, item++) {
                        new_A[item] = A[_i * n + _j];
                    }
                }
                
                evc(nnn, max_iterations - iteration, epsilon, new_A, new_E, tmp, precision);
                return 0;
            }
        }
        
        long double sdvig = A[(_n - 1) * n + _n - 1];
        if (fabsl(A[0] - sdvig) < epsilon) {
            sdvig += 0.1;
        }
        if (_DEBUG) printf("Shifting matrix (s = %.9LF)\n", sdvig);
        for (int i = 0; i < _n; i++) {
            A[i * n + i] -= sdvig;
        }
        
        if (_DEBUG) printf("Computing LR decomposition\n");
        for (int k = 0; k < _n; k++) {
            tmp[0 * n + k] = A[0 * n + k];
        }
        for (int i = 1; i < _n; i++) {
            if (fabsl(tmp[(i - 1) * n + i - 1]) < epsilon) {
                if (_DEBUG || _ERROR) {
                    printf("ERROR: division by zero");
                    return -1;
                }
            }
            tmp[i * n + i - 1] = A[i * n + i - 1] / tmp[(i - 1) * n + i - 1];
            for (int k = i; k < _n; k++) {
                tmp[i * n + k] = A[i * n + k] - tmp[i * n + i - 1] * tmp[(i - 1) * n + k];
            }
        }
        
        if (_DEBUG) printf("Computing R*L\n");
        for (int i = 0; i < _n; i++) {
            A[i * n + _n - 1] = tmp[i * n + _n - 1];
            for (int k = i; k < _n - 1; k++) {
                A[i * n + k] = tmp[i * n + k] + tmp[i * n + k + 1] * tmp[(k + 1) * n + k];
            }
            if (i != 0) {
                A[i * n + i - 1] = tmp[i * n + i] * tmp[i * n + i - 1];
            }
        }
        
        if (_DEBUG) printf("Shifting back\n");
        for (int i = 0; i < _n; i++) {
            A[i * n + i] += sdvig;
        }
        
        if (fabsl(A[(_n - 1) * n + _n - 2]) < e) {
            E[_n - 1] = A[(_n - 1) * n + _n - 1];
            _n--;
        }
        
        if (_DEBUG) {
            printf("\nА:\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.9LF\t", A[i * n + j]);
                }
                printf("\n");
            }
        }
        if (max_iterations > 0 && iteration > max_iterations) {
            if (_DEBUG) printf("The number of iterations has exceeded the maximum value\n");
            return 1;
        }
    }
    if (n < 2) {
        E[0] = A[0];
        return 0;
    }
    if (_DEBUG) printf("\n2х2 Matrix left, solving a quadratic equation\n");
    // x^2 + p*x + q = 0
    long double p = -A[0 * n + 0] - A[1 * n + 1];
    long double q = A[0 * n + 0] * A[1 * n + 1] - A[0 * n + 1] * A[1 * n + 0];
    E[0] = (-p - sqrt(p * p - 4 * q)) / 2;
    E[1] = (-p + sqrt(p * p - 4 * q)) / 2;
    
    return 0;
}
