#include "task.h"

size_t sim_memsize(int n) {
    return (n * n + n) * sizeof(double);
}

void copy_tmp_to_A(double *A, double *tmp, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = tmp[i * n + j];
        }
    }
}

int sim(int n, double* A, double* tmp, double precision) {
    double *xa = &tmp[n * n];
    for (int step = 0; step < n - 2; step++) {
        int indent = step + 1;
        if (_DEBUG) printf("\n**Step #%d**\n", indent);
        
        if (_DEBUG) printf("Extracting vector a and computing it's norm\n");
        xa[step] = 0.;
        double norm_a = 0.;
        for (int i = indent; i < n; i++) {
            xa[i] = A[i * n + step];
            norm_a += xa[i] * xa[i];
        }
        norm_a = sqrt(norm_a);
        
        if (_DEBUG) printf("Computing vector x\n");
        xa[indent] -= norm_a;
        double norm_x = 0.;
        for (int i = indent; i < n; i++) {
            norm_x += xa[i] * xa[i];
        }
        
        norm_x = sqrt(norm_x);
        if (fabs(norm_x) < precision) {
            continue;
        }
        for (int i = indent; i < n; i++) {
            xa[i] /= norm_x;
        }
        
        if (_DEBUG) printf("Computing A_%d\n", indent);
        if (_DEBUG) printf("U * A\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp[i * n + j] = 0.;
                
                for (int k = 0; k < n; k++) {
                    double U;
                    if (i == k) {
                        U = 1 - 2 * xa[i] * xa[k];
                    } else {
                        U = -2 * xa[i] * xa[k];
                    }
                    tmp[i * n + j] += U * A[k * n + j];
                }
                
            }
        }
        copy_tmp_to_A(A, tmp, n);
        
        if (_DEBUG) printf("U * A * U\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp[i * n + j] = 0.;
                
                for (int k = 0; k < n; k++) {
                    double U;
                    if (k == j) {
                        U = 1 - 2 * xa[k] * xa[j];
                    } else {
                        U = -2 * xa[k] * xa[j];
                    }
                    tmp[i * n + j] += A[i * n + k] * U;
                }
                
            }
        }
        copy_tmp_to_A(A, tmp, n);
        
    }
    
    if (_DEBUG) {
        printf("\nА:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.9lf\t", A[i * n + j]);
            }
            printf("\n");
        }
    }
    
    return 0;
}
