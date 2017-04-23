#include "task.h"

int _DEBUG = 0, _ERROR = 0;

int length(const char *str) {
    int count = 0;
    while (*str++ != '\0') {
        count++;
    }
    return count;
}

int compare (const void *a, const void *b) {
    double _a = *(double *)a;
    double _b = *(double *)b;
    if (_a > _b) {
        return 1;
    } else if (_a < _b) {
        return -1;
    }
    return 0;
}

int _compare(const char *a, const char *b, int len) {
    while (len--) {
        if (*a++ != *b++) {
            return 1;
        }
    }
    return 0;
}

void print_help() {
    printf("Usage: app.out [input_filename] [output_filename] [options]\n"
           "Where options include:\n"
           "-d                print debug messages [default OFF]\n"
           "-e                print errors [default OFF]\n"
           "-p                print matrix [default OFF]\n"
           "-t                print execution time [default OFF]\n"
           "-prec=<num>       precision [default - 1e-14]\n"
           "-eps=<num>        'epsilon' [default - 1e-10]\n"
           "-max_iter=<num>   limit number of iterations [default - 0, i.e. not limit]\n"
           "-h, -?            print this and exit\n");
}

size_t total_memsize(n) {
    size_t a = sim_memsize(n);
    size_t b = evc_memsize(n);
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

int main(int argc, const char * argv[]) {
    // проверка на количество параметров
    if (argc > 11 || (argc > 9 && argv[1][0] == '-')) {
        printf("Too much parameters\n");
        return 1;
    }
    int max_iter = 0;
    double prec = 1e-14;
    double eps = 1e-10;
    // проверка на корректность параметров
    int* parameters = (int*)malloc(9 * sizeof(int));
    for (int i = 0; i < 9; i++) {
        parameters[i] = 0;
    }
    if (argc > 1) {
        int parameter_index = 1;
        if (argv[1][0] != '-') {
            parameter_index = 3;
        }

        for (; parameter_index < argc; parameter_index++) {
            if (_compare(argv[parameter_index], "-d", 2) == 0) {
                parameters[0]++;
            } else if (_compare(argv[parameter_index], "-e", 2) == 0 && _compare(argv[parameter_index], "-eps", 4) == 1) {
                parameters[1]++;
            } else if (_compare(argv[parameter_index], "-p", 2) == 0 && _compare(argv[parameter_index], "-prec", 5) == 1) {
                parameters[2]++;
            } else if (_compare(argv[parameter_index], "-t", 2) == 0) {
                parameters[3]++;
            } else if (_compare(argv[parameter_index], "-h", 2) == 0) {
                parameters[4]++;
            } else if (_compare(argv[parameter_index], "-?", 2) == 0) {
                parameters[5]++;
            } else if (_compare(argv[parameter_index], "-prec=", 6) == 0) {
                prec = atof(argv[parameter_index] + 6);
                parameters[6]++;
            } else if (_compare(argv[parameter_index], "-eps=", 5) == 0) {
                eps = atof(argv[parameter_index] + 5);
                parameters[7]++;
            } else if (_compare(argv[parameter_index], "-max_iter=", 10) == 0) {
                max_iter = atoi(argv[parameter_index] + 10);
                parameters[8]++;
            } else {
                printf("Unknown parameter %s\n", argv[parameter_index]);
                return 2;
            }
        }
        for (int i = 0; i < 9; i++) {
            if (parameters[i] > 1) {
                printf("Duplicate parameters\n");
                return 3;
            }
        }
    }
    if (parameters[4] + parameters[5] > 0) {
        print_help();
        return 0;
    }
    if (!prec) {
        if (_ERROR) printf("Invalid value for '-prec'\n");
        return -11;
    }
    if (!eps) {
        if (_ERROR) printf("Invalid value for '-eps'\n");
        return -12;
    }
    _DEBUG = parameters[0];
    _ERROR = parameters[1];

    const char *in_filename = "infile.txt";
    const char *out_filename = "outfile.txt";

    FILE *in_file;
    if (argc > 1) {
        if (argv[1][0] == '-') {
            in_file = fopen(in_filename, "r");
        } else {
            in_file = fopen(argv[1], "r");
        }
    } else {
        in_file = fopen(in_filename, "r");
    }

    if (in_file == NULL) {
        if (_ERROR) printf("Error on opening file\n");
        return 4;
    }

    int n;
    if (fscanf(in_file, "%d", &n) != 1) {
        if (_ERROR) printf("Error on reading file\n");
        return 5;
    }
    if (n < 0) {
        if (_ERROR) printf("Invalid matrix size\n");
        return 6;
    }

    double *A, *E, *tmp;
    if (!(A = malloc(n * n * sizeof(double)))) {
        if (_ERROR) printf("Memory allocation error\n");
        return 7;
    }
    if (!(E = malloc(n * sizeof(double)))) {
        if (_ERROR) printf("Memory allocation error\n");
        return 7;
    }
    if (!(tmp = malloc(total_memsize(n)))) {
        if (_ERROR) printf("Memory allocation error\n");
        return 7;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(in_file, "%lf", &A[i * n + j]) != 1) {
                if (_ERROR) printf("Error on reading file\n");
                return 5;
            }
        }
    }

    if (parameters[2]) {
        printf("А:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.9lf\t", A[i * n + j]);
            }
            printf("\n");
        }
    }

    if (fclose(in_file)) {
        if (_ERROR) printf("Error on closing file\n");
        return 8;
    }

    FILE *out_file;
    if (argc > 1) {
        if (argv[1][0] == '-') {
            out_file = fopen(out_filename, "wb");
        } else {
            out_file = fopen(argv[2], "wb");
        }
    } else {
        out_file = fopen(out_filename, "wb");
    }
    if (out_file == NULL) {
        if (_ERROR) printf("Error on opening file\n");
        return 4;
    }

    clock_t start = clock();
    int _sim = sim(n, A, tmp, eps);
    if (_sim != 0) {
        if (fprintf(out_file, "0") < 0) {
            if (_ERROR) printf("Error on writing into file\n");
            return 9;
        }

        if (fclose(out_file)) {
            if (_ERROR) printf("Error on closing file\n");
            return 8;
        }

        free(A);
        free(E);
        free(tmp);

        return 1;
    }
    if (parameters[2]) {
        printf("А:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.9lf\t", A[i * n + j]);
            }
            printf("\n");
        }
    }

    int _evc = evc(n, max_iter, eps, A, E, tmp, prec);
    if (parameters[3] == 1) printf("\nExecution time: %lf sec(s)\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    if (_evc != 0) {
        if (fprintf(out_file, "0") < 0) {
            if (_ERROR) printf("Error on writing into file\n");
            return 9;
        }

        if (fclose(out_file)) {
            if (_ERROR) printf("Error on closing file\n");
            return 8;
        }

        free(A);
        free(E);
        free(tmp);

        return 1;
    }

    if (fprintf(out_file, "%d\n", n) < 0) {
        if (_ERROR) printf("Error on writing into file\n");
        return 9;
    }
    qsort(E, n, sizeof(double), compare);
    for (int i = 0; i < n; i++) {
        if (fprintf(out_file, "%.9lf\n", E[i]) < 0) {
            if (_ERROR) printf("Error on writing into file\n");
            return 9;
        }
    }
    if (fclose(out_file)) {
        if (_ERROR) printf("Error on closing file\n");
        return 8;
    }

    free(A);
    free(E);
    free(tmp);

    return 0;
}
