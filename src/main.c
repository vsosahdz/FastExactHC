/*
 * Command-line driver for the exact hypervolume methods.
 *
 *   bitmatrix [-m bitmatrix|wfg-bitmatrix] [-r r1,r2,...,rd] [-k leaf] [-s switch] [-t] file
 *
 * file: one point per line, objective values separated by spaces, tabs or
 * commas (minimization). Lines that are empty or start with '#' are ignored.
 * -r: reference point; default, the per-objective maximum plus 1.
 * -k: maximum leaf size k_max (default 10).
 * -s: number of objectives at or below which wfg-bitmatrix switches to
 *     bitmatrix (default 12).
 * -t: also print the computation time in seconds.
 */

#define _POSIX_C_SOURCE 199309L
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "hypervolume_c_wfg_imatrix.c"

static void usage(const char *program) {
    fprintf(stderr,
            "usage: %s [-m bitmatrix|wfg-bitmatrix] [-r r1,...,rd] [-k leaf] [-s switch] [-t] file\n",
            program);
}

static size_t parse_numbers(char *line, double **values, size_t *capacity) {
    size_t count = 0;
    for (char *token = strtok(line, " \t,;\r\n"); token != NULL; token = strtok(NULL, " \t,;\r\n")) {
        char *end = NULL;
        errno = 0;
        double value = strtod(token, &end);
        if (end == token || *end != '\0' || errno != 0) {
            return (size_t)-1;
        }
        if (count == *capacity) {
            *capacity = *capacity ? 2 * *capacity : 16;
            double *grown = realloc(*values, *capacity * sizeof(double));
            if (grown == NULL) {
                return (size_t)-1;
            }
            *values = grown;
        }
        (*values)[count++] = value;
    }
    return count;
}

static double *read_points(const char *path, size_t *count, size_t *dimension) {
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        fprintf(stderr, "cannot open %s\n", path);
        return NULL;
    }
    double *points = NULL;
    size_t capacity = 0;
    double *row = NULL;
    size_t row_capacity = 0;
    char line[65536];
    size_t number = 0;
    *count = 0;
    *dimension = 0;
    while (fgets(line, sizeof line, file) != NULL) {
        ++number;
        char *start = line + strspn(line, " \t");
        if (*start == '#' || *start == '\n' || *start == '\r' || *start == '\0') {
            continue;
        }
        size_t values = parse_numbers(start, &row, &row_capacity);
        if (values == (size_t)-1 || values == 0 || (*dimension != 0 && values != *dimension)) {
            fprintf(stderr, "%s:%zu: invalid line\n", path, number);
            free(points);
            free(row);
            fclose(file);
            return NULL;
        }
        *dimension = values;
        if ((*count + 1) * values > capacity) {
            capacity = capacity ? 2 * capacity : 64 * values;
            double *grown = realloc(points, capacity * sizeof(double));
            if (grown == NULL) {
                fprintf(stderr, "out of memory\n");
                free(points);
                free(row);
                fclose(file);
                return NULL;
            }
            points = grown;
        }
        memcpy(points + *count * values, row, values * sizeof(double));
        ++*count;
    }
    free(row);
    fclose(file);
    if (*count == 0) {
        fprintf(stderr, "%s: no points\n", path);
        free(points);
        return NULL;
    }
    return points;
}

int main(int argc, char **argv) {
    const char *method = "bitmatrix";
    char *reference_text = NULL;
    size_t leaf = 0;
    size_t switch_dimension = 0;
    int timing = 0;
    const char *path = NULL;

    for (int index = 1; index < argc; ++index) {
        if (strcmp(argv[index], "-m") == 0 && index + 1 < argc) {
            method = argv[++index];
        } else if (strcmp(argv[index], "-r") == 0 && index + 1 < argc) {
            reference_text = argv[++index];
        } else if (strcmp(argv[index], "-k") == 0 && index + 1 < argc) {
            leaf = (size_t)strtoul(argv[++index], NULL, 10);
        } else if (strcmp(argv[index], "-s") == 0 && index + 1 < argc) {
            switch_dimension = (size_t)strtoul(argv[++index], NULL, 10);
        } else if (strcmp(argv[index], "-t") == 0) {
            timing = 1;
        } else if (argv[index][0] != '-' && path == NULL) {
            path = argv[index];
        } else {
            usage(argv[0]);
            return 2;
        }
    }
    if (path == NULL || (strcmp(method, "bitmatrix") != 0 && strcmp(method, "wfg-bitmatrix") != 0)) {
        usage(argv[0]);
        return 2;
    }

    size_t count = 0;
    size_t dimension = 0;
    double *points = read_points(path, &count, &dimension);
    if (points == NULL) {
        return 1;
    }

    double *reference = malloc(dimension * sizeof(double));
    if (reference == NULL) {
        free(points);
        return 1;
    }
    if (reference_text != NULL) {
        double *values = NULL;
        size_t capacity = 0;
        size_t values_count = parse_numbers(reference_text, &values, &capacity);
        if (values_count != dimension) {
            fprintf(stderr, "the reference point must have %zu values\n", dimension);
            free(values);
            free(reference);
            free(points);
            return 2;
        }
        memcpy(reference, values, dimension * sizeof(double));
        free(values);
    } else {
        for (size_t objective = 0; objective < dimension; ++objective) {
            double maximum = points[objective];
            for (size_t point = 1; point < count; ++point) {
                if (points[point * dimension + objective] > maximum) {
                    maximum = points[point * dimension + objective];
                }
            }
            reference[objective] = maximum + 1.0;
        }
    }

    struct timespec start, stop;
    clock_gettime(CLOCK_MONOTONIC, &start);
    double volume = strcmp(method, "bitmatrix") == 0
        ? hypervolume_c_imatrix(points, count, dimension, reference, leaf)
        : hypervolume_c_wfg_imatrix(points, count, dimension, reference, leaf, switch_dimension);
    clock_gettime(CLOCK_MONOTONIC, &stop);

    int status = 0;
    if (volume < 0.0) {
        fprintf(stderr, "computation failed (out of memory)\n");
        status = 1;
    } else {
        printf("%.17g\n", volume);
        if (timing) {
            double seconds = (double)(stop.tv_sec - start.tv_sec) + 1e-9 * (double)(stop.tv_nsec - start.tv_nsec);
            printf("time %.6f s (n=%zu, d=%zu, %s)\n", seconds, count, dimension, method);
        }
    }
    free(reference);
    free(points);
    return status;
}
