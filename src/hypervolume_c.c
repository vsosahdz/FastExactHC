/*
 * Base exact hypervolume (minimization) by pivot partitioning: the pivot box is added
 * and the rest of the region is split into d disjoint children, recursively.
 * Small regions are solved by full inclusion-exclusion. Shared helpers for all versions.
 */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

enum {
    HYPERVOLUME_OK = 0,
    HYPERVOLUME_INVALID_ARGUMENT = 1,
    HYPERVOLUME_ALLOCATION_ERROR = 2,
    HYPERVOLUME_NUMERICAL_ERROR = 3,
};

#define INCLUSION_EXCLUSION_LIMIT 6

static double box_volume(const double *point, const double *reference, size_t dimension) {
    double volume = 1.0;
    for (size_t objective = 0; objective < dimension; ++objective) {
        volume *= reference[objective] - point[objective];
    }
    return volume;
}

static int dominates(const double *left, const double *right, size_t dimension) {
    for (size_t objective = 0; objective < dimension; ++objective) {
        if (left[objective] > right[objective]) {
            return 0;
        }
    }
    return 1;
}

static size_t nondominated(
    const double *points,
    size_t count,
    size_t dimension,
    double *frontier
) {
    size_t frontier_count = 0;
    for (size_t candidate = 0; candidate < count; ++candidate) {
        const double *value = points + candidate * dimension;
        int is_dominated = 0;
        for (size_t index = 0; index < frontier_count; ++index) {
            if (dominates(
                    frontier + index * dimension,
                    value,
                    dimension
            )) {
                is_dominated = 1;
                break;
            }
        }
        if (is_dominated) {
            continue;
        }

        size_t index = 0;
        while (index < frontier_count) {
            if (dominates(value, frontier + index * dimension, dimension)) {
                --frontier_count;
                for (size_t objective = 0; objective < dimension; ++objective) {
                    frontier[index * dimension + objective] =
                        frontier[frontier_count * dimension + objective];
                }
                continue;
            }
            ++index;
        }
        for (size_t objective = 0; objective < dimension; ++objective) {
            frontier[frontier_count * dimension + objective] = value[objective];
        }
        ++frontier_count;
    }
    return frontier_count;
}

static double closed_form(const double *points, size_t count, const double *reference, size_t dimension) {
    double volume = 0.0;
    if (count >= sizeof(size_t) * CHAR_BIT) {
        return -1.0;
    }
    size_t subset_count = (size_t)1 << count;
    double intersection[dimension];

    for (size_t mask = 1; mask < subset_count; ++mask) {
        for (size_t objective = 0; objective < dimension; ++objective) {
            intersection[objective] = -DBL_MAX;
        }

        size_t selected = 0;
        for (size_t point = 0; point < count; ++point) {
            if ((mask & ((size_t)1 << point)) != 0) {
                ++selected;
                for (size_t objective = 0; objective < dimension; ++objective) {
                    double value = points[point * dimension + objective];
                    if (value > intersection[objective]) {
                        intersection[objective] = value;
                    }
                }
            }
        }

        double intersection_volume = box_volume(intersection, reference, dimension);
        volume += (selected % 2 == 1) ? intersection_volume : -intersection_volume;
    }
    return volume;
}

typedef struct {
    double **buffers;
    size_t *capacities;
    size_t levels;
} Workspace;

static double *workspace_buffer(
    Workspace *workspace,
    size_t level,
    size_t count,
    size_t dimension
) {
    if (dimension != 0 && count > SIZE_MAX / dimension) {
        return NULL;
    }
    size_t required = count * dimension;
    if (level >= workspace->levels) {
        size_t new_levels = level + 1;
        double **new_buffers = malloc(new_levels * sizeof(*new_buffers));
        size_t *new_capacities = malloc(new_levels * sizeof(*new_capacities));
        if (new_buffers == NULL || new_capacities == NULL) {
            free(new_buffers);
            free(new_capacities);
            return NULL;
        }
        for (size_t index = 0; index < workspace->levels; ++index) {
            new_buffers[index] = workspace->buffers[index];
            new_capacities[index] = workspace->capacities[index];
        }
        free(workspace->buffers);
        free(workspace->capacities);
        workspace->buffers = new_buffers;
        workspace->capacities = new_capacities;
        for (size_t index = workspace->levels; index < new_levels; ++index) {
            workspace->buffers[index] = NULL;
            workspace->capacities[index] = 0;
        }
        workspace->levels = new_levels;
    }
    if (required > workspace->capacities[level]) {
        double *buffer = realloc(
            workspace->buffers[level],
            required * sizeof(*buffer)
        );
        if (buffer == NULL) {
            return NULL;
        }
        workspace->buffers[level] = buffer;
        workspace->capacities[level] = required;
    }
    return workspace->buffers[level];
}

static void workspace_destroy(Workspace *workspace) {
    for (size_t level = 0; level < workspace->levels; ++level) {
        free(workspace->buffers[level]);
    }
    free(workspace->buffers);
    free(workspace->capacities);
}

typedef struct {
    const double *points;
    size_t count;
    const double *lower;
    const double *reference;
    size_t dimension;
} Region;

static double partition_region(
    const Region *region,
    Workspace *workspace,
    size_t level
) {
    size_t dimension = region->dimension;
    double *local_points = workspace_buffer(
        workspace,
        level,
        region->count,
        dimension
    );
    if (local_points == NULL) {
        return -1.0;
    }

    size_t local_count = 0;
    for (size_t point = 0; point < region->count; ++point) {
        const double *source = region->points + point * dimension;
        int inside = 1;
        for (size_t objective = 0; objective < dimension; ++objective) {
            if (source[objective] >= region->reference[objective]) {
                inside = 0;
                break;
            }
        }
        if (inside) {
            for (size_t objective = 0; objective < dimension; ++objective) {
                double value = source[objective] > region->lower[objective]
                    ? source[objective]
                    : region->lower[objective];
                local_points[local_count * dimension + objective] = value;
            }
            ++local_count;
        }
    }

    if (local_count == 0) {
        return 0.0;
    }
    if (local_count == 1) {
        double volume = box_volume(local_points, region->reference, dimension);
        return volume;
    }
    if (local_count <= INCLUSION_EXCLUSION_LIMIT) {
        double volume = closed_form(
            local_points,
            local_count,
            region->reference,
            dimension
        );
        return volume;
    }

    size_t pivot_index = 0;
    double pivot_volume = box_volume(local_points, region->reference, dimension);
    for (size_t point = 1; point < local_count; ++point) {
        double value = box_volume(
            local_points + point * dimension,
            region->reference,
            dimension
        );
        if (value > pivot_volume) {
            pivot_index = point;
            pivot_volume = value;
        }
    }

    double *pivot = local_points + pivot_index * dimension;
    double volume = pivot_volume;

    for (size_t objective = 0; objective < dimension; ++objective) {
        double child_lower[dimension];
        double child_reference[dimension];

        for (size_t index = 0; index < dimension; ++index) {
            child_lower[index] = region->lower[index];
            child_reference[index] = region->reference[index];
        }
        child_reference[objective] = pivot[objective];
        for (size_t previous = 0; previous < objective; ++previous) {
            if (child_lower[previous] < pivot[previous]) {
                child_lower[previous] = pivot[previous];
            }
        }

        int valid = 1;
        for (size_t index = 0; index < dimension; ++index) {
            if (child_lower[index] >= child_reference[index]) {
                valid = 0;
                break;
            }
        }
        if (valid) {
            Region child = {
                local_points,
                local_count,
                child_lower,
                child_reference,
                dimension,
            };
            double child_volume = partition_region(&child, workspace, level + 1);
            if (child_volume < 0.0) {
                return -1.0;
            }
            volume += child_volume;
        }
    }

    return volume;
}

static int valid_input(
    const double *points,
    size_t count,
    size_t dimension,
    const double *reference
) {
    if (points == NULL || reference == NULL || count == 0 || dimension == 0) {
        return 0;
    }
    for (size_t objective = 0; objective < dimension; ++objective) {
        if (!isfinite(reference[objective])) {
            return 0;
        }
    }
    for (size_t point = 0; point < count; ++point) {
        for (size_t objective = 0; objective < dimension; ++objective) {
            double value = points[point * dimension + objective];
            if (!isfinite(value) || value >= reference[objective]) {
                return 0;
            }
        }
    }
    return 1;
}

/* Returns -1.0 only when memory allocation or arithmetic fails. */
double hypervolume_c(
    const double *input_points,
    size_t count,
    size_t dimension,
    const double *reference
) {
    if (!valid_input(input_points, count, dimension, reference)) {
        return 0.0;
    }

    double *candidates = malloc(count * dimension * sizeof(double));
    double *frontier = malloc(count * dimension * sizeof(double));
    if (candidates == NULL || frontier == NULL) {
        free(candidates);
        free(frontier);
        return -1.0;
    }

    size_t candidate_count = 0;
    for (size_t point = 0; point < count; ++point) {
        int inside = 1;
        for (size_t objective = 0; objective < dimension; ++objective) {
            double value = input_points[point * dimension + objective];
            if (value >= reference[objective]) {
                inside = 0;
                break;
            }
        }
        if (inside) {
            for (size_t objective = 0; objective < dimension; ++objective) {
                candidates[candidate_count * dimension + objective] =
                    input_points[point * dimension + objective];
            }
            ++candidate_count;
        }
    }

    if (candidate_count == 0) {
        free(candidates);
        free(frontier);
        return 0.0;
    }

    double result;
    if (candidate_count <= INCLUSION_EXCLUSION_LIMIT) {
        result = closed_form(candidates, candidate_count, reference, dimension);
    } else {
        size_t frontier_count = nondominated(
            candidates,
            candidate_count,
            dimension,
            frontier
        );
        double lower[dimension];
        for (size_t objective = 0; objective < dimension; ++objective) {
            lower[objective] = frontier[objective];
            for (size_t point = 1; point < frontier_count; ++point) {
                double value = frontier[point * dimension + objective];
                if (value < lower[objective]) {
                    lower[objective] = value;
                }
            }
        }

        if (frontier_count <= INCLUSION_EXCLUSION_LIMIT) {
            result = closed_form(frontier, frontier_count, reference, dimension);
        } else {
            Region root = {frontier, frontier_count, lower, reference, dimension};
            Workspace workspace = {NULL, NULL, 0};
            result = partition_region(&root, &workspace, 0);
            workspace_destroy(&workspace);
        }
    }

    free(candidates);
    free(frontier);
    return result;
}

int hypervolume_c_status(
    const double *input_points,
    size_t count,
    size_t dimension,
    const double *reference,
    double *result
) {
    if (result == NULL || !valid_input(input_points, count, dimension, reference)) {
        return HYPERVOLUME_INVALID_ARGUMENT;
    }

    double value = hypervolume_c(input_points, count, dimension, reference);
    if (value == -1.0) {
        return HYPERVOLUME_ALLOCATION_ERROR;
    }
    if (!isfinite(value)) {
        return HYPERVOLUME_NUMERICAL_ERROR;
    }

    double tolerance = 64.0 * DBL_EPSILON * fmax(1.0, fabs(value));
    if (value < 0.0) {
        if (fabs(value) <= tolerance) {
            value = 0.0;
        } else {
            return HYPERVOLUME_NUMERICAL_ERROR;
        }
    }
    *result = value;
    return HYPERVOLUME_OK;
}
