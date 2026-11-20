/*
 * Pivot partitioning with a pruned inclusion-exclusion in the leaves: when a later
 * point covers the join of a subset, the whole subtree of terms cancels and is skipped.
 * Points are sorted by ascending box volume so that large boxes can cover the joins.
 */
#include "hypervolume_c.c"

/*
 * Exact HV with a pruned inclusion-exclusion base case.
 *
 * Subsets are enumerated depth-first in index order with a running join
 * (componentwise maximum). If a later point p_j satisfies p_j <= join(S),
 * every subset in the subtree of S pairs with the same subset plus j: both
 * have the same intersection box and opposite signs, so the subtree sums
 * to zero and is skipped.
 */

static void sort_by_volume_ascending(
    double *points,
    size_t count,
    const double *reference,
    size_t dimension
) {
    double current[dimension];
    for (size_t point = 1; point < count; ++point) {
        for (size_t objective = 0; objective < dimension; ++objective) {
            current[objective] = points[point * dimension + objective];
        }
        double current_volume = box_volume(current, reference, dimension);

        size_t index = point;
        while (index > 0 && box_volume(
                points + (index - 1) * dimension,
                reference,
                dimension
        ) > current_volume) {
            for (size_t objective = 0; objective < dimension; ++objective) {
                points[index * dimension + objective] =
                    points[(index - 1) * dimension + objective];
            }
            --index;
        }
        for (size_t objective = 0; objective < dimension; ++objective) {
            points[index * dimension + objective] = current[objective];
        }
    }
}

static double pruned_inclusion_exclusion_subtree(
    const double *points,
    size_t count,
    const double *reference,
    size_t dimension,
    const double *join,
    size_t next,
    double sign
) {
    for (size_t point = next; point < count; ++point) {
        if (dominates(points + point * dimension, join, dimension)) {
            return 0.0;
        }
    }

    double volume = sign * box_volume(join, reference, dimension);
    double child_join[dimension];
    for (size_t point = next; point < count; ++point) {
        for (size_t objective = 0; objective < dimension; ++objective) {
            double value = points[point * dimension + objective];
            child_join[objective] = value > join[objective] ? value : join[objective];
        }
        volume += pruned_inclusion_exclusion_subtree(
            points,
            count,
            reference,
            dimension,
            child_join,
            point + 1,
            -sign
        );
    }
    return volume;
}

/* Sorts points in place; callers pass buffers they no longer need ordered. */
static double pruned_closed_form(
    double *points,
    size_t count,
    const double *reference,
    size_t dimension
) {
    sort_by_volume_ascending(points, count, reference, dimension);
    double volume = 0.0;
    for (size_t point = 0; point < count; ++point) {
        volume += pruned_inclusion_exclusion_subtree(
            points,
            count,
            reference,
            dimension,
            points + point * dimension,
            point + 1,
            1.0
        );
    }
    return volume;
}

static double pruned_partition_region(
    const Region *region,
    Workspace *workspace,
    size_t level,
    size_t exact_limit
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
        return box_volume(local_points, region->reference, dimension);
    }
    if (local_count <= exact_limit) {
        return pruned_closed_form(
            local_points,
            local_count,
            region->reference,
            dimension
        );
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
            double child_volume = pruned_partition_region(
                &child,
                workspace,
                level + 1,
                exact_limit
            );
            if (child_volume < 0.0) {
                return -1.0;
            }
            volume += child_volume;
        }
    }

    return volume;
}

/* Same contract as hypervolume_c; exact_limit = 0 uses INCLUSION_EXCLUSION_LIMIT. */
double hypervolume_c_ie_pruned(
    const double *input_points,
    size_t count,
    size_t dimension,
    const double *reference,
    size_t exact_limit
) {
    if (!valid_input(input_points, count, dimension, reference)) {
        return 0.0;
    }
    if (exact_limit == 0) {
        exact_limit = INCLUSION_EXCLUSION_LIMIT;
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
    size_t frontier_count = nondominated(
        candidates,
        candidate_count,
        dimension,
        frontier
    );
    if (frontier_count <= exact_limit) {
        result = pruned_closed_form(frontier, frontier_count, reference, dimension);
    } else {
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

        Region root = {frontier, frontier_count, lower, reference, dimension};
        Workspace workspace = {NULL, NULL, 0};
        result = pruned_partition_region(&root, &workspace, 0, exact_limit);
        workspace_destroy(&workspace);
    }

    free(candidates);
    free(frontier);
    return result;
}
