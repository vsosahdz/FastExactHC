/*
 * bitmatrix: pivot partitioning with an inherited comparison matrix W, built once,
 * whose bit masks answer dominance and cover tests in every node and leaf, plus
 * children ordered by size. Entry point: hypervolume_c_imatrix.
 */
#include "hypervolume_c_ie_pruned.c"

#include <string.h>

/*
 * Exact HV with the pivot partition plus four changes:
 *
 * 1. One-pass distribution: a point goes to child j only when p_j < pivot_j.
 * 2. Objective order per node: children are generated in ascending order of
 *    how many points they receive, so the largest children get the most
 *    objectives clamped to the pivot and lose the most dominated points.
 * 3. Inherited comparison matrix: worse[l][a] is the mask of objectives in
 *    which p_l is worse than p_a. It is built once and every descendant
 *    reuses it through point ids. In a node with lower bound L, clamping
 *    gives p_l' <= p_a' in objective o iff p_l <= p_a or p_l <= L_o, so with
 *    low[l] = objectives where p_l was clamped to L:
 *        p_l dominates p_a          iff (worse[l][a] & ~low[l]) == 0
 *        p_l covers join(p_a, p_b)  iff (worse[l][a] & worse[l][b] & ~low[l]) == 0
 *    Dominance filtering and the base-case matrix become bit operations.
 * 4. Intersection matrix in the base case: cell (a, b) is the mask of later
 *    points covering join(p_a, p_b); a subset containing a and b has a larger
 *    join, so its subtree is pruned with bit operations.
 */

#define IMATRIX_DEFAULT_LIMIT 10
#define IMATRIX_MAX_LIMIT 64
#define IMATRIX_MATRIX_CAP 2048

typedef uint64_t CoverMask;
typedef uint64_t ObjectiveMask;

typedef struct {
    uint64_t **buffers;
    size_t *capacities;
    size_t levels;
} MaskWorkspace;

static uint64_t *mask_buffer(MaskWorkspace *workspace, size_t level, size_t count) {
    if (level >= workspace->levels) {
        size_t new_levels = level + 1;
        uint64_t **new_buffers = realloc(
            workspace->buffers,
            new_levels * sizeof(*new_buffers)
        );
        if (new_buffers == NULL) {
            return NULL;
        }
        workspace->buffers = new_buffers;
        size_t *new_capacities = realloc(
            workspace->capacities,
            new_levels * sizeof(*new_capacities)
        );
        if (new_capacities == NULL) {
            return NULL;
        }
        workspace->capacities = new_capacities;
        for (size_t index = workspace->levels; index < new_levels; ++index) {
            workspace->buffers[index] = NULL;
            workspace->capacities[index] = 0;
        }
        workspace->levels = new_levels;
    }
    if (count > workspace->capacities[level]) {
        uint64_t *buffer = realloc(workspace->buffers[level], count * sizeof(*buffer));
        if (buffer == NULL) {
            return NULL;
        }
        workspace->buffers[level] = buffer;
        workspace->capacities[level] = count;
    }
    return workspace->buffers[level];
}

static void mask_workspace_destroy(MaskWorkspace *workspace) {
    for (size_t level = 0; level < workspace->levels; ++level) {
        free(workspace->buffers[level]);
    }
    free(workspace->buffers);
    free(workspace->capacities);
}

typedef struct {
    const ObjectiveMask *worse;
    size_t width;
} WorseMatrix;

static ObjectiveMask worse_mask(const WorseMatrix *matrix, uint64_t left, uint64_t right) {
    return matrix->worse[left * matrix->width + right];
}

typedef struct {
    const double *points;
    size_t count;
    const double *reference;
    size_t dimension;
    const CoverMask *pair_cover;
    size_t *members;
} MatrixLeaf;

static double matrix_subtree(
    const MatrixLeaf *leaf,
    const double *join,
    CoverMask cover,
    size_t member_count,
    size_t next,
    double sign
) {
    CoverMask later = next >= IMATRIX_MAX_LIMIT ? 0 : ~(CoverMask)0 << next;
    if ((cover & later) != 0) {
        return 0.0;
    }

    size_t dimension = leaf->dimension;
    double volume = sign * box_volume(join, leaf->reference, dimension);
    double child_join[dimension];
    for (size_t point = next; point < leaf->count; ++point) {
        CoverMask child_cover = cover;
        for (size_t member = 0; member < member_count; ++member) {
            child_cover |= leaf->pair_cover[
                leaf->members[member] * leaf->count + point
            ];
        }
        for (size_t objective = 0; objective < dimension; ++objective) {
            double value = leaf->points[point * dimension + objective];
            child_join[objective] = value > join[objective] ? value : join[objective];
        }
        leaf->members[member_count] = point;
        volume += matrix_subtree(
            leaf,
            child_join,
            child_cover,
            member_count + 1,
            point + 1,
            -sign
        );
    }
    return volume;
}

static void swap_rows(
    double *points,
    uint64_t *ids,
    ObjectiveMask *low,
    size_t first,
    size_t second,
    size_t dimension
) {
    for (size_t objective = 0; objective < dimension; ++objective) {
        double value = points[first * dimension + objective];
        points[first * dimension + objective] = points[second * dimension + objective];
        points[second * dimension + objective] = value;
    }
    uint64_t id = ids[first];
    ids[first] = ids[second];
    ids[second] = id;
    ObjectiveMask mask = low[first];
    low[first] = low[second];
    low[second] = mask;
}

/* count must not exceed IMATRIX_MAX_LIMIT. Reorders points, ids and low. */
static double matrix_closed_form(
    double *points,
    uint64_t *ids,
    ObjectiveMask *low,
    size_t count,
    const double *reference,
    size_t dimension,
    const WorseMatrix *matrix
) {
    double volumes[count];
    for (size_t point = 0; point < count; ++point) {
        volumes[point] = box_volume(points + point * dimension, reference, dimension);
    }
    for (size_t point = 1; point < count; ++point) {
        for (size_t index = point; index > 0 && volumes[index - 1] > volumes[index]; --index) {
            double volume = volumes[index];
            volumes[index] = volumes[index - 1];
            volumes[index - 1] = volume;
            swap_rows(points, ids, low, index, index - 1, dimension);
        }
    }

    CoverMask pair_cover[count * count];
    for (size_t first = 0; first < count; ++first) {
        for (size_t second = first + 1; second < count; ++second) {
            CoverMask mask = 0;
            for (size_t point = second + 1; point < count; ++point) {
                ObjectiveMask failing =
                    worse_mask(matrix, ids[point], ids[first]) &
                    worse_mask(matrix, ids[point], ids[second]) &
                    ~low[point];
                if (failing == 0) {
                    mask |= (CoverMask)1 << point;
                }
            }
            pair_cover[first * count + second] = mask;
        }
    }

    size_t members[count];
    MatrixLeaf leaf = {points, count, reference, dimension, pair_cover, members};
    double volume = 0.0;
    for (size_t point = 0; point < count; ++point) {
        members[0] = point;
        volume += matrix_subtree(
            &leaf,
            points + point * dimension,
            0,
            1,
            point + 1,
            1.0
        );
    }
    return volume;
}

/* In-place nondominated filter using the inherited matrix. */
static size_t matrix_nondominated(
    double *points,
    uint64_t *ids,
    ObjectiveMask *low,
    size_t count,
    size_t dimension,
    const WorseMatrix *matrix
) {
    size_t kept = 0;
    for (size_t candidate = 0; candidate < count; ++candidate) {
        int is_dominated = 0;
        for (size_t index = 0; index < kept; ++index) {
            if ((worse_mask(matrix, ids[index], ids[candidate]) & ~low[index]) == 0) {
                is_dominated = 1;
                break;
            }
        }
        if (is_dominated) {
            continue;
        }

        size_t index = 0;
        while (index < kept) {
            if ((worse_mask(matrix, ids[candidate], ids[index]) & ~low[candidate]) == 0) {
                --kept;
                swap_rows(points, ids, low, index, kept, dimension);
                continue;
            }
            ++index;
        }
        if (kept != candidate) {
            swap_rows(points, ids, low, kept, candidate, dimension);
        }
        ++kept;
    }
    return kept;
}

typedef struct {
    Workspace points;
    Workspace scratch;
    MaskWorkspace ids;
    MaskWorkspace low;
    size_t exact_limit;
} MatrixWorkspace;

/* points are already inside [lower, reference) and clamped to lower. */
static double imatrix_partition_region(
    double *points,
    uint64_t *ids,
    size_t count,
    const double *lower,
    const double *reference,
    size_t dimension,
    MatrixWorkspace *workspace,
    size_t level,
    const WorseMatrix *matrix
) {
    if (count == 0) {
        return 0.0;
    }

    ObjectiveMask *owned = NULL;
    WorseMatrix local_matrix;
    if (matrix == NULL && count <= IMATRIX_MATRIX_CAP) {
        owned = malloc(count * count * sizeof(*owned));
        if (owned != NULL) {
            for (size_t left = 0; left < count; ++left) {
                ids[left] = left;
                for (size_t right = 0; right < count; ++right) {
                    ObjectiveMask mask = 0;
                    for (size_t objective = 0; objective < dimension; ++objective) {
                        if (points[left * dimension + objective] >
                                points[right * dimension + objective]) {
                            mask |= (ObjectiveMask)1 << objective;
                        }
                    }
                    owned[left * count + right] = mask;
                }
            }
            local_matrix.worse = owned;
            local_matrix.width = count;
            matrix = &local_matrix;
        }
    }

    ObjectiveMask *low = NULL;
    if (matrix != NULL) {
        low = mask_buffer(&workspace->low, level, count);
        if (low == NULL) {
            free(owned);
            return -1.0;
        }
        for (size_t point = 0; point < count; ++point) {
            ObjectiveMask mask = 0;
            for (size_t objective = 0; objective < dimension; ++objective) {
                if (points[point * dimension + objective] <= lower[objective]) {
                    mask |= (ObjectiveMask)1 << objective;
                }
            }
            low[point] = mask;
        }
    }

    double volume = 0.0;
    if (count > workspace->exact_limit) {
        if (matrix != NULL) {
            count = matrix_nondominated(points, ids, low, count, dimension, matrix);
        } else {
            double *frontier = workspace_buffer(&workspace->scratch, level, count, dimension);
            if (frontier == NULL) {
                free(owned);
                return -1.0;
            }
            count = nondominated(points, count, dimension, frontier);
            memcpy(points, frontier, count * dimension * sizeof(double));
        }
    }
    if (count == 1) {
        volume = box_volume(points, reference, dimension);
        free(owned);
        return volume;
    }
    if (count <= workspace->exact_limit) {
        if (matrix != NULL) {
            volume = matrix_closed_form(points, ids, low, count, reference, dimension, matrix);
        } else {
            volume = pruned_closed_form(points, count, reference, dimension);
        }
        free(owned);
        return volume;
    }

    size_t pivot_index = 0;
    double pivot_volume = box_volume(points, reference, dimension);
    for (size_t point = 1; point < count; ++point) {
        double value = box_volume(points + point * dimension, reference, dimension);
        if (value > pivot_volume) {
            pivot_index = point;
            pivot_volume = value;
        }
    }
    double pivot[dimension];
    memcpy(pivot, points + pivot_index * dimension, dimension * sizeof(double));

    size_t child_count[dimension];
    size_t order[dimension];
    for (size_t objective = 0; objective < dimension; ++objective) {
        child_count[objective] = 0;
        order[objective] = objective;
    }
    for (size_t point = 0; point < count; ++point) {
        for (size_t objective = 0; objective < dimension; ++objective) {
            if (points[point * dimension + objective] < pivot[objective]) {
                ++child_count[objective];
            }
        }
    }
    for (size_t position = 1; position < dimension; ++position) {
        size_t objective = order[position];
        size_t index = position;
        while (index > 0 && child_count[order[index - 1]] > child_count[objective]) {
            order[index] = order[index - 1];
            --index;
        }
        order[index] = objective;
    }

    volume = pivot_volume;
    double child_lower[dimension];
    double child_reference[dimension];
    for (size_t objective = 0; objective < dimension; ++objective) {
        child_lower[objective] = lower[objective];
        child_reference[objective] = reference[objective];
    }

    for (size_t position = 0; position < dimension; ++position) {
        size_t split = order[position];
        if (position > 0) {
            size_t previous = order[position - 1];
            child_reference[previous] = reference[previous];
            if (child_lower[previous] < pivot[previous]) {
                child_lower[previous] = pivot[previous];
            }
        }
        child_reference[split] = pivot[split];
        if (child_count[split] == 0) {
            continue;
        }

        double *child_points = workspace_buffer(
            &workspace->points,
            level + 1,
            child_count[split],
            dimension
        );
        uint64_t *child_ids = mask_buffer(&workspace->ids, level + 1, child_count[split]);
        if (child_points == NULL || child_ids == NULL) {
            free(owned);
            return -1.0;
        }
        size_t local_count = 0;
        for (size_t point = 0; point < count; ++point) {
            const double *source = points + point * dimension;
            if (!(source[split] < pivot[split])) {
                continue;
            }
            for (size_t objective = 0; objective < dimension; ++objective) {
                child_points[local_count * dimension + objective] =
                    source[objective] > child_lower[objective]
                    ? source[objective]
                    : child_lower[objective];
            }
            child_ids[local_count] = ids[point];
            ++local_count;
        }

        double child_volume = imatrix_partition_region(
            child_points,
            child_ids,
            local_count,
            child_lower,
            child_reference,
            dimension,
            workspace,
            level + 1,
            matrix
        );
        if (child_volume < 0.0) {
            free(owned);
            return -1.0;
        }
        volume += child_volume;
    }

    free(owned);
    return volume;
}

/*
 * Same contract as hypervolume_c. exact_limit = 0 uses IMATRIX_DEFAULT_LIMIT;
 * values above IMATRIX_MAX_LIMIT are clamped because leaf masks use 64 bits.
 * More than 64 objectives falls back to hypervolume_c_ie_pruned.
 */
double hypervolume_c_imatrix(
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
        exact_limit = IMATRIX_DEFAULT_LIMIT;
    }
    if (exact_limit > IMATRIX_MAX_LIMIT) {
        exact_limit = IMATRIX_MAX_LIMIT;
    }
    if (dimension > sizeof(ObjectiveMask) * CHAR_BIT) {
        return hypervolume_c_ie_pruned(input_points, count, dimension, reference, exact_limit);
    }

    MatrixWorkspace workspace = {
        {NULL, NULL, 0},
        {NULL, NULL, 0},
        {NULL, NULL, 0},
        {NULL, NULL, 0},
        exact_limit,
    };
    double *frontier = workspace_buffer(&workspace.points, 0, count, dimension);
    uint64_t *ids = mask_buffer(&workspace.ids, 0, count);
    double result = -1.0;
    if (frontier != NULL && ids != NULL) {
        size_t frontier_count = nondominated(input_points, count, dimension, frontier);
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
        result = imatrix_partition_region(
            frontier,
            ids,
            frontier_count,
            lower,
            reference,
            dimension,
            &workspace,
            0,
            NULL
        );
    }
    workspace_destroy(&workspace.points);
    workspace_destroy(&workspace.scratch);
    mask_workspace_destroy(&workspace.ids);
    mask_workspace_destroy(&workspace.low);
    return result;
}
