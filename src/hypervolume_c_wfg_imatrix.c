/*
 * wfg-bitmatrix: WFG recursion at high dimension with limit sets filtered by the
 * comparison matrix W, switching to the bitmatrix engine at switch_dimension or fewer
 * objectives. Entry point: hypervolume_c_wfg_imatrix.
 */
#include "hypervolume_c_imatrix.c"

/*
 * Hybrid exact HV: WFG recursion at high dimension, imatrix partition at low
 * dimension, both driven by one comparison matrix built at the root.
 *
 * WFG step. Points sorted by the last active objective, worst first:
 *     HV = sum_i (ref_last - q_i,last) * (vol(q_i) - HV(limit_i))
 * over the remaining objectives, with limit_i = nd{ max(q_i, q_j) : j > i }.
 *
 * Every point of a limit set has the form max(C, p_j), where C is the context
 * (the max of the points chosen in the enclosing levels) and p_j an original
 * point. With low_C[j] = objectives where p_j <= C:
 *     max(C, p_j) <= max(C, p_k)  iff  (worse[j][k] & ~low_C[j] & active) == 0
 * which is the same test imatrix uses after clamping to a lower bound. Limit
 * sets are filtered with bit operations on the root matrix, and the imatrix
 * engine reuses that matrix through an objective map once the dimension is
 * small enough.
 */

#define WFG_IMATRIX_DEFAULT_SWITCH 12

typedef struct {
    double key;
    uint64_t id;
} WfgOrder;

static int wfg_order_descending(const void *left, const void *right) {
    double a = ((const WfgOrder *)left)->key;
    double b = ((const WfgOrder *)right)->key;
    return (a < b) - (a > b);
}

/* imatrix_partition_region over the objectives listed in objective_map. */
static double mapped_partition_region(
    double *points,
    uint64_t *ids,
    size_t count,
    const double *lower,
    const double *reference,
    size_t dimension,
    const size_t *objective_map,
    ObjectiveMask inactive,
    MatrixWorkspace *workspace,
    size_t level,
    const WorseMatrix *matrix
) {
    if (count == 0) {
        return 0.0;
    }

    ObjectiveMask *low = mask_buffer(&workspace->low, level, count);
    if (low == NULL) {
        return -1.0;
    }
    for (size_t point = 0; point < count; ++point) {
        ObjectiveMask mask = inactive;
        for (size_t objective = 0; objective < dimension; ++objective) {
            if (points[point * dimension + objective] <= lower[objective]) {
                mask |= (ObjectiveMask)1 << objective_map[objective];
            }
        }
        low[point] = mask;
    }

    if (count > workspace->exact_limit) {
        count = matrix_nondominated(points, ids, low, count, dimension, matrix);
    }
    if (count == 1) {
        return box_volume(points, reference, dimension);
    }
    if (count <= workspace->exact_limit) {
        return matrix_closed_form(points, ids, low, count, reference, dimension, matrix);
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

    double volume = pivot_volume;
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

        double child_volume = mapped_partition_region(
            child_points,
            child_ids,
            local_count,
            child_lower,
            child_reference,
            dimension,
            objective_map,
            inactive,
            workspace,
            level + 1,
            matrix
        );
        if (child_volume < 0.0) {
            return -1.0;
        }
        volume += child_volume;
    }

    return volume;
}

typedef struct {
    const double *points;
    size_t dimension;
    const double *reference;
    const WorseMatrix *matrix;
    size_t switch_dimension;
    MatrixWorkspace engine;
    Workspace level_points;
    MaskWorkspace level_ids;
    MaskWorkspace level_order;
    MaskWorkspace level_low;
} WfgContext;

/*
 * HV of { max(context, p_id) : id in ids } over the active objectives.
 * active lists original objective indices; ids is not modified.
 */
static double wfg_imatrix_level(
    WfgContext *context_data,
    const uint64_t *ids,
    size_t count,
    const double *context,
    const size_t *active,
    size_t active_count,
    size_t level
) {
    if (count == 0) {
        return 0.0;
    }
    size_t dimension = context_data->dimension;
    const double *points = context_data->points;
    const double *reference = context_data->reference;

    double active_reference[active_count];
    double active_lower[active_count];
    ObjectiveMask active_mask = 0;
    for (size_t objective = 0; objective < active_count; ++objective) {
        active_reference[objective] = reference[active[objective]];
        active_lower[objective] = context[active[objective]];
        active_mask |= (ObjectiveMask)1 << active[objective];
    }

    WfgOrder *order = (WfgOrder *)mask_buffer(
        &context_data->level_order,
        level,
        2 * count
    );
    double *values = workspace_buffer(
        &context_data->level_points,
        level,
        count,
        active_count
    );
    if (order == NULL || values == NULL) {
        return -1.0;
    }
    size_t last = active[active_count - 1];
    for (size_t point = 0; point < count; ++point) {
        const double *source = points + ids[point] * dimension;
        order[point].key = source[last] > context[last] ? source[last] : context[last];
        order[point].id = ids[point];
    }
    if (count == 1) {
        const double *source = points + order[0].id * dimension;
        double volume = 1.0;
        for (size_t objective = 0; objective < active_count; ++objective) {
            size_t original = active[objective];
            double value = source[original] > context[original]
                ? source[original]
                : context[original];
            volume *= reference[original] - value;
        }
        return volume;
    }

    if (active_count <= context_data->switch_dimension ||
            count <= context_data->engine.exact_limit) {
        uint64_t *engine_ids = mask_buffer(&context_data->level_ids, level, count);
        if (engine_ids == NULL) {
            return -1.0;
        }
        for (size_t point = 0; point < count; ++point) {
            const double *source = points + ids[point] * dimension;
            for (size_t objective = 0; objective < active_count; ++objective) {
                size_t original = active[objective];
                values[point * active_count + objective] =
                    source[original] > context[original]
                    ? source[original]
                    : context[original];
            }
            engine_ids[point] = ids[point];
        }
        return mapped_partition_region(
            values,
            engine_ids,
            count,
            active_lower,
            active_reference,
            active_count,
            active,
            ~active_mask,
            &context_data->engine,
            0,
            context_data->matrix
        );
    }

    qsort(order, count, sizeof(*order), wfg_order_descending);

    size_t sliced_count = active_count - 1;
    ObjectiveMask sliced_mask = active_mask & ~((ObjectiveMask)1 << last);
    uint64_t *limit_ids = mask_buffer(&context_data->level_ids, level, count);
    ObjectiveMask *limit_low = mask_buffer(&context_data->level_low, level, count);
    if (limit_ids == NULL || limit_low == NULL) {
        return -1.0;
    }

    double total = 0.0;
    double child_context[dimension];
    for (size_t index = 0; index < count; ++index) {
        const double *chosen = points + order[index].id * dimension;
        double inclusive = 1.0;
        for (size_t objective = 0; objective < dimension; ++objective) {
            child_context[objective] = chosen[objective] > context[objective]
                ? chosen[objective]
                : context[objective];
        }
        for (size_t objective = 0; objective < sliced_count; ++objective) {
            size_t original = active[objective];
            inclusive *= reference[original] - child_context[original];
        }

        size_t kept = 0;
        int covered = 0;
        for (size_t later = index + 1; later < count && !covered; ++later) {
            uint64_t candidate = order[later].id;
            const double *source = points + candidate * dimension;
            ObjectiveMask low = ~sliced_mask;
            for (size_t objective = 0; objective < sliced_count; ++objective) {
                size_t original = active[objective];
                if (source[original] <= child_context[original]) {
                    low |= (ObjectiveMask)1 << original;
                }
            }
            if (low == ~(ObjectiveMask)0) {
                covered = 1;
                break;
            }

            int is_dominated = 0;
            for (size_t other = 0; other < kept; ++other) {
                if ((worse_mask(context_data->matrix, limit_ids[other], candidate) &
                        ~limit_low[other]) == 0) {
                    is_dominated = 1;
                    break;
                }
            }
            if (is_dominated) {
                continue;
            }
            size_t other = 0;
            while (other < kept) {
                if ((worse_mask(context_data->matrix, candidate, limit_ids[other]) &
                        ~low) == 0) {
                    --kept;
                    limit_ids[other] = limit_ids[kept];
                    limit_low[other] = limit_low[kept];
                    continue;
                }
                ++other;
            }
            limit_ids[kept] = candidate;
            limit_low[kept] = low;
            ++kept;
        }

        double exclusive = 0.0;
        if (!covered) {
            double limit_volume = wfg_imatrix_level(
                context_data,
                limit_ids,
                kept,
                child_context,
                active,
                sliced_count,
                level + 1
            );
            if (limit_volume < 0.0) {
                return -1.0;
            }
            exclusive = inclusive - limit_volume;
        }
        total += (reference[last] - order[index].key) * exclusive;
    }
    return total;
}

/*
 * Same contract as hypervolume_c. exact_limit = 0 uses IMATRIX_DEFAULT_LIMIT;
 * switch_dimension = 0 uses WFG_IMATRIX_DEFAULT_SWITCH: subproblems with that
 * many objectives or fewer go to the imatrix engine. Fronts larger than
 * IMATRIX_MATRIX_CAP or more than 64 objectives fall back to imatrix.
 */
double hypervolume_c_wfg_imatrix(
    const double *input_points,
    size_t count,
    size_t dimension,
    const double *reference,
    size_t exact_limit,
    size_t switch_dimension
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
    if (switch_dimension == 0) {
        switch_dimension = WFG_IMATRIX_DEFAULT_SWITCH;
    }
    if (switch_dimension < 2) {
        switch_dimension = 2;
    }
    if (dimension <= switch_dimension || dimension > sizeof(ObjectiveMask) * CHAR_BIT) {
        return hypervolume_c_imatrix(input_points, count, dimension, reference, exact_limit);
    }

    double *frontier = malloc(count * dimension * sizeof(double));
    if (frontier == NULL) {
        return -1.0;
    }
    size_t frontier_count = nondominated(input_points, count, dimension, frontier);
    if (frontier_count > IMATRIX_MATRIX_CAP) {
        free(frontier);
        return hypervolume_c_imatrix(input_points, count, dimension, reference, exact_limit);
    }

    ObjectiveMask *worse = malloc(frontier_count * frontier_count * sizeof(*worse));
    uint64_t *ids = malloc(frontier_count * sizeof(*ids));
    if (worse == NULL || ids == NULL) {
        free(worse);
        free(ids);
        free(frontier);
        return -1.0;
    }
    for (size_t left = 0; left < frontier_count; ++left) {
        ids[left] = left;
        for (size_t right = 0; right < frontier_count; ++right) {
            ObjectiveMask mask = 0;
            for (size_t objective = 0; objective < dimension; ++objective) {
                if (frontier[left * dimension + objective] >
                        frontier[right * dimension + objective]) {
                    mask |= (ObjectiveMask)1 << objective;
                }
            }
            worse[left * frontier_count + right] = mask;
        }
    }
    WorseMatrix matrix = {worse, frontier_count};

    WfgContext context_data = {
        frontier,
        dimension,
        reference,
        &matrix,
        switch_dimension,
        {{NULL, NULL, 0}, {NULL, NULL, 0}, {NULL, NULL, 0}, {NULL, NULL, 0}, exact_limit},
        {NULL, NULL, 0},
        {NULL, NULL, 0},
        {NULL, NULL, 0},
        {NULL, NULL, 0},
    };
    double context[dimension];
    size_t active[dimension];
    for (size_t objective = 0; objective < dimension; ++objective) {
        context[objective] = -DBL_MAX;
        active[objective] = objective;
    }

    double result = wfg_imatrix_level(
        &context_data,
        ids,
        frontier_count,
        context,
        active,
        dimension,
        0
    );

    workspace_destroy(&context_data.engine.points);
    workspace_destroy(&context_data.engine.scratch);
    mask_workspace_destroy(&context_data.engine.ids);
    mask_workspace_destroy(&context_data.engine.low);
    workspace_destroy(&context_data.level_points);
    mask_workspace_destroy(&context_data.level_ids);
    mask_workspace_destroy(&context_data.level_order);
    mask_workspace_destroy(&context_data.level_low);
    free(worse);
    free(ids);
    free(frontier);
    return result;
}
