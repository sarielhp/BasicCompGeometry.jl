#!/usr/bin/env julia

using QuickEnv

#
# examples/prune_halfplanes_depth.jl
#
# Prunes a set of halfplanes to a minimal subset whose arrangement preserves
# depth at least k at all cell centers (centroids of the vertical trapezoids).
#
# 1. Takes a parameter k (default k=5).
# 2. Generates 6*k random halfplanes (or reads an initial set from a text file).
# 3. Computes the default view and vertical decomposition inside the view.
# 4. Computes the vector of trapezoid centers and their initial depths.
# 5. Repeatedly attempts to remove halfplanes from the set. A removal is legal if
#    all center depths remain >= k after removal.
# 6. Stops when all remaining halfplanes have been tried in a full pass and none
#    could be legally removed.
# 7. Saves the remaining halfplanes to a text file using write_halfplanes.
# 8. Prints the final size of the set and depth statistics.
#

using BasicCompGeometry
using LinearAlgebra
using Printf
using Random
using Statistics

function run_prune_halfplanes_depth(
    k::Int = 5;
    input_file::Union{String,Nothing} = nothing,
    output_file::Union{String,Nothing} = nothing,
    seed::Int = 42,
    multiplier::Int = 6,
)
    rng = Random.Xoshiro(seed)

    println("="^70)
    @printf("Program: prune_halfplanes_depth\n")
    @printf("Minimal Subset of Halfplanes Preserving Center Depth >= k\n")
    println("="^70)
    @printf("Target minimum depth (k):       %d\n", k)

    # 1. Load initial set or generate 6*k random halfplanes
    hps = Halfplane2F[]
    if input_file !== nothing && isfile(input_file)
        println("\n[1/5] Loading initial halfplanes from file: ", input_file)
        hps = read_halfplanes(input_file)
        @printf("      Loaded %d halfplanes from file.\n", length(hps))
    else
        m = multiplier * k
        @printf("\n[1/5] Generating %d random halfplanes (6*k, seed = %d)...\n", m, seed)
        # Ensure initial min depth >= k by retrying if needed
        for attempt in 1:10
            candidate_hps = [rand_halfplane(rng) for _ = 1:m]
            cand_view = default_view(candidate_hps; factor = 1.1)
            cand_decomp = vertical_decomposition(candidate_hps, cand_view)
            cand_centers = [centroid(t.corners) for t in cand_decomp]
            cand_min_d = minimum(depth(candidate_hps, c) for c in cand_centers)
            if cand_min_d >= k
                hps = candidate_hps
                @printf("      Generated %d halfplanes (initial arrangement min depth = %d >= k).\n", m, cand_min_d)
                break
            elseif attempt == 10
                hps = candidate_hps
                @printf("      Generated %d halfplanes (initial min depth = %d).\n", m, cand_min_d)
            end
        end
    end

    m = length(hps)
    if m == 0
        error("No halfplanes to process.")
    end

    # 2. Compute default view
    println("\n[2/5] Computing default view rectangle...")
    t_view = @elapsed view = default_view(hps; factor = 1.1)
    view_w = width(view)
    view_h = height(view)
    @printf("      View box:                 [%.4f, %.4f] × [%.4f, %.4f]\n",
            view.mini[1], view.maxi[1], view.mini[2], view.maxi[2])
    @printf("      Dimensions:               width = %.4f, height = %.4f, area = %.4f\n",
            view_w, view_h, view_w * view_h)

    # 3. Compute vertical decomposition inside the view
    println("\n[3/5] Computing vertical decomposition inside view...")
    t_decomp = @elapsed decomp = vertical_decomposition(hps, view)
    num_traps = length(decomp)
    num_verts = length(vertices(decomp))
    @printf("      Decomposition time:       %.3f s\n", t_decomp)
    @printf("      Arrangement vertices:     %d\n", num_verts)
    @printf("      Vertical trapezoids:      %d\n", num_traps)

    # 4. Explicitly compute centers and initial depths
    println("\n[4/5] Computing trapezoid centers and incidence matrix...")
    centers = [centroid(trap.corners) for trap in decomp]
    num_centers = length(centers)

    # Precompute incidence: which halfplanes contain which centers
    # inside_matrix[i, j] = is center j inside halfplane i?
    t_inc = @elapsed begin
        inside_matrix = BitMatrix(undef, m, num_centers)
        for i = 1:m
            h = hps[i]
            for j = 1:num_centers
                inside_matrix[i, j] = is_inside(centers[j], h)
            end
        end
        # Map each halfplane to the list of centers it contains
        halfplane_centers = [findall(j -> inside_matrix[i, j], 1:num_centers) for i = 1:m]
        # Current depth at each center
        depths = [count(i -> inside_matrix[i, j], 1:m) for j in 1:num_centers]
    end

    init_min_depth = minimum(depths)
    init_max_depth = maximum(depths)
    init_mean_depth = mean(depths)
    @printf("      Incidence computed in:    %.3f s\n", t_inc)
    @printf("      Initial min depth:        %d\n", init_min_depth)
    @printf("      Initial max depth:        %d\n", init_max_depth)
    @printf("      Initial mean depth:       %.2f\n", init_mean_depth)

    if init_min_depth < k
        @warn @sprintf("Initial arrangement has min depth %d < target k=%d. Pruning will maintain depths >= %d.",
                       init_min_depth, k, init_min_depth)
    end

    # 5. Greedy pruning loop
    println("\n[5/5] Pruning halfplanes while keeping all center depths >= k...")
    active = trues(m)
    pass_num = 0
    total_removed = 0

    t_prune = @elapsed while true
        pass_num += 1
        removed_in_pass = 0

        for i = 1:m
            !active[i] && continue

            # Removal is legal iff for all centers contained in halfplane i,
            # their current depth > k (so after removal depth >= k).
            c_indices = halfplane_centers[i]
            can_remove = true
            for j in c_indices
                if depths[j] <= k
                    can_remove = false
                    break
                end
            end

            if can_remove
                active[i] = false
                for j in c_indices
                    depths[j] -= 1
                end
                removed_in_pass += 1
                total_removed += 1
            end
        end

        @printf("      Pass %d: removed %d halfplanes (%d remaining)\n",
                pass_num, removed_in_pass, count(active))

        # Stop when a full pass over all currently active halfplanes fails to remove any
        if removed_in_pass == 0
            break
        end
    end

    final_hps = hps[active]
    final_size = length(final_hps)
    final_min_d = minimum(depths)
    final_max_d = maximum(depths)
    final_mean_d = mean(depths)

    println("\n" * "="^70)
    @printf("PRUNING COMPLETE (took %.3f s across %d passes)\n", t_prune, pass_num)
    println("="^70)
    @printf("  • Target depth parameter (k):        %d\n", k)
    @printf("  • Initial halfplanes count:           %d\n", m)
    @printf("  • Halfplanes removed:                 %d (%.1f%%)\n",
            total_removed, (total_removed / m) * 100)
    @printf("  • Final pruned set size:              %d\n", final_size)
    @printf("  • Final min depth over all centers:   %d (>= k = %d)\n", final_min_d, k)
    @printf("  • Final max depth over all centers:   %d\n", final_max_d)
    @printf("  • Final mean depth over all centers:  %.2f\n", final_mean_d)
    println("="^70)

    # 6. Save pruned halfplanes to text file
    out_path = output_file !== nothing ? output_file : "output/pruned_halfplanes_k$(k).txt"
    write_halfplanes(out_path, final_hps)
    println("\nSaved pruned set of $(final_size) halfplanes to: ", out_path)

    # Explicit output of the final set size
    @printf("\nFinal size of the halfplane set: %d\n", final_size)

    return final_hps
end

function main()
    k = 5
    input_file = nothing
    output_file = nothing
    seed = 42

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--input" || arg == "-i"
            i += 1
            input_file = ARGS[i]
        elseif arg == "--output" || arg == "-o"
            i += 1
            output_file = ARGS[i]
        elseif arg == "--seed" || arg == "-s"
            i += 1
            seed = parse(Int, ARGS[i])
        elseif !startswith(arg, "-")
            # Try to parse as integer k, or treat as input file if it exists
            if isfile(arg)
                input_file = arg
            else
                try
                    k = parse(Int, arg)
                catch
                    input_file = arg
                end
            end
        end
        i += 1
    end

    run_prune_halfplanes_depth(
        k;
        input_file = input_file,
        output_file = output_file,
        seed = seed,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
