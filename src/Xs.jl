module Xs
using StaticArrays, Distributed, ProgressBars, Printf
using ..GEN

# ---------------------------------------------------------------------------
# Auxiliary functions: bin values and range
# ---------------------------------------------------------------------------
function binx(i::Int64, bin, iaxis::Int64)::Float64
    return bin.min[iaxis] + (i - 0.5) / bin.Nbin * (bin.max[iaxis] - bin.min[iaxis])
end

# Direct calculation, axis can be specified by labels, e.g., "p1:p2"
function binrange(axis::String, tecm, proc, stype)
    laxes = [findfirst(==(element), proc.pf) for element in split(axis, ":")]
    msum = sum(proc.mf[laxes])
    min = msum^stype
    max = (tecm - sum(proc.mf) + msum)^stype
    return min, max
end

# Single combination -> (min, max) where axis refers to the indices of final-state particles
function binrange(axis::Vector{Int64}, tecm, proc, stype)
    msum = sum(proc.mf[axis])
    return msum^stype, (tecm - sum(proc.mf) + msum)^stype
end

# Multiple combinations -> the smallest min and the largest max among all combinations
function binrange(laxes::Vector{Vector{Int64}}, tecm, proc, stype)
    mins = Float64[]
    maxs = Float64[]
    for axis in laxes
        mn, mx = binrange(axis, tecm, proc, stype)
        push!(mins, mn)
        push!(maxs, mx)
    end
    return minimum(mins), maximum(maxs)
end

# Multi-axis version -> one scalar (min, max) per axis, returning two vectors
function binrange(laxes::Vector{Vector{Vector{Int64}}}, tecm, proc, stype; Range=[])
    min_all = Float64[]
    max_all = Float64[]
    for i in eachindex(laxes)
        min, max = binrange(laxes[i], tecm, proc, stype)   # 标量版本
        push!(min_all, min)
        push!(max_all, max)
    end
    # Preset values
    if !isempty(Range)
        for i in eachindex(min_all)
            if !isempty(Range[i])
                !ismissing(Range[i][1]) && (min_all[i] = Range[i][1])
                !ismissing(Range[i][2]) && (max_all[i] = Range[i][2])
            end
        end
    end
    return min_all, max_all
end

# bin index
function Nsij(kijs, min::Float64, max::Float64, Nbin)
    return convert(Int64, cld((kijs - min) * Nbin, (max - min)))
end

# Retrieve s value by index (single axis)
function Nsum3(laxes::Vector{Vector{Int64}}, i, bin::NamedTuple, kf, stype)
    Ns = Int64[]
    sij = Float64[]
    for combo in laxes
        kij = sum(kf[combo])
        kijs = stype == 2 ? kij * kij : sqrt(kij * kij)
        push!(sij, kijs)
        push!(Ns, Nsij(kijs, bin.min[i], bin.max[i], bin.Nbin))
    end
    return Ns, sij
end

# Multi-axis 
function Nsum3(laxes::Vector{Vector{Vector{Int64}}}, bin::NamedTuple, kf, stype)
    Ns = Vector{Vector{Int64}}(undef, length(laxes))
    sij = Vector{Vector{Float64}}(undef, length(laxes))
    for i in eachindex(laxes)
        Ns[i], sij[i] = Nsum3(laxes[i], i, bin, kf, stype)
    end
    return Ns, sij
end

# ---------------------------------------------------------------------------
# Kinematics auxiliary functions
# ---------------------------------------------------------------------------
function plab2pcm(p::Float64, mi::Vector{Float64})
    ebm = sqrt(p^2 + mi[1]^2)
    return sqrt((ebm + mi[2])^2 - p^2)
end

function getkf(p::Float64, kf::Vector{SVector{5,Float64}}, proc::NamedTuple)
    ebm = sqrt(p^2 + proc.mi[1]^2)
    tecm = sqrt((ebm + proc.mi[2])^2 - p^2)
    GAM = (ebm + proc.mi[2]) / tecm
    ETA = p / tecm
    return [SVector{5,Float64}(
        GAM * k[1] + ETA * k[4],
        k[2], k[3],
        ETA * k[1] + GAM * k[4],
        k[5]
    ) for k in kf]
end

function plab(p::Float64, mi::Vector{Float64})
    m1, m2 = mi[1], mi[2]
    p1 = SVector{5,Float64}([p, 0.0, 0.0, sqrt(p^2 + m1^2), m1])
    p2 = SVector{5,Float64}([0.0, 0.0, 0.0, m2, m2])
    return p1, p2
end

function pcm(tecm::Float64, mi::Vector{Float64})
    m1, m2 = mi[1], mi[2]
    E1 = (tecm^2 + m1^2 - m2^2) / (2.0 * tecm)
    E2 = (tecm^2 + m2^2 - m1^2) / (2.0 * tecm)
    p = sqrt(E2^2 - m2^2)
    p1 = SVector{5,Float64}([p, 0.0, 0.0, E1, m1])
    p2 = SVector{5,Float64}([-p, 0.0, 0.0, E2, m2])
    return p1, p2
end

# ---------------------------------------------------------------------------
# main function for cross section and Dalitz plot
# ---------------------------------------------------------------------------
function Xsection(tecm, proc, callback; axes=[], Range=[], nevtot=Int64(1e6),
    Nbin=1000, para=(l=1.0), p0=[], stype=1, fixed=true, symmetrize=true)

    # ---------- 1. Generate all unique particle index combinations ----------
    laxes_full = Vector{Vector{Int64}}[]
    for axis in axes
        particles = split(axis, ":")
        positions = [findall(==(p), proc.pf) for p in particles]
        combos = Vector{Int64}[]
        for inds in Iterators.product(positions...)
            if allunique(inds)
                push!(combos, collect(inds))
            end
        end
        if length(unique(particles)) == 1
            seen = Set{Vector{Int64}}()
            unique_combos = Vector{Int64}[]
            for c in combos
                sc = sort(c)
                if !(sc in seen)
                    push!(seen, sc)
                    push!(unique_combos, sc)
                end
            end
            combos = unique_combos
        end
        push!(laxes_full, combos)
    end

    Naxes = length(axes)

    # ---------- 2. Build all valid Dalitz pairs ----------
    fill_pairs = Tuple{Vector{Int64},Vector{Int64}}[]
    if Naxes >= 2
        function is_dup(p)
            s = sort([p[1], p[2]])
            return any(x -> sort([x[1], x[2]]) == s, fill_pairs)
        end

        for c1 in laxes_full[1], c2 in laxes_full[2]
            if isempty(intersect(c1, c2)) && !is_dup((c1, c2))
                push!(fill_pairs, (c1, c2))
            end
        end

        if isempty(fill_pairs)
            for c1 in laxes_full[1], c2 in laxes_full[2]
                if length(intersect(c1, c2)) == 1 && !is_dup((c1, c2))
                    push!(fill_pairs, (c1, c2))
                end
            end
        end

        if isempty(fill_pairs)
            @warn "No standard Dalitz pairs found; using ALL combinations."
            for c1 in laxes_full[1], c2 in laxes_full[2]
                if !is_dup((c1, c2))
                    push!(fill_pairs, (c1, c2))
                end
            end
        end
    end

    # ---------- 3. Determine the combination list used for 1D projections ----------
    if symmetrize
        laxes = laxes_full
    else
        laxes = Vector{Vector{Int64}}[]
        if Naxes >= 2
            push!(laxes, [fill_pairs[1][1]])
            push!(laxes, [fill_pairs[1][2]])
            for i in 3:Naxes
                push!(laxes, [laxes_full[i][1]])
            end
        elseif Naxes == 1
            push!(laxes, [laxes_full[1][1]])
        end
    end

    # ---------- 4. Common bin ranges ----------
    Nf = length(proc.pf)
    axesV = []
    min_vals, max_vals = Float64[], Float64[]
    bin = (Nbin=Nbin, min=min_vals, max=max_vals)

    if Nf > 2
        min_vals, max_vals = binrange(laxes_full, tecm, proc, stype, Range=Range)
        bin = (Nbin=Nbin, min=min_vals, max=max_vals)
        axesV = [[binx(ix, bin, iaxis) for ix in 1:Nbin] for iaxis in eachindex(laxes_full)]
    end

    # ---------- 5. Initialize accumulators ----------
    kf, wt = GENEV(tecm, proc.mf)
    amp_init = proc.amps(tecm, kf, proc, para, p0)
    leng = length(amp_init)

    zsum = zeros(Float64, leng)


    dim_axes = Nf > 2 ? Naxes : 1
    dim_bin = Nf > 2 ? Nbin : 1

    zsumt_3d = zeros(Float64, leng, dim_axes, dim_bin)
    zsumd_3d = zeros(Float64, leng, dim_bin, dim_bin)


    wtamp_buffer = zeros(Float64, leng)


    pair_indices = Vector{Tuple{Int,Int}}()
    if Nf > 2 && Naxes >= 2
        for (c1, c2) in fill_pairs
            idx1 = findfirst(==(c1), laxes[1])
            idx2 = findfirst(==(c2), laxes[2])
            if idx1 !== nothing && idx2 !== nothing
                push!(pair_indices, (idx1, idx2))
            end
        end
    end

    # ---------- 6. Event loop ----------
    for ine in 1:nevtot
        kf, wt = GENEV(tecm, proc.mf, fixed=fixed)

        if Nf > 2
            Nsij, sij = Nsum3(laxes, bin, kf, stype)


            is_in_range = true
            @inbounds for i in 1:Naxes
                axis_ok = false
                for val in sij[i]
                    if min_vals[i] <= val <= max_vals[i]
                        axis_ok = true
                        break
                    end
                end
                if !axis_ok
                    is_in_range = false
                    break
                end
            end

            if is_in_range
                amp0 = proc.amps(tecm, kf, proc, para, p0)


                @inbounds for k in 1:leng
                    wtamp_buffer[k] = wt * amp0[k]
                    zsum[k] += wtamp_buffer[k]
                end

                # ---- 1D filling  ----
                for iaxis in 1:Naxes
                    if symmetrize
                        n_comb = length(Nsij[iaxis])
                        inv_n_comb = 1.0 / n_comb
                        for isij in Nsij[iaxis]
                            if 1 < isij <= Nbin

                                @inbounds for k in 1:leng
                                    zsumt_3d[k, iaxis, isij] += wtamp_buffer[k] * inv_n_comb
                                end
                            end
                        end
                    else
                        isij = Nsij[iaxis][1]
                        if 1 < isij <= Nbin
                            @inbounds for k in 1:leng
                                zsumt_3d[k, iaxis, isij] += wtamp_buffer[k]
                            end
                        end
                    end
                end

                # ---- 2D filling  ----
                if Naxes >= 2
                    if symmetrize
                        n_fill = length(pair_indices)
                        inv_n_fill = 1.0 / n_fill
                        for idx in 1:n_fill
                            @inbounds idx1, idx2 = pair_indices[idx]
                            @inbounds isij = Nsij[1][idx1]
                            @inbounds jsij = Nsij[2][idx2]
                            if 1 < isij <= Nbin && 1 < jsij <= Nbin

                                @inbounds for k in 1:leng
                                    zsumd_3d[k, isij, jsij] += wtamp_buffer[k] * inv_n_fill
                                end
                            end
                        end
                    else
                        if !isempty(pair_indices)
                            @inbounds idx1, idx2 = pair_indices[1]
                            @inbounds isij = Nsij[1][idx1]
                            @inbounds jsij = Nsij[2][idx2]
                            if 1 < isij <= Nbin && 1 < jsij <= Nbin
                                @inbounds for k in 1:leng
                                    zsumd_3d[k, isij, jsij] += wtamp_buffer[k]
                                end
                            end
                        end
                    end
                end
            end
        elseif Nf == 2
            amp0 = proc.amps(tecm, kf, proc, para, p0)
            @inbounds for k in 1:leng
                zsum[k] += wt * amp0[k]
            end
        end

        callback(ine)

    end

    # ---------- 7. Normalization ----------
    cs0 = zsum / nevtot
    cs1 = []
    cs2 = []

    if Nf > 2
        binwidths = [(maximum(axesV[i]) - minimum(axesV[i])) / Nbin for i in 1:Naxes]

        if Naxes > 1
            cs2 = [Matrix{Float64}(undef, Nbin, Nbin) for _ in 1:leng]
            for k in 1:leng
                for j in 1:Nbin, i in 1:Nbin
                    cs2[k][i, j] = zsumd_3d[k, i, j] / (nevtot * binwidths[1] * binwidths[2])
                end
            end
        end


        cs1 = [Matrix{Float64}(undef, Naxes, Nbin) for _ in 1:leng]
        for k in 1:leng
            for j in 1:Nbin, i in 1:Naxes
                cs1[k][i, j] = zsumt_3d[k, i, j] / (nevtot * binwidths[i])
            end
        end
    end

    return (cs0=cs0, cs1=cs1, cs2=cs2, axesV=axesV, laxes=laxes, proc=proc, stype=stype)
end

# Parallel worker
function worker_Xsection(tecm, proc, axes, Range, nevt, Nbin, para, p0, stype, progressbar, fixed, symmetrize)
    if myid() == 2 && progressbar
        pb = ProgressBar(1:nevt)
        callback = i -> ProgressBars.update(pb)
    else
        callback = _ -> nothing
    end
    return Xsection(tecm, proc, callback, axes=axes, Range=Range, nevtot=nevt, Nbin=Nbin,
        para=para, p0=p0, stype=stype, fixed=fixed, symmetrize=symmetrize)
end

# Parallel entry 
function Xsection(tecm, proc; axes=[], Range=[], nevtot=Int64(1e6), Nbin=100,
    para=(l=1.0), p0=[], stype=1, progressbar=true, fixed=true, symmetrize=true)
    if fixed
        GEN.reset_genev_rngs!()
    end

    num_workers = nworkers()
    nevt_per_worker = div(nevtot, num_workers)
    ranges = [(i * nevt_per_worker + 1, min((i + 1) * nevt_per_worker, nevtot)) for i in 0:(num_workers-1)]
    GC.gc(false)

    results = pmap(r -> worker_Xsection(tecm, proc, axes, Range, r[2] - r[1] + 1, Nbin,
            para, p0, stype, progressbar, fixed, symmetrize), ranges)

    zsum = nothing;
    zsumt = nothing;
    zsumd = nothing
    for res in results
        if zsum === nothing
            zsum = res.cs0;
            zsumt = res.cs1;
            zsumd = res.cs2
        else
            zsum = zsum .+ res.cs0
            if !isempty(res.cs1)
                zsumt .+= res.cs1
            end
            if !isnothing(res.cs2)
                zsumd .+= res.cs2
            end
        end
    end

    cs0 = zsum / num_workers
    cs1 = isempty(zsumt) ? [] : (zsumt / num_workers)
    cs2 = isnothing(zsumd) ? nothing : (zsumd / num_workers)

    return (cs0=cs0, cs1=cs1, cs2=cs2, axesV=results[1].axesV, laxes=results[1].laxes, proc=proc, stype=stype)
end

end  # module Xs