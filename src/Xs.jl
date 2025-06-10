module Xs
using StaticArrays, Distributed, ProgressBars
#include("GEN.jl")
using ..GEN
#############################################################################

function LorentzBoost(k::SVector{5,Float64}, p::SVector{5,Float64})
    kp = k[1] * p[1] + k[2] * p[2] + k[3] * p[3]  # 点积 k ⋅ p
    p5 = sqrt(p[4]^2 - p[1]^2 - p[2]^2 - p[3]^2)
    beta_factor = kp / (p[4] + p5) + k[4]       # β-related factor
    inv_p5 = 1.0 / p5


    k1 = k[1] + p[1] * beta_factor * inv_p5
    k2 = k[2] + p[2] * beta_factor * inv_p5
    k3 = k[3] + p[3] * beta_factor * inv_p5
    k4 = (p[4] * k[4] + kp) * inv_p5
    k5 = k[5]

    return @SVector [k1, k2, k3, k4, k5]
end

function LorentzBoost(momenta::Vector{SVector{5,Float64}}, p::SVector{5,Float64})
    return [LorentzBoost(i, p) for i in momenta]
end

function Rotation(k::SVector{5,Float64}, ct::Float64, st::Float64, cp::Float64, sp::Float64)
    k1 = ct * cp * k[1] + ct * sp * k[2] - st * k[3]
    k2 = -sp * k[1] + cp * k[2]
    k3 = st * cp * k[1] + st * sp * k[2] + ct * k[3]
    return @SVector [k1, k2, k3, k[4], k[5]]
end

function Rotation(momenta::Vector{SVector{5,Float64}}, ct::Float64, st::Float64, cp::Float64, sp::Float64)
    return [Rotation(k, ct, st, cp, sp) for k in momenta]
end

#############################################################################
function binx(i::Int64, bin, iaxis::Int64)::Float64
    return bin.min[iaxis][1] + (i - 0.5) / bin.Nbin * (bin.max[iaxis][1] - bin.min[iaxis][1])
end
function binrange(laxes::Vector{Int64}, tecm, ch, stype)
    min = (ch.mf[laxes[1]] + ch.mf[laxes[2]])^stype
    max = (tecm - sum(ch.mf) + ch.mf[laxes[1]] + ch.mf[laxes[2]])^stype
    return min, max
end

function binrange(laxes::Vector{Vector{Int64}}, tecm, ch, stype)
    min = Float64[(ch.mf[axes0[1]] + ch.mf[axes0[2]])^stype for axes0 in laxes]
    max = Float64[(tecm - sum(ch.mf) + ch.mf[axes0[1]] + ch.mf[axes0[2]])^stype for axes0 in laxes]
    return min, max
end

function binrange(laxes::Vector{Vector{Vector{Int64}}}, tecm, ch, stype; Range=[])
    minn = Vector{Float64}[]
    maxn = Vector{Float64}[]

    for i in eachindex(laxes)

        min0, max0 = binrange(laxes[i], tecm, ch, stype)

        # If the user provides min and max parameters, override the default calculated values
        if !isempty(Range)
            if !isempty(Range[i])
                if !ismissing(Range[i][1])
                    min0 .= Range[i][1]
                end

                if !ismissing(Range[i][2])
                    max0 .= Range[i][2]
                end
            end
        end
        push!(minn, min0)
        push!(maxn, max0)
    end


    return minn, maxn
end
function Nsij(kijs, min, max, Nbin)
    return convert(Int64, cld((kijs - min) * Nbin, (max - min)))
end

function Nsum3(laxes::Vector{Vector{Int64}}, i, bin::NamedTuple, kf, stype)
    Ns = Int64[]
    sij = Float64[]

    for iaxes in eachindex(laxes)
        axes0 = laxes[iaxes]

        kij = kf[axes0[1]] + kf[axes0[2]]
        kijs = stype == 2 ? kij * kij : sqrt(kij * kij)
        sij = push!(sij, kijs)

        push!(Ns, Nsij(kijs, bin.min[i][iaxes], bin.max[i][iaxes], bin.Nbin))
    end
    return Ns, sij
end

function Nsum3(laxes::Vector{Vector{Vector{Int64}}}, bin::NamedTuple, kf, stype)
    Ns = Vector{Vector{Int64}}(undef, length(laxes))
    sij = Vector{Vector{Float64}}(undef, length(laxes))
    for i in eachindex(laxes)
        Ns[i], sij[i] = Nsum3(laxes[i], i, bin, kf, stype)
    end
    return Ns, sij
end
#############################################################################
function plab2pcm(p::Float64, mi::Vector{Float64})
    ebm = sqrt(p^2 + mi[1]^2) #
    W = sqrt((ebm + mi[2])^2 - p^2)
    return W
end
function getkf(p::Float64, kf::Vector{SVector{5,Float64}}, ch::NamedTuple)

    ebm = sqrt(p^2 + ch.mi[1]^2)
    tecm = sqrt((ebm + ch.mi[2])^2 - p^2)

    GAM = (ebm + ch.mi[2]) / tecm
    ETA = p / tecm


    klab = SVector{5,Float64}[SVector{5,Float64}(
        GAM * k[1] + ETA * k[4],  # k1
        k[2],                    # k2
        k[3],                    # k3
        ETA * k[1] + GAM * k[4], # k4
        k[5]                     # k5
    ) for k in kf]

    return klab
end
function plab(p::Float64, mi::Vector{Float64})

    m1, m2 = mi[1], mi[2]
    p1 = SVector{5,Float64}([p 0.0 0.0 sqrt(p^2 + m1^2) m1]) #f0
    p2 = SVector{5,Float64}([0.0 0.0 0.0 m2 m2]) #p

    return p1, p2
end
function pcm(tecm::Float64, mi::Vector{Float64})
    m1, m2 = mi[1], mi[2]
    E1 = (tecm^2 + m1^2 - m2^2) / (2.0 * tecm)
    E2 = (tecm^2 + m2^2 - m1^2) / (2.0 * tecm)
    p = sqrt(E2^2 - m2^2)
    p1 = SVector{5,Float64}([p 0.0 0.0 E1 m1]) #f0
    p2 = SVector{5,Float64}([-p 0.0 0.0 E2 m2]) #p

    return p1, p2
end
#############################################################################
function Xsection(tecm, ch, callback; axes=[], Range=[], nevtot=Int64(1e6),
    Nbin=1000, para=(l = 1.0), p0=[], stype=1)
    Nf = length(ch.pf)
    laxes = Vector{Vector{Int64}}[]
    for axis in axes
        positions = [findall(==(element), ch.pf) for element in axis]
        laxes0 = Vector{Int64}[]
        for p1 in positions[1], p2 in positions[2]
            push!(laxes0, [p1, p2])
        end
        push!(laxes, laxes0)
    end


    Naxes = length(axes)
    axesV = []
    if Nf > 2
        min, max = binrange(laxes, tecm, ch, stype, Range=Range)
        bin = (Nbin=Nbin, min=min, max=max)

        axesV = [[binx(ix, bin, iaxes) for ix in 1:Nbin] for iaxes in eachindex(laxes)]
    end
    zsum = 0e0
    zsumt = zeros(Float64, Naxes, Nbin)
    zsumd = zeros(Float64, Nbin, Nbin)

    for ine in 1:nevtot

        kf, wt = GENEV(tecm, ch.mf)
        if Nf > 2
            Nsij, sij = Nsum3(laxes, bin, kf, stype)

            if all([any(min[i] .<= sij[i] .&& sij[i] .<= max[i]) for i in eachindex(min)])

                amp0 = ch.amps(tecm, kf, ch, para, p0)
                wt = wt * amp0
                if Nf > 2
                    for i in 1:Naxes
                        for isij in Nsij[i]
                            if 1 < isij <= Nbin
                                if i == 2
                                    zsumt[i, isij] += wt
                                else
                                    zsumt[i, isij] += wt
                                end

                            end
                        end
                    end
                    if Naxes >= 2
                        for isij in Nsij[1], jsij in Nsij[2]
                            if 1 < isij <= Nbin && 1 < jsij <= Nbin
                                zsumd[isij, jsij] += wt
                            end
                        end
                    end
                end
            end
        elseif Nf == 2
            amp0 = ch.amps(tecm, kf, ch, para, p0)
            wt = wt * amp0
        end

        zsum += wt
        callback(ine)
    end

    cs0 = zsum / nevtot
    cs1 = zsumt / nevtot
    cs2 = zsumd / nevtot
    res = (cs0=cs0, cs1=cs1, cs2=cs2, axesV=axesV, laxes=laxes, ch=ch)
    return res
end

function worker_Xsection(tecm, ch, axes, Range, nevt, Nbin, para, p0, stype, progressbar)
    # Define progress bar update function
    function progress_callback(pb)
        ProgressBars.update(pb)
    end

    if myid() == 2 && progressbar # Only the first process displays the progress bar
        pb = ProgressBar(1:nevt)  # Create a progress bar ranging from 1 to n
        callback = i -> progress_callback(pb)   # Callback function to update the progress bar
    else
        callback = _ -> nothing
    end
    return Xsection(tecm, ch, callback, axes=axes, Range=Range, nevtot=nevt, Nbin=Nbin, para=para, p0=p0, stype=stype)
end

function Xsection(tecm, ch; axes=[23, 21], Range=[], nevtot=Int64(1e6), Nbin=100, para=(l = 1.0), p0=[], stype=1, progressbar=true)
    num_workers = nworkers()
    nevt_per_worker = div(nevtot, num_workers)
    ranges = [(i * nevt_per_worker + 1, Base.min((i + 1) * nevt_per_worker, nevtot)) for i in 0:(num_workers-1)]
    GC.gc(false)
    results = pmap(r -> worker_Xsection(tecm, ch, axes, Range, r[2] - r[1] + 1, Nbin, para, p0, stype, progressbar), ranges)

    zsum = 0.0
    zsumt = nothing
    zsumd = nothing

    for res in results
        zsum += res.cs0
        if isnothing(zsumt)
            zsumt = res.cs1
        else
            zsumt .+= res.cs1
        end
        if !isnothing(res.cs2)
            if isnothing(zsumd)
                zsumd = res.cs2
            else
                zsumd .+= res.cs2
            end
        end
    end

    cs0 = zsum / num_workers
    cs1 = zsumt / num_workers
    cs2 = isnothing(zsumd) ? nothing : zsumd / num_workers

    return (cs0=cs0, cs1=cs1, cs2=cs2, axesV=results[1].axesV, laxes=results[1].laxes, ch=ch)
end


end
