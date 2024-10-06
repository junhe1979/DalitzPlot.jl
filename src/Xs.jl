
module Xs
using StaticArrays
include("GEN.jl")
using .GEN

function binx(i::Int64, bin, iaxis::Int64)::Float64
    return bin.min[iaxis] + (i - 0.5) / bin.Nbin * (bin.max[iaxis] - bin.min[iaxis])
end

function binrange(laxes, tecm, ch)
    min = Float64[(ch.mf[axes0[1]] + ch.mf[axes0[2]])^2  for axes0 in laxes]
    max = Float64[(tecm - sum(ch.mf) + ch.mf[axes0[1]] + ch.mf[axes0[2]])^2  for axes0 in laxes]
    return min, max
end

#bin
function Nsum3(laxes, bin::NamedTuple, kf::MMatrix{5,18,Float64,90})
    Nsum = zeros(Int64, 2)
    for iaxes in eachindex(laxes)
        axes0 = laxes[iaxes]
        kf1 = SVector{5,Float64}(kf[1:5, axes0[1]])
        kf2 = SVector{5,Float64}(kf[1:5, axes0[2]])

        k12 = kf1 + kf2
        k12s = k12*k12
        Nsum[iaxes] = convert(Int64, cld((k12s - bin.min[iaxes]) * bin.Nbin, (bin.max[iaxes] - bin.min[iaxes])))
    end
    return Nsum
end
#############################################################################
#Transfer
#############################################################################

function plab2pcm(p::Float64, mi::Vector{Float64})
    ebm = sqrt(p^2 + mi[1]^2) #
    W = sqrt((ebm + mi[2])^2 - p^2)
    return W
end

function getkf(p, kf::MMatrix{5,18,Float64,90}, ch::NamedTuple)

    ebm = sqrt(p^2 + ch.mi[1]^2) #
    tecm = sqrt((ebm + ch.mi[2])^2 - p^2)

    GAM = (ebm + ch.mi[2]) / tecm
    ETA = p / tecm
    for i in 1:3
        EI = kf[4, i]
        PX = kf[1, i]
        kf[1, i] = GAM * PX + ETA * EI
        kf[4, i] = ETA * PX + GAM * EI
        #kf[5, i] = ch.mf[i] #注意GENEV产生的kf[5]是三动量的大小
    end

    k1 = SVector{5,Float64}(kf[:, 1])
    k2 = SVector{5,Float64}(kf[:, 2])
    k3 = SVector{5,Float64}(kf[:, 3])

    return k1, k2, k3
end

function getkf(kf::MMatrix{5,18,Float64,90})

    k1 = SVector{5,Float64}(kf[:, 1])
    k2 = SVector{5,Float64}(kf[:, 2])
    k3 = SVector{5,Float64}(kf[:, 3])

    return k1, k2, k3
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

function Xsection(tecm, ch, callback; axes=[23, 21], nevtot=Int64(1e6), Nbin=100, para=(l = 1.0))
    Nf = length(ch.mf)
    axesV = []

    if Nf > 2
        laxes = [[div(axes[1], 10), mod(axes[1], 10)], [div(axes[2], 10), mod(axes[2], 10)]]

        min, max = binrange(laxes, tecm, ch)
        bin = (Nbin=Nbin, min=min, max=max)
        axesV = [[binx(ix, bin, iaxes) for ix in 1:Nbin] for iaxes in eachindex(laxes)]
    end
    zsum = 0e0
    zsumt = zeros(Float64, 2, Nbin)
    zsumd = zeros(Float64, Nbin , Nbin)

    for ine in 1:nevtot

        kf, wt = GENEV(tecm, ch.mf)
        amp0 = ch.amp(tecm, kf, ch, para)
        wt = wt * amp0
        if Nf > 2
            Nsum = Nsum3(laxes, bin, kf)
            for i in 1:2
                zsumt[i, Nsum[i]] += wt
            end
            zsumd[Nsum[1], Nsum[2]] += wt
        end
        zsum += wt
        callback(ine)
    end

    cs0 = zsum / nevtot
    cs1 = zsumt / nevtot
    cs2 = zsumd / nevtot
    res = (cs0=cs0, cs1=cs1, cs2=cs2, axesV=axesV, axes=axes,ch=ch)
    return res
end

end
