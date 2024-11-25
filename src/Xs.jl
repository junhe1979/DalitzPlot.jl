module Xs
using StaticArrays, Distributed, ProgressBars
#include("GEN.jl")
using ..GEN
#############################################################################

function LorentzBoost(k::SVector{5, Float64}, p::SVector{5, Float64})
    kp = k[1] * p[1] + k[2] * p[2] + k[3] * p[3]  # 点积 k ⋅ p
    beta_factor = kp / (p[4] + p[5]) + k[4]       # β-related factor

    # 计算新的动量和能量分量
    k1 = k[1] + p[1] * beta_factor / p[5]
    k2 = k[2] + p[2] * beta_factor / p[5]
    k3 = k[3] + p[3] * beta_factor / p[5]
    k4 = (p[4] * k[4] + kp) / p[5]
    k5 = k[5]  # 不变的标量

    # 返回新的 SVector
    return @SVector [k1, k2, k3, k4, k5]
end

function LorentzBoost(momenta::Vector{SVector{5, Float64}}, p::SVector{5, Float64})
    return [LorentzBoost(i, p) for i in momenta]
end

function Rotation(k::SVector{5, Float64}, ct::Float64, st::Float64, cp::Float64, sp::Float64)
    k1 = ct * cp * k[1] + ct * sp * k[2] - st * k[3]
    k2 = -sp * k[1] + cp * k[2]
    k3 = st * cp * k[1] + st * sp * k[2] + ct * k[3]
    return @SVector [k1, k2, k3, k[4], k[5]]
end

function Rotation(momenta::Vector{SVector{5, Float64}}, ct::Float64, st::Float64, cp::Float64, sp::Float64)
    # 对每个动量应用旋转，返回新数组
    return [Rotation(k, ct, st, cp, sp) for k in momenta]
end

#############################################################################
function binx(i::Int64, bin, iaxis::Int64)::Float64
    return bin.min[iaxis] + (i - 0.5) / bin.Nbin * (bin.max[iaxis] - bin.min[iaxis])
end
function binrange(laxes, tecm, ch, stype)
    min = Float64[(ch.mf[axes0[1]] + ch.mf[axes0[2]])^stype for axes0 in laxes]
    max = Float64[(tecm - sum(ch.mf) + ch.mf[axes0[1]] + ch.mf[axes0[2]])^stype for axes0 in laxes]
    return min, max
end
function Nsij(kijs, min, max, Nbin)
    return convert(Int64, cld((kijs - min) * Nbin, (max - min)))
end

function Nsum3(laxes, bin::NamedTuple, kf, stype)
    Ns = zeros(Int64, 2)
    sij = zeros(Float64, 2)
    for iaxes in eachindex(laxes)
        axes0 = laxes[iaxes]

        kij = kf[axes0[1]] + kf[axes0[2]]
        kijs = stype == 2 ? kij * kij : sqrt(kij * kij)
        sij[iaxes] = kijs
        #Nsij[iaxes] = convert(Int64, cld((kijs - bin.min[iaxes]) * bin.Nbin, (bin.max[iaxes] - bin.min[iaxes])))
        Ns[iaxes] = Nsij(kijs, bin.min[iaxes], bin.max[iaxes], bin.Nbin)
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
    # 预计算常量
    ebm = sqrt(p^2 + ch.mi[1]^2)
    tecm = sqrt((ebm + ch.mi[2])^2 - p^2)

    GAM = (ebm + ch.mi[2]) / tecm
    ETA = p / tecm

    # 使用静态数组和广播计算
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
function Xsection(tecm, ch, callback; axes=[23, 21], min=[], max=[], nevtot=Int64(1e6),
    Nbin=100, para=(l = 1.0), stype=1)
    Nf = length(ch.mf)
    Naxes = length(axes)
    Nmax = length(max)
    Nmin = length(min)
    axesV = []

    if Nf > 2
        laxes = [[div(axes0, 10), mod(axes0, 10)] for axes0 in axes]
        mins, maxs = binrange(laxes, tecm, ch, stype)

        if Nmin == 0
            min = mins
        end
        if Nmax == 0
            max = maxs
        end

        bin = (Nbin=Nbin, min=min, max=max)
        axesV = [[binx(ix, bin, iaxes) for ix in 1:Nbin] for iaxes in eachindex(laxes)]
    end
    zsum = 0e0
    zsumt = zeros(Float64, Naxes, Nbin)
    if Naxes >= 2
        zsumd = zeros(Float64, Nbin, Nbin)
    end

    for ine in 1:nevtot

        kf, wt = GENEV(tecm, ch.mf)
        #@show wt
        if Nf > 2
            Nsij, sij = Nsum3(laxes, bin, kf, stype)
        end
        if all(min .<= sij) && all(sij .<= max)
            amp0 = ch.amp(tecm, kf, ch, para)
            wt = wt * amp0
            if Nf > 2
                for i in 1:Naxes
                    zsumt[i, Nsij[i]] += wt
                end
                if Naxes >= 2
                    zsumd[Nsij[1], Nsij[2]] += wt
                end
            end
        end

        zsum += wt
        callback(ine)
    end

    cs0 = zsum / nevtot
    cs1 = zsumt / nevtot
    cs2 = zsumd / nevtot
    res = (cs0=cs0, cs1=cs1, cs2=cs2, axesV=axesV, axes=axes, ch=ch)
    return res
end

function worker_Xsection(tecm, ch, axes, min_vals, max, nevt, Nbin, para, stype)
    # 定义进度条更新函数
    function progress_callback(pb)
        ProgressBars.update(pb)  # 更新进度条
    end
    # 检查是否是主进程
    if myid() == 2  # 仅第一个进程显示进度条
        pb = ProgressBar(1:nevt)  # 创建进度条，范围从1到n
        callback = i -> progress_callback(pb)   # 回调函数更新进度条
    else
        callback = _ -> nothing  # 其他进程不显示进度条
    end
    return Xsection(tecm, ch, callback, axes=axes, min=min_vals, max=max, nevtot=nevt, Nbin=Nbin, para=para, stype=stype)
end



function Xsection(tecm, ch; axes=[23, 21], min_vals=[], max=[], nevtot=Int64(1e6), Nbin=100, para=(l = 1.0), stype=1)
    num_workers = nworkers()
    nevt_per_worker = div(nevtot, num_workers)
    ranges = [(i * nevt_per_worker + 1, Base.min((i + 1) * nevt_per_worker, nevtot)) for i in 0:(num_workers-1)]

    # 使用 pmap 执行任务并获取结果
    results = pmap(r -> worker_Xsection(tecm, ch, axes, min_vals, max, r[2] - r[1] + 1, Nbin, para, stype), ranges)

    # 汇总结果
    zsum = 0.0
    zsumt = nothing
    zsumd = nothing

    for res in results
        zsum += res.cs0 * nevtot
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

    cs0 = zsum / nevtot
    cs1 = zsumt / nevtot
    cs2 = isnothing(zsumd) ? nothing : zsumd / nevtot

    return (cs0=cs0, cs1=cs1, cs2=cs2, axesV=results[1].axesV, axes=results[1].axes, ch=ch)
end


end
