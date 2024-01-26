
module AuxiliaryFunction 

export binx, binrange,Nsum2,Nsum3,plab2pcm,getkf,plab,pcm
using StaticArrays


function cdot(Q::SVector{5,Float64}, W::SVector{5,Float64})::Float64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

function binx(i::Int64, bin, iaxis::Int64)::Float64
    return bin.min[iaxis] + (i - 0.5) / bin.Nbin * (bin.max[iaxis] - bin.min[iaxis])
end

function binrange(tecm, ch)
    min = Float64[(ch.mf[2]+ch.mf[3])^2*0.99,
        (ch.mf[1]+ch.mf[3])^2*0.99,
        (ch.mf[1]+ch.mf[2])^2*0.99]
    max = Float64[(tecm-ch.mf[1])^2*1.01,
        (tecm-ch.mf[2])^2*1.01,
        (tecm-ch.mf[3])^2*1.01]
    return min, max
end
#bin
function Nsum3(bin::NamedTuple, kf::MMatrix{5,18,Float64,90})
    Nsum = zeros(Int64, 3)
    kf1=SVector{5,Float64}(kf[1:5, 1])
    kf2=SVector{5,Float64}(kf[1:5, 2])
    kf3=SVector{5,Float64}(kf[1:5, 3])

    k23 = kf2 + kf3
    k1 = cdot(k23, k23)
    k13 = kf1 + kf3
    k2 = cdot(k13, k13)
    k12 = kf1 + kf2
    k3 = cdot(k12, k12)

    Nsum[1] = convert(Int64, cld((k1 - bin.min[1]) * bin.Nbin, (bin.max[1] - bin.min[1])))
    Nsum[2] = convert(Int64, cld((k2 - bin.min[2]) * bin.Nbin, (bin.max[2] - bin.min[2])))
    Nsum[3] = convert(Int64, cld((k3 - bin.min[3]) * bin.Nbin, (bin.max[3] - bin.min[3])))
    return Nsum
end

#############################################################################
#Transfer
#############################################################################

function plab2pcm(p::Float64,mi::Vector{Float64})
    ebm =  sqrt(p^2 + mi[1]^2) #
    W = sqrt((ebm + mi[2])^2 - p^2)
    return W
end

function getkf(p, kf::MMatrix{5,18,Float64,90},ch::NamedTuple)

    ebm =  sqrt(p^2 + ch.mi[1]^2) #
    tecm = sqrt((ebm + ch.mi[2])^2 - p^2)

    GAM = (ebm + ch.mi[2]) / tecm
    ETA = p / tecm
    for i in 1:3
        EI = kf[4, i]
        PX = kf[1, i]
        kf[1, i] = GAM * PX + ETA * EI
        kf[4, i] = ETA * PX + GAM * EI
        kf[5, i] = ch.mf[i] #注意GENEV产生的kf[5]是三动量的大小
    end

    k1 = SVector{5,Float64}(kf[:, 1]) 
    k2 = SVector{5,Float64}(kf[:, 2]) 
    k3 = SVector{5,Float64}(kf[:, 3]) 

    return k1,k2,k3
end

function getkf(kf::MMatrix{5,18,Float64,90})

    k1 = SVector{5,Float64}(kf[:, 1]) 
    k2 = SVector{5,Float64}(kf[:, 2]) 
    k3 = SVector{5,Float64}(kf[:, 3]) 

    return k1,k2,k3
end

function plab(p::Float64, mi::Vector{Float64})

    m1,m2=mi[1],mi[2]
    p1 = SVector{5,Float64}([p 0.0 0.0 sqrt(p^2 + m1^2) m1]) #f0
    p2 = SVector{5,Float64}([0.0 0.0 0.0 m2 m2]) #p
  
    return p1,p2   
end

function pcm(tecm::Float64, mi::Vector{Float64})
    m1,m2=mi[1],mi[2]
    E1=(tecm^2+m1^2-m2^2)/(2. *tecm)
    E2=(tecm^2+m2^2-m1^2)/(2. *tecm)
    p=sqrt(E2^2-m2^2)
    p1 = SVector{5,Float64}([ p 0.0 0.0 E1 m1]) #f0
    p2 = SVector{5,Float64}([-p 0.0 0.0 E2 m2]) #p
  
    return p1,p2   
end

end
