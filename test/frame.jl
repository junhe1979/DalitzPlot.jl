include("../src/DalitzPlot.jl")
using .DalitzPlot, .DalitzPlot.qBSE, .DalitzPlot.FR, .DalitzPlot.GEN, .DalitzPlot.Xs
using StaticArrays

function LorentzBoost(k, p)
    kp = k[1] * p[1] + k[2] * p[2] + k[3] * p[3]
    k1 = k[1] + p[1] * ((kp / (p[5] + p[4]) + k[4]) / p[5])
    k2 = k[2] + p[2] * ((kp / (p[5] + p[4]) + k[4]) / p[5])
    k3 = k[3] + p[3] * ((kp / (p[5] + p[4]) + k[4]) / p[5])
    k4 = (p[4] * k[4] + kp) / p[5]
    k5 = k[5]
    return @SVector [k1, k2, k3, k4, k5]
end

function LorentzRotation(k, ct, st, cp, sp)
    k1 = ct*cp*k[1] + ct*sp*k[2]-st*k[3]
    k2 = -sp*k[1]+cp*k[2]
    k3 = st*cp*k[1] + st*sp*k[2]+ct*k[3]
    return @SVector [k1, k2, k3, k[4], k[5]]
end


function main()
    mf = [1.0, 2.0, 1.5]
    kf, wt = GENEV(5.0, mf)
    k1, k2, k3 = Xs.getkf(kf)
    @show k2
    @show k3


    p = @MVector [k1[1], k1[2], k1[3], k2[4] + k3[4], sqrt((k2 + k3) * (k2 + k3))]
    #U2=FR.U(k2, 1)
    #U3=FR.U(k3, 1,bar=true)
    #U2m=FR.U(k1, -1)
    #U3m=FR.U(k1, -1,bar=true)

    amp2 = 0.0
    for l2 in -1:1, l3 in -1:1
        U2 = FR.eps(k2, l2)
        U3 = FR.eps(k3, l3, star=true)
        amp2 += (U3 * U2)
    end
    @show amp2
    @show k2 * k3

    amp2 = 0.0
    k2new = LorentzBoost(k2, p)
    k3new = LorentzBoost(k3, p)
    ct, st, cp, sp=FR.kph(k2new)
    k2new= LorentzRotation(k2new,ct, st, cp, sp)
    k3new= LorentzRotation(k3new,ct, st, cp, sp)
    @show k2new
    @show k3new
    for l2 in -1:1, l3 in -1:1
        U2 = FR.eps(k2new, l2)
        U3 = FR.eps(k3new, l3, star=true)
        amp2 += (U3 * U2)
    end
    @show amp2
    @show k2new * k3new

end
main()