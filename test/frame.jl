include("../src/DalitzPlot.jl")
using .DalitzPlot, .DalitzPlot.qBSE, .DalitzPlot.FR, .DalitzPlot.GEN, .DalitzPlot.Xs
using StaticArrays

function LorentzBoost(k,p)
kp=k[1]*p[1]+k[2]*p[2]+k[3]*p[3]
knew1=  k[1]+p[1]*((kp/(p[5]+p[4])+k[4])/p[5])
knew2=  k[2]+p[2]*((kp/(p[5]+p[4])+k[4])/p[5])
knew3=  k[3]+p[3]*((kp/(p[5]+p[4])+k[4])/p[5])
knew4=(p[4]*k[4]+kp)/p[5]
knew5=k[5]
return @SVector [knew1,knew2,knew3,knew4,knew5  ]
end

function main()
mf=[1.,2.,1.5]
kf, wt = GENEV(5., mf)
k1,k2,k3=Xs.getkf(kf)
p=@MVector [k1[1],k1[2],k1[3],k2[4]+k3[4],sqrt((k2+k3)*(k2+k3))]
#U2=FR.U(k2, 1)
#U3=FR.U(k3, 1,bar=true)
#U2m=FR.U(k1, -1)
#U3m=FR.U(k1, -1,bar=true)
U2=FR.eps(k2, 0)
U3=FR.eps(k3, 0,star=true)
U21=FR.eps(k2, 1)
U31=FR.eps(k3, 1,star=true)
U2m=FR.eps(k2, -1)
U3m=FR.eps(k3, -1,star=true)
@show U3m*k3

k2new=LorentzBoost(k2,p)
k3new=LorentzBoost(k3,p)
k1new=LorentzBoost(k1,p)
U2new=FR.eps(k2new, 0)
U3new=FR.eps(k3new, 0,star=true)
U2new1=FR.eps(k2new, 1)
U3new1=FR.eps(k3new, 1,star=true)
U2newm=FR.eps(k2new, -1)
U3newm=FR.eps(k3new, -1,star=true)
@show U3newm*k3new

end
main()