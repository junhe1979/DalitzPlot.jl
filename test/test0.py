#!/usr/bin/env python
# coding=utf-8
from juliacall import Main as jl
jl.seval("using StaticArrays")
jl.seval("using DalitzPlot")
jl.seval("using DalitzPlot.FR")
jl.seval("using DalitzPlot.GEN")
jl.seval("using DalitzPlot.Xs")
jl.seval("using DalitzPlot.qBSE")
jl.seval("using DalitzPlot.plot")

mf=jl.Vector[jl.Float64]([1.0, 2.0, 1.5])
kf, wt = jl.GENEV(5.0, mf)
k1, k2, k3 = jl.Xs.getkf(kf)

p = jl.MVector(k1[1], k1[2], k1[3], k2[4] + k3[4], jl.sqrt((k2 + k3) * (k2 + k3)))

print(k1,p)