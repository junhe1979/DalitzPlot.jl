# DalitzPlot

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/DalitzPlot.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/DalitzPlot.jl/dev/)
[![Coverage](https://codecov.io/gh/gridap/DalitzPlot.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gridap/DalitzPlot.jl)

This package is used to plot the Dalitz plot with a given amplitude.

## Requirements
Julia 1.6+ https://julialang.org/downloads/

## Installation

```julia
Pkg.add("DalitzPlot")
using DalitzPlot
```

## Examples

Provide the amplitudes. Here we take it as 1
```julia
amp(tecm, kf, ch, para)=1.
```
Provide the mass of initial and final particles.
```julia
ch = (mi=[1.0, 1.0], mf=[1.0, 1.0, 1.0],namei=["p^i_{1}", "p^i_{2}"], namef=["p^f_{1}", "p^f_{2}", "p^f_{3}"], amp=amp) 
```
Provide the momentum
```julia
p = 10.0
```
Calculate 
```julia
res = Xsection(plab2pcm(p, ch.mi), ch, nevtot=Int64(1e7), para=(p=p, l=1.0), ProgressBars=true)
```
Plot Dalitz Plot.
```julia
plotD(res, ch, axes=[1, 3])
```

![ex1.png](test/DP.png)



