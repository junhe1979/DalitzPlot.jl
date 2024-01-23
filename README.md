# DalitzPlot

This Julia package, designed for plotting Dalitz plots, allows users to visualize particle three-body decays with a specified amplitude. The package offers flexibility in amplitude settings, enabling users to customize the plot according to their specific requirements. 

## Installation

To install the package, use the standard Julia package manager procedure:

```julia
Pkg.add("DalitzPlot")
using DalitzPlot
```

## Examples

Provide the amplitudes. Here we take it as 1
```julia
amp(tecm, kf, ch, para)=1.
```

Define more complicated amplitudes for a 2->3 process. The "tecm" is the total energy in the center-of-mass frame. The "kf" is the momenta generated. The "ch" is the information about the process. The "para" contains extra parameters.

```julia
function amp(tecm, kf, ch, para)
    # get kf as momenta in the center-of-mass ,
    #k1,k2,k3=getkf(kf)       
    #get kf as momenta in laboratory frame
    k1, k2, k3 = getkf(para.p, kf, ch)

    # Incoming particle momentum
    # Center-of-mass frame: p1 = [p 0.0 0.0 E1]
    #p1, p2 = pcm(tecm, ch.mi)
    # Laboratory frame
    p1, p2 = plab(para.p, ch.mi)

    #flux
    #flux factor for cross section
    fac = 1 / (4 * para.p * ch.mi[2] * (2 * pi)^5)

    k12 = k1 + k2
    s12 = cdot(k12, k12)
    m = 3.0
    A = 1 / (s12 - m^2 + im * m * 0.1)

    total = abs2(A) * fac* 0.389379e-3

    return total

```
Provide the mass of initial and final particles. The name of the particles can be also provided for PlotD.
```julia
ch = (mi=[1.0, 1.0], mf=[1.0, 1.0, 1.0],namei=["p^i_{1}", "p^i_{2}"], namef=["p^f_{1}", "p^f_{2}", "p^f_{3}"], amp=amp) 
```
Provide the momentum. Here the "p" is the momentum in Laboratory frame.
```julia
p = 10.0
```
Calculate:  The function `plab2pcm` transforms the momentum of the incoming particle from the Laboratory frame to the total energy in the center-of-mass frame, a crucial step in the subsequent calculations. The variable "nevtot" represents the total number of events produced.

```julia
res = Xsection(plab2pcm(p, ch.mi), ch, nevtot=Int64(1e7), para=(p=p, l=1.0), ProgressBars=true)
```
The results of the calculations are stored in the variable `res` as a NamedTuple. Specifically, `res.cs0` corresponds to the total cross section, `res.cs1` represents the invariant mass spectrum, and `res.cs2` captures the data for the Dalitz plot.

Plot Dalitz Plot.
```julia
plotD(res, ch, axes=[1, 3])
```

![ex1.png](test/DP.png)
