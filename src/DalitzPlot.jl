#*******************************************************************************************
# Code for  Qusipotential Bethe-Salpter equation  
# by Jun He @NNU   July 2024  
# @ https://github.com/junhe1979/DalitzPlot.jl
# Version 0.1.1: 26 Jul 2024 
# Version 0.2.0:  6 Oct 2024
# Version 0.2.7: 27 Nov 2024
# Version 0.3.0: 16 Feb 2025 current
#*******************************************************************************************

module DalitzPlot
include("GEN.jl")
include("Xs.jl")
include("FR.jl")
include("qBSE.jl")
include("PLOT.jl")
include("AUXs.jl")
using .GEN
using .Xs
using .FR
using .qBSE
using .PLOT
using .AUXs
# 顶层导出
export GEN, Xs, FR, qBSE, PLOT, AUXs
end
