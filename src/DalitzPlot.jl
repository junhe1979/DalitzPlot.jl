#*******************************************************************************************
# Code for Dalitz Plot and Qusipotential Bethe-Salpter equation  
# by Jun He @NNU   July 2024  
# @ https://github.com/junhe1979/DalitzPlot.jl
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

export GEN, Xs, FR, qBSE, PLOT, AUXs
end
