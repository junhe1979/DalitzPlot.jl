
module DalitzPlot
export GENEV, Xsection, Xsection2, plotD, III, GA, GS, epsilon, Uc, Ubc, LCV, cdot, plab2pcm, getkf, plab, pcm
using StaticArrays, ProgressBars, Distributed, Plots, LaTeXStrings, Colors, Compose, DelimitedFiles

include("GEN.jl")
include("QFT.jl")
include("AuxiliaryFunction.jl")
using .GEN, .QFT, .AuxiliaryFunction


function Xsection(tecm, ch; nevtot=Int64(1e6), Nbin=100, para=(l = 1.0), ProgressBars=false)
    Nf = length(ch.mf)
    axes=[]
    
    if Nf > 2
        min, max = binrange(tecm, ch)
        bin = (Nbin=Nbin, min=min, max=max)
        axes = [[binx(ix, bin, iaxes) for ix in 1:Nbin] for iaxes in 1:3]
    end
    zsum = 0e0
    zsumt = zeros(Float64, 3, Nbin + 1)
    zsumd = zeros(Float64, 3, 3, Nbin + 1, Nbin + 1)
    if ProgressBars == true
        if nprocs() > 1
            ne = 1:nevtot
            if myid() == 2
                ne = ProgressBar(1:nevtot)
            end
        else
            ne = ProgressBar(1:nevtot)
        end
    else
        ne = 1:nevtot
    end

    for ine in ne

        kf, wt = GENEV(tecm, ch.mf)
        amp0 = ch.amp(tecm, kf, ch, para)
        wt = wt * amp0
        if Nf > 2
            Nsum = Nsum3(bin, kf)
            for i in 1:Nf
                zsumt[i, Nsum[i]] += wt
                for j in 1:Nf
                    zsumd[i, j, Nsum[i], Nsum[j]] += wt
                end
            end
        end
        zsum += wt
    end

    cs0 = zsum / nevtot
    cs1 = zsumt / nevtot
    cs2 = zsumd / nevtot
    res = (cs0=cs0, cs1=cs1, cs2=cs2, axes=axes)
    return res
end

function plotD(res, ch; axes=[1, 2], cg=cgrad([:white, :green, :blue, :red], [0, 0.01, 0.1, 0.5, 1.0]))

    Laxes = [ch.namef[2] * ch.namef[3], ch.namef[1] * ch.namef[3], ch.namef[1] * ch.namef[2]]
    cs0 = res[1]
    cs1 = res[2]
    cs2 = res[3]
    axesV = res[4]
    Nbin = length(axesV[1])
    x1 = [axesV[axes[1]][i] for i in 1:Nbin]
    y1 = [cs1[axes[1], i]  for i in 1:Nbin]
    xlims = (minimum(x1), maximum(x1))
    ylims = (minimum(y1), maximum(y1))
    dx=(maximum(x1)-minimum(x1))/Nbin 
    ylims = (minimum(y1)/dx, maximum(y1)/dx)
    @show dx
    p1 = Plots.plot(x1, y1/dx, xlims=xlims, ylims=ylims, xticks=:auto, ylabel=latexstring("d\\sigma/m^2_{" * Laxes[axes[1]] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, xmirror=true, legend=:none, linetype=:steppre)


    y2 = [axesV[axes[2]][i] for i in 1:Nbin]
    x2 = [cs1[axes[2], i]  for i in 1:Nbin]
    xlims = (minimum(x2), maximum(x2))
    ylims = (minimum(y2), maximum(y2))
    dy=(maximum(y2)-minimum(y2))/Nbin 
    xlims = (minimum(x2)/dy, maximum(x2)/dy)
    p2 = Plots.plot(x2/dy, y2, xlims=xlims, ylims=ylims, xlabel=latexstring("d\\sigma/m^2_{" * Laxes[axes[2]] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, ymirror=true, legend=:none, linetype=:steppre)

    x = [axesV[axes[1]][ix] for ix in 1:Nbin]
    y = [axesV[axes[2]][iy] for iy in 1:Nbin]
    xlims = (minimum(x), maximum(x))
    ylims = (minimum(y), maximum(y))
    dx=(maximum(x)-minimum(x))/Nbin 
    dy=(maximum(y)-minimum(y))/Nbin 
    z= [cs2[axes[1], axes[2], ix, iy]/(dx*dy)  for iy in 1:Nbin, ix in 1:Nbin]
    p3 = Plots.heatmap(x, y, z, xlims=xlims, ylims=ylims, c=cg, xlabel=latexstring(Laxes[axes[1]]), ylabel=latexstring(Laxes[axes[2]]), framestyle=:box, cb=:none)


    l = @layout [a _
        b{0.8w,0.8h} c]
    DP = Plots.plot(p1, p3, p2, layout=l, titleloc=:left, titlefont=10, size=(1000, 1000), left_margin=1mm, right_margin=1mm, bottom_margin=1mm, top_margin=1mm, link=:all)

    Plots.savefig("DP.png")
    return DP
end


end
