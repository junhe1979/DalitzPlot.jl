module PLOT
using Plots, LaTeXStrings, Colors, Compose, DelimitedFiles
function readdata(filename)
    data_groups = []              # Stores all data groups
    current_group = []           # Temporary storage for current data group

    open(filename, "r") do file
        for line in eachline(file)
            # Skip lines starting with #
            if startswith(strip(line), "#")
                continue
            end

            if isempty(strip(line))  # Check for blank line
                # If current group is not empty, we've reached the end of a group
                if !isempty(current_group)
                    push!(data_groups, hcat(current_group...)')  # Save current group as matrix
                    empty!(current_group)  # Clear temporary array for the next group
                end
            else
                # Parse the line into Float64 numbers and add to current group
                row = [parse(Float64, x) for x in split(line)]
                push!(current_group, row)
            end
        end
        # At the end, check if there is an unfinished group to be saved
        if !isempty(current_group)
            push!(data_groups, hcat(current_group...)')
        end
    end
    return data_groups
end



function plotD(res; cg=cgrad([:white, :green, :blue, :red], [0, 0.01, 0.1, 0.5, 1.0]), xx=[], xy=[], yx=[], yy=[], filename="DP.pdf")
    ENV["GKSwstype"] = "100"
    cs1 = res.cs1
    cs2 = res.cs2
    axesV = res.axesV
    laxes = res.laxes
    ch = res.ch

    laxes1 = [laxes0[1] for laxes0 in laxes]
    Laxes = [ch.namef[laxes1[1][1]] * ch.namef[laxes1[1][2]], ch.namef[laxes1[2][1]] * ch.namef[laxes1[2][2]]]


    Nbin = length(axesV[1])
    x1 = axesV[1]
    y1 = cs1[1, :]
    xlims = (minimum(x1), maximum(x1))
    ylims = (minimum(y1), maximum(y1))
    dx = (maximum(x1) - minimum(x1)) / Nbin
    ylims = (minimum(y1) / dx, maximum(y1) / dx * 1.1)

    p0 = Plots.plot(x1, y1 / dx, xlims=xlims, ylims=ylims, xticks=:auto, ylabel=latexstring("d\\sigma/m^2_{" * Laxes[1] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, xmirror=true, legend=:none, linetype=:steppre)
    if !isempty(xx)
        p1 = Plots.scatter!(xx, xy, markersize=1)
    else
        p1 = p0
    end

    y2 = axesV[2]
    x2 = cs1[2, :]

    xlims = (minimum(x2), maximum(x2))
    ylims = (minimum(y2), maximum(y2))
    dy = (maximum(y2) - minimum(y2)) / Nbin
    xlims = (minimum(x2) / dy, maximum(x2) / dy * 1.1)
    p0 = Plots.plot(x2 / dy, y2, xlims=xlims, ylims=ylims, xlabel=latexstring("d\\sigma/m^2_{" * Laxes[2] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, ymirror=true, legend=:none, linetype=:steppre)
    if !isempty(yy)
        p2 = Plots.scatter!(yy, yx, markersize=1)
    else
        p2 = p0
    end

    x = axesV[1]
    y = axesV[2]
    xlims = (minimum(x), maximum(x))
    ylims = (minimum(y), maximum(y))
    dx = (maximum(x) - minimum(x)) / Nbin
    dy = (maximum(y) - minimum(y)) / Nbin
    z = [cs2[ix, iy] / (dx * dy) for iy in 1:Nbin, ix in 1:Nbin]
    p3 = Plots.heatmap(x, y, z, xlims=xlims, ylims=ylims, c=cg, xlabel=latexstring(Laxes[1]), ylabel=latexstring(Laxes[2]), framestyle=:box, cb=:none)


    l = @layout [a _
        b{0.8w,0.8h} c]
    DP = Plots.plot(p1, p3, p2, layout=l, titleloc=:left, titlefont=10, size=(1000, 1000), left_margin=1mm, right_margin=1mm, bottom_margin=1mm, top_margin=1mm, link=:all)

    Plots.savefig(filename)
    return DP
end
function plotWeb(res)
    cs0, cs1, cs2 = res.cs0, res.cs1, res.cs2
    Nbin = length(res.axesV[1])
    x, y = res.axesV
    xlims, ylims = (minimum(x), maximum(x)), (minimum(y), maximum(y))
    dx, dy = (maximum(x) - minimum(x)) / Nbin, (maximum(y) - minimum(y)) / Nbin
    y1, y2 = cs1[1, :] / dx, cs1[2, :] / dy
    z = [[cs2[ix, iy] / (dx * dy) for iy in 1:Nbin] for ix in 1:Nbin]
    return x, y, z, y1, y2
end
end
