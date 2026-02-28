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
function plotD(res; cg=cgrad([:white, :green, :blue, :red], [0, 0.01, 0.1, 0.5, 1.0]), topx=[], topy=[], topye=[], rightx=[], righty=[], rightye=[], toprightx=[], toprighty=[], toprightye=[], filename="DP.pdf")
    ENV["GKSwstype"] = "100"
    cs1 = res.cs1
    cs2 = res.cs2
    axesV = res.axesV
    laxes = res.laxes
    proc = res.proc

    laxes1 = [laxes0[1] for laxes0 in laxes]
    Laxes = [proc.namef[laxes1[1][1]] * proc.namef[laxes1[1][2]], proc.namef[laxes1[2][1]] * proc.namef[laxes1[2][2]]]


    Nbin = length(axesV[1])
    x1 = axesV[1]
    xlims = (minimum(x1), maximum(x1))
    dx = (maximum(x1) - minimum(x1)) / Nbin

    p0 = Plots.plot(xticks=:auto, ylabel=latexstring("d\\sigma/m^2_{" * Laxes[1] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, xmirror=true, legend=:none, linetype=:steppre)
    y1min, y1max = 1e20, 0.
    for i in eachindex(cs1)
        y1 = cs1[i][1, :]
        if minimum(y1) / dx < y1min
            y1min = minimum(y1) / dx
        end
        if maximum(y1) / dx > y1max
            y1max = maximum(y1) / dx
        end
        Plots.plot!(p0, x1, y1 / dx, xlims=xlims, ylims=(y1min, y1max * 1.1))
    end

    if !isempty(topx)
        p1 = if isempty(topye)
            Plots.scatter!(topx, topy, markersize=1)
        else
            Plots.scatter!(topx, topy, yerr=topye, markersize=1)
        end
    else
        p1 = p0
    end



    y2 = axesV[2]
    ylims = (minimum(y2), maximum(y2))
    dy = (maximum(y2) - minimum(y2)) / Nbin
    p0 = Plots.plot(xlabel=latexstring("d\\sigma/m^2_{" * Laxes[2] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, ymirror=true, legend=:none, linetype=:steppre)
    x2min, x2max = 1e20, 0.
    for i in eachindex(cs1)
        x2 = cs1[i][2, :]

        if minimum(x2) / dx < x2min
            x2min = minimum(x2) / dy
        end
        if maximum(x2) / dx > x2max
            x2max = maximum(x2) / dy
        end
        p0 = Plots.plot!(x2 / dy, y2, xlims=(x2min, x2max * 1.1), ylims=ylims)
    end


    if !isempty(topx)
        p2 = if isempty(rightye)
            Plots.scatter!(righty, rightx, markersize=1)
        else
            Plots.scatter!(righty, rightx, xerr=rightye, markersize=1)
        end
    else
        p2 = p0
    end



    x = axesV[1]
    y = axesV[2]
    xlims = (minimum(x), maximum(x))
    ylims = (minimum(y), maximum(y))
    dx = (maximum(x) - minimum(x)) / Nbin
    dy = (maximum(y) - minimum(y)) / Nbin
    z = [cs2[1][ix, iy] / (dx * dy) for iy in 1:Nbin, ix in 1:Nbin]
    p = Plots.heatmap(x, y, z, xlims=xlims, ylims=ylims, c=cg, xlabel=latexstring(Laxes[1]), ylabel=latexstring(Laxes[2]), framestyle=:box, cb=:none)

    if length(axesV) >= 3
        Nbin = length(axesV[3])
        x3 = axesV[3]
        xlims = (minimum(x3), maximum(x3))
        dx = (maximum(x3) - minimum(x3)) / Nbin

        p0 = Plots.plot(xticks=:auto, ylabel=latexstring("d\\sigma/m^2_{" * Laxes[1] * "} (\\textrm{ barn/GeV^2})"), framestyle=:box, xmirror=true, ymirror=true, legend=:none, linetype=:steppre)
        y3min, y3max = 1e20, 0.
        for i in eachindex(cs1)
            y3 = cs1[i][3, :]
            if minimum(y3) / dx < y3min
                y3min = minimum(y3) / dx
            end
            if maximum(y3) / dx > y3max
                y3max = maximum(y3) / dx
            end
            Plots.plot!(p0, x3, y3 / dx, xlims=xlims, ylims=(y3min, y3max * 1.1))
        end

        if !isempty(toprightx)
            p3 = Plots.scatter!(toprightx, toprighty, yerr=toprightye, markersize=1)
        else
            p3 = p0
        end

        l = @layout [a b
            b{0.7w,0.7h} c]
        DP = Plots.plot(p1, p3, p, p2, layout=l, titleloc=:left, titlefont=10, size=(880, 800), left_margin=1mm, right_margin=1mm, bottom_margin=1mm, top_margin=1mm, link=:all)
    else
        l = @layout [
            a _
            b{0.7w,0.7h} c
        ]
        DP = Plots.plot(p1, p, p2, layout=l, titleloc=:left, titlefont=10, size=(800, 800), left_margin=1mm, right_margin=1mm, bottom_margin=1mm, top_margin=1mm, link=:all)
    end

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
