#*******************************************************************************************
# Code for  Qusipotential Bethe-Salpter equation
# by Jun He @NNU   July 2024  Version: 29 July 2024
#*******************************************************************************************
#paprameters and defintions
#*******************************************************************************************
using LinearAlgebra, Printf,Plots, PythonCall
using DalitzPlot
using .DalitzPlot.qBSE, .DalitzPlot.QFT
qBSE.particles!(qBSE.p, "data/particlesjl.txt")
#*******************************************************************************************
function Output(Nih, Vc, Gc, II, Ec, Pole, IH)

    VGI = II - Vc * Gc
    # Compute determinant of (1 - VGc)
    detVGI = det(VGI)
    # Inverse of (1 - VGc)
    VGI0_inv = inv(VGI)
    # Define T matrix
    T = VGI0_inv * Vc
    # Write the results
    if abs(imag(Ec)) < 0.000001
        @printf("Output: Ec=%.4f, logD=%.4f\n", real(Ec), log(abs(detVGI)^2))
    end
    open("outputjl.txt", "a") do file
        write(file, @sprintf("%.4f, %.4f, %.4f\n", real(Ec), imag(Ec) * 1e3, log(abs(detVGI)^2)))
    end
    if Pole.Ampmin > log(abs(detVGI)^2)
        Pole.Ampmin = log(abs(detVGI)^2)
        Pole.Ampminx = real(Ec)
        Pole.Ampminy = imag(Ec) * 1e3
    end
    return Pole
end
function fV(k, l, ilV, LLi, LLf, lregu, lFFex, qt) #慢

    U2i = QFT.Uc(k.i2, ilV * l.i2)  # Define srUc_h based on the original function
    U2f = QFT.Ubc(k.f2, l.f2)  # Define srUc_h based on the original function
    GSk1 = QFT.GS(k.i1 + k.f1)

    FF::ComplexF64 = qBSE.fPropFF(k, 1, (LLi + LLf) / 2.0, LLi, LLf, lregu, lFFex, qt)  # Define fPropFF based on the original function

    fV = -1.0 / 4.0 / (0.85 * 0.093)^2 * FF * (U2f * GSk1 * U2i)
    #@show FF, chk.qd2

    return fV
end
#*******************************************************************************************
function main()
    #----------------------------------------------------------------
    dashline = repeat('-', 100)
    println(dashline)
    println("Start project: ")
    #________________________________________________________________   
    Np = 10 #for gausslaguerre
    lRm = 1 #Riemann sheet
    lregu = 1 #regulation
    lFFex = 0 #type of form factors
    qn = qBSE.QuantumNumberType(1, 2, 1, 2, -1, -1) #I,I_h,J,J_h,P,C               
    Range = qBSE.RangeType(1.2, 1.6, 50, 0.201, 10) #ERmin,ERmax,NER, EIt,   NEI
    #----------------------------------------------------------------
    println(dashline)
    println(dashline)
    println("I(J,P)=$(qn.lII)($(qn.lJJ),$(qn.lPP))  regu: exp      FFex=$lFFex")
    #________________________________________________________________   
    channels = (
        ("K", "N", qBSE.p[qBSE.pkey["K"]].m + 0.22 * 1.63 * 0.5),
        ("pi", "Sigma", qBSE.p[qBSE.pkey["pi"]].m + 0.22 * 1.63 * 0.5))
    CC = Dict(
        "K:N-->K:N" => (1, [qBSE.pkey["V"] 1], 3.0),
        "K:N-->pi:Sigma" => (1, [qBSE.pkey["V"] 1], -sqrt(1.5)),
        "pi:Sigma-->K:N" => (1, [qBSE.pkey["V"] 1], -sqrt(1.5)),
        "pi:Sigma-->pi:Sigma" => (1, [qBSE.pkey["V"] 1], 4.0))
    Nih, Nc, CB, IH, IA = qBSE.Independent_amp(channels, CC, qn.lJJ, Np)
    #----------------------------------------------------------------
    for iNc in 1:Nih
        println("")
        println("$(CB[IH[iNc].ich].name[3]) cutoff = $(CB[IH[iNc].ich].cutoff) GeV: $(IH[iNc].hel[1]) / $(IH[iNc].hel_lh[1])  $(IH[iNc].hel[2]) / $(IH[iNc].hel_lh[2])")
    end

    println(dashline)
    println("channel        ", join([CB[i].name[3] for i in 1:Nc], " "))
    for i1 in 1:Nc
        println(CB[i1].name[3], join([IA[i1, i].Nex for i in 1:Nc], " "))
    end
    println(dashline)
    println("ER=$(Range.ERmax) to $(Range.ERmin) GeV, NER=$(Range.NER)  ;   EI=$(Range.EIt*1e3) MeV, NEI=$(Range.NEI)")
    println(dashline)
    labelname = [CB[itemp].name[1] for itemp in 1:Nc]
    println("E  D= :: ", join(labelname[1:Nc], " "))
    #__________________calculate the bound energy____________________
    # We will use two loops to scan all real and image part of the energy,
    if isfile("outputjl.txt")
        run(`rm outputjl.txt`)
    end
    Pole = qBSE.PoleType(1e5, 0.0, 0.0)
    for iER in 0:Range.NER
        ER = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / Range.NER
        for iEI in -Range.NEI:Range.NEI-1
            EI = iEI * Range.EIt / Range.NEI
            Ec = ER + EI * Complex{Float64}(0, 1)
            # Calculate the V, G, and unit matrix II
            Vc, Gc, II = qBSE.srAB(Ec, qn, lregu, lFFex, lRm, Np, Nih, CB, IH, IA, fV)
            # Print the results.
            Pole = Output(Nih, Vc, Gc, II, Ec, Pole, IH)
        end
    end
    #----------------------------------------------------------------
    println("I(J,P)=$(qn.lII)($(qn.lJJ),$(qn.lPP)) reg=  ex= $(lFFex) pole= $(Pole.Ampmin) at $((Pole.Ampminx) * 1e3), $(Pole.Ampminy)")
    println(dashline)
    println("END A LOOP")
    println(dashline)

end
@time main()
include("plot.jl")
println("END PROGRAM")