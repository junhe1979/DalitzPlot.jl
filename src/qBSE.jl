# module for quasipotential Bethe-Salpeter equation
module qBSE
using Distributed, FastGaussQuadrature, StaticArrays, WignerD, LinearAlgebra, Printf
using ..Xs, ..FR
#*******************************************************************************************
#store informations about system 
struct structSys #SYS
    Sys::String # label of system
    kv::Vector{Float64} #momenta of discreting
    wv::Vector{Float64} #weights of discreting
    xv::Vector{Float64} #cos theta of discreting
    wxv::Vector{Float64} #weights of discreting
    d::Vector{Matrix{Float64}}   #Wignerd
    pv::Vector{Float64} #phi of discreting
    wpv::Vector{Float64} #weights of discreting
    sp::Vector{Float64} #sin phi
    cp::Vector{Float64} #cos phi
    cutoff_re_type::Symbol
    cutoff_ex::Float64
    cutoff_ex_type::Symbol
    FF_ex_type::Int64
end
# store the information of a interaction
struct structInterAction #IA[]
    Nex::Int64  #total number of exchanges
    name_ex::Vector{String}  #exchange
    key_ex::Vector{Int64} #label
    J_ex::Vector{Int64}  #J for integral spin, 2J for half  integral spin
    Jh_ex::Vector{Int64} #1 for integral spin, 2 for half  integral spin
    P_ex::Vector{Int64} #parity
    m_ex::Vector{Float64} #mass
    dc::Vector{Int64}  #direct or cross
    Ff::Vector{Float64}  #flavor factors
end
# store the information of a channel.
struct structChannel #CH[] 
    p::Vector{Int64} #particles of channel
    m::Vector{Float64} #particle masses of channel
    J::Vector{Int64} #particle spins of channel
    Jh::Vector{Int64} #particle  spins of channel
    P::Vector{Int64} #particle parities of channel
    cutoff::Float64 #cutoff of channel
    name::Vector{String}  #name of channel
    name0::Vector{String} #name of channel (without charge)
    IHb::Int64 # rank of first independent helicity
    IHe::Int64 # rank of last independent helicity
    IHn::Int64 #total number of helicity amplitudes of a channel
end
# store the information of an independent helicity for matrix
mutable struct structIndependentHelicity #IH[]
    iCH::Int64 # which channel this independent helicity belongs to 
    hel::Vector{Int64} # helicities
    helh::Vector{Int64} # fermion or boson
    Dimb::Int64 #the rank of begin dimension in this independent helicities
    Dime::Int64 #the rank of end dimension in this independent helicities
    Dimn::Int64 #the number of dimensions in this independent helicities
    Dimo::Int64 # +1 or not for onshell W>sum m. Np+Dimo is the total dimensions of a independent helicities
end
# store the information of each dimension for matrix
mutable struct structDimension #Dim[]
    iIH::Int64 # which independent helicity this dimension belongs to 
    k::ComplexF64 #momenta of discreting
    w::Float64 #weights of discreting
    Dimo::Int64 # +1 or not for onshell W>sum m. Np+Dimo is the total dimensions of a independent helicities
end
#*******************************************************************************************
# momenta of partilces of a 2-2 interaction.
mutable struct structMomentum
    i1::SVector{5,ComplexF64} #initial 1
    f1::SVector{5,ComplexF64} #final 1
    i2::SVector{5,ComplexF64}
    f2::SVector{5,ComplexF64}
    q::SVector{5,ComplexF64} #exchange 
    q2::Complex{Float64} # q^2
    qt::Complex{Float64} #
end
# helicities of partilces of a 2-2 interaction.
mutable struct structHelicity #l
    i1::Int64 # intilal 1
    i1h::Int64 #for fermion, the helicity is obtianed by i1/i1h
    f1::Int64
    f1h::Int64
    i2::Int64
    i2h::Int64
    f2::Int64
    f2h::Int64
end
#*******************************************************************************************
#store the information of particles involved in the work. 
struct structParticle
    name::String #
    name0::String #without charge, used for the cases with symmetry.
    nameL::String #Latex.
    m::Float64 #mass
    J::Int64 #spin
    Jh::Int64  #for fermion, the helicity is obtianed by J/Jh
    P::Int64 #parity
end
# function to read the information of partilces from a file with "filename".
function particles!(particles::Vector{structParticle}, pkey::Dict{String,Int64}, filename::String) #read the information 
    open(filename, "r") do file
        readline(file)
        i = 1
        for line in eachline(file)
            parts = split(line)
            particle = structParticle(parts[1], parts[2], parts[3], parse(Float64, parts[4]),
                parse(Int64, parts[5]), parse(Int64, parts[6]), parse(Int64, parts[7]))
            push!(particles, particle)
            pkey[parts[1]] = i
            i += 1
        end
    end
end
const p = structParticle[] #store of information of particles in this global vector
const pkey = Dict{String,Int64}() #store the dictionary for key and the number of particles
#*******************************************************************************************
function ch(pf, pin, amps)
    mf = [qBSE.p[qBSE.pkey[p]].m for p in pf]
    namef = ["\\" * qBSE.p[qBSE.pkey[p]].nameL for p in pf]
    mi = [qBSE.p[qBSE.pkey[p]].m for p in pin]
    namei = ["\\" * qBSE.p[qBSE.pkey[p]].nameL for p in pin]
    ch = (pf=pf, namef=namef, mf=mf, pin=pin, namei=namei, mi=mi, amps=amps)
end
#*******************************************************************************************
# Calculate the \sum|M|^2 for all channels 
function M2_channel(T, CH, IH)

    NCH = length(CH)
    resM2 = zeros(Float64, NCH, NCH)
    for iCHf in 1:NCH
        for iCHi in 1:NCH
            for f in CH[iCHf].IHb:CH[iCHf].IHe, i in CH[iCHi].IHb:CH[iCHi].IHe
                if IH[f].Dimo == 1 && IH[i].Dimo == 1
                    resM2[iCHf, iCHi] += abs2(T[IH[f].Dime, IH[i].Dime])
                end
            end
        end
    end

    return resM2
end
@inline function lambda(m1, m2, m3)
    m1_sq = m1 * m1
    m2_plus_m3 = m2 + m3
    m2_minus_m3 = m2 - m3
    l = (m1_sq - m2_plus_m3 * m2_plus_m3) * (m1_sq - m2_minus_m3 * m2_minus_m3)
    return l > 0.0 ? sqrt(l) : 0.0
end
# Frame transition from the static frame of total system  to  cms of two partilce.
function LorentzBoostRotation(k, tecm, p1, p2)
    kij = k[p1] + k[p2]
    sij = kij * kij
    P = @SVector [0.0, 0.0, 0.0, tecm, tecm]
    pLB = @SVector [-kij[1], -kij[2], -kij[3], kij[4], sqrt(sij)]
    knew = Xs.LorentzBoost(k, pLB)
    Pnew = Xs.LorentzBoost(P, pLB)
    ct, st, cp, sp = FR.kph(knew[p2])
    Pnew = Xs.Rotation(Pnew, ct, st, cp, sp)
    knew = Xs.Rotation(knew, ct, st, cp, sp)
    return knew, Pnew
end
# Calculate cross section for all channels , only for 2-2 process
function simpleXsection(ER, M2, CH, qn; Ep="cm")
    Xsection = Matrix{Float64}[]
    NCH = length(CH)
    N = length(ER)
    cons = 0.3894 / (256.0 * pi^3)

    for i in 1:N
        Xs0 = zeros(size(M2[1]))  # Preallocate `Xs0` and overwrite it in each loop iteration instead of reallocating
        W = ER[i]

        for iM2 in 1:NCH
            pi1 = p[CH[iM2].p[1]]
            pi2 = p[CH[iM2].p[2]]
            mi1 = pi1.m
            mi2 = pi2.m
            jitilde = (2.0 * pi1.J + 1.0) * (2.0 * pi2.J + 1.0)
            fac2mi1 = (pi1.Jh == 1) ? 1 : 2.0 * mi1
            fac2mi2 = (pi2.Jh == 1) ? 1 : 2.0 * mi2


            if Ep == "L"
                W = sqrt((sqrt(W^2 + mi1^2) + mi2)^2 - W^2)
            end

            for fM2 in 1:NCH
                pf1 = p[CH[fM2].p[1]]
                pf2 = p[CH[fM2].p[2]]
                mf1 = pf1.m
                mf2 = pf2.m
                fac2mf1 = (pf1.Jh == 1) ? 1 : 2.0 * mf1
                fac2mf2 = (pf2.Jh == 1) ? 1 : 2.0 * mf2

                lam = lambda(W, mf1, mf2) / lambda(W, mi1, mi2) / W^2
                fac2m = fac2mi1 * fac2mi2 * fac2mf1 * fac2mf2
                fac = lam * fac2m * cons * (2.0 * qn.J + 1.0) / jitilde

                @inbounds Xs0[iM2, fM2] = M2[i][iM2, fM2] * fac
            end
        end

        push!(Xsection, Xs0)
    end

    return Xsection
end
# Evaluate log-determinant |1 - V G| to locate poles
function res(Range, iER, qn, SYS, IA, CH, IH, fV)
    Ect = ComplexF64[] # Array to store complex total energy Ec=ER+EI*im
    reslogt = Float64[] # Array to store log|1-VG|
    Dim = nothing
    TG = nothing
    NCH = length(CH)
    resM2 = zeros(Float64, NCH, NCH)
    ER = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / (Range.NER - 1) # Calculate ER


    EI = 0.0
    if Range.Ep == "L" #here, we use the pL, so should be transfered to ER
        PL = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / (Range.NER - 1)# Calculate ER
        ER = sqrt((sqrt(PL^2 + qBSE.p[qBSE.pkey["K_m"]].m^2) + qBSE.p[qBSE.pkey["N_p"]].m)^2 - PL^2)
    end
    NEI = Range.NEI
    for iEI in -NEI:NEI
        if NEI != 0
            EI = iEI * Range.EIt / NEI
        end
        Ec = ER + EI * im
        Vc, Gc, II, IH, Dim = qBSE.VGI(Ec, qn, SYS, IA, CH, IH, fV, lRm=qn.lRm) # Calculate the V, G, and unit matrix II
        VGI = II - Vc * Gc
        detVGI = det(VGI)   # Compute determinant of (1 - VGc)
        logdetVGI = log(abs(detVGI)^2)
        push!(Ect, Ec)
        push!(reslogt, logdetVGI)
        if iEI == 0
            T = inv(VGI) * Vc
            TG = T * Gc

            resM2 = qBSE.M2_channel(T, CH, IH)
        end
    end

    if Range.Ep == "L"
        Ect[1] = PL
    end
    return Ect, reslogt, resM2, IH, Dim, TG
end
#*******************************************************************************************
# Get the corresponding IH and Dim for the given case.
# Note: Interpolation is not performed here simultaneously, as it may cause memory leaks for unknown reasons.
function IHDim(E, Et, Range, Tt, IHt, Dimt)
    ii = Range.NER - Xs.Nsij(E, Range.ERmin, Range.ERmax, Range.NER - 1)
    Emin, Emax = Et[ii+1], Et[ii]
    Tmin, Tmax = Tt[ii+1], Tt[ii]
    if size(Tmin) == size(Tmax)
        return IHt[ii], Dimt[ii]
    else
        mid = 0.5 * (Emin + Emax)
        idx = (E < mid) + 1
        return IHt[ii+idx-1], Dimt[ii+idx-1]
    end
end
#set the frame and other things for calculating TGA
function setTGA(par, sij, k, tecm, i, j)
    IHt, Dimt, TGt, Et, Range = par.IHt, par.Dimt, par.TGt, par.Et, par.Range
    IH, Dim = qBSE.IHDim(sqrt(sij), Et, Range, TGt, IHt, Dimt)
    kn, Pn = qBSE.LorentzBoostRotation(k, tecm, i, j) #reference frame thansformation
    return (E=sqrt(sij), par=par, IH=IH, Dim=Dim, k=kn, P=Pn)
end
#calculate TGA
function TGA(cfinal, cinter, Vert, para)

    E, par, IH, Dim, k, P = para.E, para.par, para.IH, para.Dim, para.k, para.P
    Et, TGt, qn, CH, SYS = par.Et, par.TGt, par.qn, par.CH, par.SYS

    # Interpolation: placing interpolation here helps prevent memory leaks
    ii = par.Range.NER - Xs.Nsij(E, par.Range.ERmin, par.Range.ERmax, par.Range.NER - 1)
    Emin, Emax = Et[ii+1], Et[ii]
    Tmin, Tmax = TGt[ii+1], TGt[ii]
    TeT = size(Tmin) == size(Tmax)
    if TeT
        inv_dE = 1.0 / (Emax - Emin)
        ww = (E - Emin) * inv_dE

    else
        mid = 0.5 * (Emin + Emax)
        idx = (E < mid) + 1
        TG = TGt[ii+idx-1]
    end


    # Independent helicities for final states
    IHb = CH[cfinal].IHb
    IHe = CH[cfinal].IHe
    # Rank of dimensions for intermediate state
    Dimb = IH[CH[cinter].IHb].Dimb
    Dime = IH[CH[cinter].IHe].Dime
    lJJ = qn.J / qn.Jh
    NJ2d2 = (2.0 * lJJ + 1.0) / (8.0 * pi)
    TGAdic = Dict()
    ndx = length(SYS.xv)
    ndphi = length(SYS.pv)
    dx = 2.0 / ndx
    dphi = 2.0 * pi / ndphi
    #loop for all Independent helicities of final state
    @inbounds for il in IHb:IHe
        m1, m2 = CH[il].m[1], CH[il].m[2]
        IHl = IH[il]
        l21 = IHl.hel[2] / IHl.helh[2] - IHl.hel[1] / IHl.helh[1]
        dl21 = Int64(lJJ + l21) + 1
        #loop for dimensions for intermediate state
        TGA = 0.0
        eta = 0.0
        @inbounds for iDim in Dimb:Dime
            p20 = real(Dim[iDim].k)
            p20m1sq = sqrt(p20^2 + m1^2)
            p20m2sq = sqrt(p20^2 + m2^2)
            IHp = IH[Dim[iDim].iIH]
            l21p = IHp.hel[2] / IHp.helh[2] - IHp.hel[1] / IHp.helh[1]
            CH0 = CH[IHp.iCH]
            etap = CH0.P[1] * CH0.P[2] * qn.P * (-1)^(qn.J / qn.Jh - CH0.J[1] / CH0.Jh[1] - CH0.J[2] / CH0.Jh[2])
            dl21pp = Int64(lJJ + l21p) + 1
            dl21pm = Int64(lJJ - l21p) + 1
            # loop for integration of angles
            A = 0.0
            la, lb = IHp.hel[1], IHp.hel[2]
            for j in 1:ndphi
                #phi
                phi = SYS.pv[j]
                cosphi, sinphi = SYS.cp[j], SYS.sp[j]
                exppl = exp(-im * phi * l21)
                for i in 1:ndx
                    #theta
                    ct = SYS.xv[i]
                    st = sqrt(1.0 - ct^2)
                    #momentum for intermediate state. 
                    px, py, pz = p20 * st * cosphi, p20 * st * sinphi, p20 * ct
                    p1 = SVector{5,Float64}(-px, -py, -pz, p20m1sq, m1)
                    p2 = SVector{5,Float64}(px, py, pz, p20m2sq, m2)
                    Ap = Vert.Vert(p1, p2, la, lb, Vert, k, P)
                    Am = Vert.Vert(p1, p2, -la, -lb, Vert, k, P)
                    Dp = SYS.d[i][dl21, dl21pp]
                    Dm = SYS.d[i][dl21, dl21pm] * etap
                    A -= exppl * (Dp * Ap + Dm * Am) * SYS.wxv[i] * SYS.wpv[i]

                end
            end

            if TeT
                TG0 = Tmin[IH[il].Dime, iDim] * (1.0 - ww) + Tmax[IH[il].Dime, iDim] * ww #插值
            else
                TG0 = TG[IH[il].Dime, iDim]
            end

            TGA += A * TG0
        end

        h1, h2 = IHl.hel[1], IHl.hel[2]
        TGA *= NJ2d2
        TGAdic[(h1, h2)] = TGA
        # extend from independent helicities to all helicities
        CHf = CH[il]
        eta = CHf.P[1] * CHf.P[2] * qn.P * (-1)^(qn.J / qn.Jh - CHf.J[1] / CHf.Jh[1] - CHf.J[2] / CHf.Jh[2])
        TGAdic[(-h1, -h2)] = TGA * eta
        if h1 == 0 && h2 == 0
            TGAdic[(h1, h2)] = TGA * sqrt(2.0)
        end

    end
    return TGAdic
end
#*******************************************************************************************
#显示
function showSYSInfo(Range, qn, IA, CH, IH)
    dashline = repeat('-', 70)
    println(dashline)
    @printf("I(J,P)= %1d(%1d,%2d)\n", qn.I, qn.J, qn.P)
    for ih in eachindex(IH)
        @printf("%-15s: %2d/%1d, %2d/%1d;  cutoff = %5.3f GeV\n", CH[IH[ih].iCH].name[3], IH[ih].hel[1],
            IH[ih].helh[1], IH[ih].hel[2], IH[ih].helh[2], CH[IH[ih].iCH].cutoff)
    end
    Nc = length(CH)
    println(dashline)
    @printf("%-15s\n", "channel")

    for i1 in 1:Nc
        @printf("%-15s", CH[i1].name[3])
        for i in 1:Nc
            @printf("%2d ", IA[i1, i].Nex)
        end
        @printf("\n")
    end
    println(dashline)
    println("ER=$(Range.ERmax) to $(Range.ERmin) GeV, NER=$(Range.NER)  ;   EI=$(Range.EIt*1e3) MeV, NEI=$(Range.NEI)")
    println(dashline)
end
function showPoleInfo(qn, Ec, reslog, filename)
    Ampmin, Ampminx, Ampminy = 0.0, 0.0, 0.0
    for i in eachindex(Ec)
        Eci = Ec[i]
        logdetVGI = reslog[i]
        open(filename, "a") do file
            write(file, @sprintf("%.4f %.4f %.4f\n", real(Eci), imag(Eci) * 1e3, logdetVGI))
        end
        if Ampmin > logdetVGI
            Ampmin = logdetVGI
            Ampminx = real(Eci)
            Ampminy = imag(Eci) * 1e3
        end
    end


    dashline = repeat('-', 70)
    println(dashline)
    println("I(J,P)=$(qn.I)($(qn.J),$(qn.P))  pole= $Ampmin at $(Ampminx * 1e3), $Ampminy")
    println(dashline)
    return Ampmin, Ampminx, Ampminy
end
#*******************************************************************************************
#form factors
function FFre(k, cutoffi, cutofff; cutoff_re_type=:Lambda, CHi=nothing, CHf=nothing, key_ex=0)

    if cutoff_re_type == :alpha_light
        cutoffi = CHi.m[1] + 0.22 * cutoffi
        cutofff = CHf.m[1] + 0.22 * cutofff
    elseif cutoff_re_type == :alpha
        cutoffi = p[key_ex].m + 0.22 * cutoffi
        cutofff = p[key_ex].m + 0.22 * cutofff
    end

    FFre = 0.0 + 0im #form factor for constituent particles

    mi1 = real(k.i1[5])
    mf1 = real(k.f1[5])
    mi2 = real(k.i2[5])
    mf2 = real(k.f2[5])
    if mi1 <= mi2
        FFre += -(mi1^2 - k.i1 * k.i1)^2 / cutoffi^4

    end
    if mi1 > mi2
        FFre += -(mi2^2 - k.i2 * k.i2)^2 / cutoffi^4
    end
    if mf1 <= mf2
        FFre += -(mf1^2 - k.f1 * k.f1)^2 / cutofff^4
    end
    if mf1 > mf2
        FFre += -(mf2^2 - k.f2 * k.f2)^2 / cutofff^4
    end

    return exp(FFre)

end
function propFFex(k, key_ex, cutoff; cutoff_ex_type=:Lambda, FF_ex_type=3)

    m = p[key_ex].m

    if cutoff_ex_type == :alpha
        cutoff = m + 0.22 * cutoff
    end

    prop = 1.0 + 0im #propagator
    FFex = 1.0 + 0im #form factor for exchanged particles
    if p[key_ex].name != "V" #not contact interaction
        prop *= (1.0 + 0im) / (k.q2 - ComplexF64(m^2))
        if FF_ex_type == 1
            FFex *= ((cutoff^2 - m^2) / (cutoff^2 - k.q2))^2  #type 1
        elseif FF_ex_type == 2
            FFex *= (cutoff^4 / ((m^2 - k.q2)^2 + cutoff^4))^2 #type 2
        elseif FF_ex_type == 3
            FFex *= exp(-2.0 * (m^2 - k.q2)^2 / cutoff^4)
        elseif FF_ex_type == 4
            FFex *= ((cutoff^4 + (k.qt - m^2)^2 / 4) / ((k.q2 - (qt + m^2) / 2)^2 + cutoff^4))^2
        elseif FF_ex_type == 5
            FFex *= (cutoff^2 / (cutoff^2 - k.q2))^2  #type 1
        elseif FF_ex_type == 6
            FFex *= exp(-(m^2 - k.q2)^2 / cutoff^4)
        end
    end

    return prop * FFex

end

#potential kernel
function kernel(kf, ki, Ec, qn, SYS, IA, CH, IHf, IHi, fV)::ComplexF64 # Calculating kernel
    ichi = IHi.iCH # channel
    ichf = IHf.iCH
    IA0 = IA[ichf, ichi]
    CHi = CH[ichi]
    CHf = CH[ichf]
    if IA[ichi, ichf].Nex == 0 # no exchange, returen 0
        return 0 + 0im
    end

    mi1, mi2 = CHi.m[1], CHi.m[2]
    mf1, mf2 = CHf.m[1], CHf.m[2]

    Ei1, Ei2 = sqrt(ki^2 + mi1^2), sqrt(ki^2 + mi2^2)
    Ef1, Ef2 = sqrt(kf^2 + mf1^2), sqrt(kf^2 + mf2^2)

    ki10, ki20, qti = (mi1 + 1e-7 <= mi2) ? (Ec - Ei2, Ei2, mi2) : (Ei1, Ec - Ei1, Ec - mi1)
    kf10, kf20, qtf = (mf1 + 1e-7 <= mf2) ? (Ec - Ef2, Ef2, mf2) : (Ef1, Ec - Ef1, Ec - mf1)

    qt = (qti - qtf)^2 + 0im
    
    l = structHelicity(  #helicities
        IHi.hel[1], IHi.helh[1],
        IHf.hel[1], IHf.helh[1],
        IHi.hel[2], IHi.helh[2],
        IHf.hel[2], IHf.helh[2]
    )

    eta = CHi.P[1] * CHi.P[2] * qn.P *
          (-1)^(qn.J / qn.Jh - CHi.J[1] / CHi.Jh[1] - CHi.J[2] / CHi.Jh[2])

    lJJ = qn.J / qn.Jh
    l21i = -l.i2 / l.i2h + l.i1 / l.i1h
    l21f = l.f2 / l.f2h - l.f1 / l.f1h
    lf = Int64(lJJ + l21f) + 1
    ndx = length(SYS.xv)
    Ker0 = 0 + 0im
    for i in 1:ndx
        x = SYS.xv[i]
        sqrt1_x2 = sqrt(1 - x^2)

        ki1 = @SVector [0.0 + 0im, 0.0 + 0im, -ki + 0im, ki10 + 0im, mi1 + 0im]
        kf1 = @SVector [-kf * sqrt1_x2 + 0im, 0.0 + 0im, -kf * x + 0im, kf10 + 0im, mf1 + 0im]
        ki2 = @SVector [0.0 + 0im, 0.0 + 0im, ki + 0im, ki20 + 0im, mi2 + 0im]
        kf2 = @SVector [kf * sqrt1_x2 + 0im, 0.0 + 0im, kf * x + 0im, kf20 + 0im, mf2 + 0im]

        q = kf2 - ki2
        q2 = q * q
        k = structMomentum(ki1, kf1, ki2, kf2, q, q2, qt)


        l.f1, l.i2 = -l.f1, -l.i2  #helicity to spin and the minus for fixed parity
        Kerm1 = fV(k, l, SYS, IA0, CHf, CHi) * SYS.d[i][Int64(lJJ + l21i)+1, lf] * eta
        l.f1, l.i2 = -l.f1, -l.i2
        l.i1, l.f1 = -l.i1, -l.f1
        Kerp1 = fV(k, l, SYS, IA0, CHf, CHi) * SYS.d[i][Int64(lJJ - l21i)+1, lf]
        l.i1, l.f1 = -l.i1, -l.f1
        Ker0 += (Kerp1 + Kerm1) * SYS.wxv[i]
        #@show i,Ker0
    end

    kernel = Ker0 * 2pi
    if l.f1 == 0 && l.f2 == 0  ## factors from fixed parity 
        kernel /= sqrt(2.0)
    end
    if l.i1 == 0 && l.i2 == 0
        kernel /= sqrt(2.0)
    end

    return kernel
end
#*******************************************************************************************
#propagator
function propagator(k, kv, w, wv, Ec, Np, CH, IH0, Dimo, lRm)
    propagator = Complex{Float64}(0, 0)
    ich = IH0.iCH
    mi1, mi2 = CH[ich].m[1], CH[ich].m[2]
    lRm0 = isa(lRm, Int) ? lRm :
           isa(lRm, Tuple{Vararg{Int}}) ? lRm[ich] :
           error("lRm should be Int or Tuple{Vararg{Int}}")
    Rm = lRm0 == 1 ? -1 : 1

    # Compute mi2p2 based on helicity conditions
    mi2p2 = 1.0
    mi2p2 *= (IH0.helh[1] == 2) ? 2.0 * mi1 : 1.0
    mi2p2 *= (IH0.helh[2] == 2) ? 2.0 * mi2 : 1.0

    # Precompute constants
    pi_factor = (2 * pi)^3
    complex_factor = Complex{Float64}(0, 1) / (32 * pi^2 * Ec)

    # For off-momenta
    if Dimo == 0
        cE1, cE2 = sqrt(k^2 + mi1^2), sqrt(k^2 + mi2^2)
        denominator = (mi1 <= mi2) ? 2 * cE2 * ((Ec - cE2)^2 - cE1^2) : 2 * cE1 * ((Ec - cE1)^2 - cE2^2)

        if denominator != 0  # To avoid division by zero
            propagator = mi2p2 * k^2 * w / pi_factor / denominator
        end
    end

    # For on-shell
    if Dimo == 1
        delta1 = Ec^2 - (mi1 + mi2)^2
        delta2 = Ec^2 - (mi2 - mi1)^2
        konc = sqrt(delta1 * delta2) / (2 * Ec)
        konc2 = delta1 * delta2 / (4 * Ec^2)

        # Add imaginary component  and ensure correct sign
        # + for II Riemann sheet, the sign for I sheet is added by Rm. 
        if imag(konc) >= 0
            propagator += mi2p2 * konc * complex_factor * Rm
        else
            propagator -= mi2p2 * konc * complex_factor * Rm #conjugate 
        end

        # Loop over Np elements
        for i3 in 1:Np
            kv_i3 = kv[i3]
            denominator = 2 * Ec * (kv_i3^2 - konc2)
            propagator += mi2p2 * wv[i3] * konc2 / pi_factor / denominator
        end
    end

    return propagator
end

#*******************************************************************************************


# calculate V G I
function VGI(Ec, qn, SYS, IA, CH, IH, fV; lRm=0)
    # Determine the dimension of the work matrix. 
    Np = length(SYS.kv)
    SYS, IH, Dim, = workSpace(Ec, lRm, Np, SYS, CH, IH)
    Nt = length(Dim)

    # Initialize matrices
    Vc = zeros(Complex{Float64}, Nt, Nt)
    Gc = zeros(Complex{Float64}, Nt, Nt)
    II = zeros(Float64, Nt, Nt)

    # Populate Vc, Gc, and II
    for i_f = 1:Nt
        Dimf = Dim[i_f]
        kf = Dimf.k

        for i_i = i_f:Nt
            Dimi = Dim[i_i]
            ki = Dimi.k

            # Compute Vc for upper triangular part
            Vc[i_f, i_i] = kernel(kf, ki, Ec, qn, SYS, IA, CH, IH[Dimf.iIH], IH[Dimi.iIH], fV)::ComplexF64
            #@show Vc[i_f, i_i] 
            #exit()

            # Compute Gc and II only on the diagonal
            if i_i == i_f
                Gc[i_f, i_i] = propagator(kf, SYS.kv, Dimf.w, SYS.wv, Ec, Np, CH, IH[Dimf.iIH], Dimf.Dimo, lRm)
                II[i_f, i_i] = 1.0
            end
        end
    end

    # Populate lower triangular part of Vc based on conjugate symmetry
    for i_f = 2:Nt
        Dimf_IH = IH[Dim[i_f].iIH]

        for i_i = 1:i_f-1
            Dimi_IH = IH[Dim[i_i].iIH]

            if Dimf_IH.iCH != Dimi_IH.iCH
                Vc[i_f, i_i] = conj(Vc[i_i, i_f])
            else
                Vc[i_f, i_i] = Vc[i_i, i_f]
            end
        end
    end
    #@show Vc[:,1]
    #exit()
    return Vc, Gc, II, IH, Dim
end

function workSpace(Ec, lRm, Np, SYS, CH, IH)
    #produce the dimensions in a independent helicity amplitudes, especially for the onshell dimension.
    Ec += Complex{Float64}(0.0, 1e-15)
    Dim = structDimension[]
    Nih = length(IH)
    Nt = 0
    for i1 in 1:Nih
        IH[i1].Dimo = 0
        IH[i1].Dimb = 0
        IH[i1].Dime = 0
        IH[i1].Dimn = 0
    end
    for i1 in 1:Nih
        IH[i1].Dimb = Nt + 1
        mass1 = p[CH[IH[i1].iCH].p[1]].m
        mass2 = p[CH[IH[i1].iCH].p[2]].m

        lRm0 = if isa(lRm, Int)
            lRm
        elseif isa(lRm, Tuple{Vararg{Int}})
            lRm[IH[i1].iCH]
        else
            error("lRm should be Int or Tuple of Int")
        end

        temp = real(Ec) < (mass1 + mass2)
        Rm = !((lRm0 == 0 || lRm0 == 1) && temp)

        if Rm
            kon = sqrt((Ec^2 - (mass1 + mass2)^2) * (Ec^2 - (mass1 - mass2)^2)) / (2.0 * Ec)
            IH[i1].Dimo = 1
            Nt += Np + IH[i1].Dimo
            IH[i1].Dime = Nt
            IH[i1].Dimn = Np + IH[i1].Dimo
            for idim in 1:Np
                Dim0 = structDimension(i1, SYS.kv[idim], SYS.wv[idim], 0)
                push!(Dim, Dim0)
            end
            Dim0 = structDimension(i1, kon, 0.0, 1)
            push!(Dim, Dim0)
        else
            Nt += Np
            IH[i1].Dime = Nt
            IH[i1].Dimn = Np

            for idim in 1:Np
                Dim0 = structDimension(i1, SYS.kv[idim], SYS.wv[idim], 0)
                push!(Dim, Dim0)

            end
        end

    end


    return SYS, IH, Dim
end

#*******************************************************************************************
function preprocessing(Sys, channels, Ff, qn, cutoff; Np=10, Nx=5, Nphi=5)
    #store the channel information, interaction information in CH,IH,IA
    CH = structChannel[]
    IH = structIndependentHelicity[]
    Nih, Nc = 1, 1

    for channel in channels
        p1, p2 = pkey[channel[1]], pkey[channel[2]]
        m1, m2 = p[p1].m, p[p2].m
        J1, J2 = p[p1].J, p[p2].J
        Jh1, Jh2 = p[p1].Jh, p[p2].Jh
        P1, P2 = p[p1].P, p[p2].P
        IHn = 0
        #independent helicities
        Nih0 = Nih
        IH0 = structIndependentHelicity[]
        IH00 = structIndependentHelicity(0, Int64[], Int64[], 0, 0, 0, 0)
        for i1 in -p[p1].J:p[p1].Jh:p[p1].J
            for i2 in -p[p2].J:p[p2].Jh:p[p2].J

                Jh = 1.0
                if p[p1].Jh == 1 && p[p2].Jh == 2 || p[p1].Jh == 2 && p[p2].Jh == 1
                    Jh = 2.0
                end
                if abs(float(i1) / float(p[p1].Jh) - float(i2) / float(p[p2].Jh)) <= float(qn.J) / float(qn.Jh) + 0.01

                    lind = 1
                    for i3 in 1:Nih-Nih0 #Nih0:(Nih-1)
                        if IH0[i3].hel[1] == -i1 && IH0[i3].hel[2] == -i2
                            lind = 0
                            break
                        end
                    end
                    if lind == 1
                        IH00.hel = [i1, i2]

                        IH00.helh = [p[p1].Jh, p[p2].Jh]
                        IH00.iCH = Nc
                        push!(IH0, deepcopy(IH00))

                        Nih += 1
                        IHn += 1
                    end
                end
            end
        end
        Nc += 1
        #channel information
        CH0 = structChannel(
            [p1, p2],
            [m1, m2],
            [J1, J2],
            [Jh1, Jh2],
            [P1, P2],
            channel[3],
            [p[p1].name, p[p2].name, "$(p[p1].name):$(p[p2].name)"],
            [p[p1].name0, p[p2].name0, "$(p[p1].name0):$(p[p2].name0)"],
            Nih0,
            Nih0 + IHn - 1,
            IHn
        )
        push!(CH, deepcopy(CH0))
        append!(IH, deepcopy(IH0))
    end

    Nc -= 1
    Nih -= 1

    # interaction information.
    chnamei = Vector{String}(undef, 1)
    chnamef = Vector{String}(undef, 1)
    IA = Matrix{structInterAction}(undef, Nc, Nc)
    for i1 in 1:Nc
        for i2 in 1:Nc
            chnamei[1] = CH[i1].name[3]
            chnamef[1] = CH[i2].name[3]
            chname = "$(chnamei[1])-->$(chnamef[1])"

            if i1 <= i2
                if haskey(Ff, chname) == true
                    Nex = length(Ff[chname])
                    IA0 = structInterAction(
                        Nex,
                        [Ff[chname][i][1][1] for i in 1:Nex],  #name
                        [pkey[Ff[chname][i][1][1]] for i in 1:Nex],  #key
                        [p[pkey[Ff[chname][i][1][1]]].J for i in 1:Nex],  #J
                        [p[pkey[Ff[chname][i][1][1]]].Jh for i in 1:Nex],  #Jh
                        [p[pkey[Ff[chname][i][1][1]]].P for i in 1:Nex],  #P
                        [p[pkey[Ff[chname][i][1][1]]].m for i in 1:Nex],  #m
                        [Ff[chname][i][1][2] for i in 1:Nex],  #dc
                        [Ff[chname][i][2] for i in 1:Nex]     #Ff
                    )
                else
                    Nex = 0
                    IA0 = structInterAction(
                        Nex,
                        [],  #name
                        [],  #key
                        [],  #J
                        [],  #Jh
                        [],  #P
                        [],  #m
                        [],  #dc
                        []     #Ff
                    )
                end
                IA[i1, i2] = IA0
                IA[i2, i1] = IA0
            end

        end
    end


    # generate Gauss–Legendre nodes and weights in[-1, 1]  
    nodes, weights = gausslegendre(Np)
    # to [0, π/2]
    k0 = 0.5 * (nodes .+ 1.0) * (π / 2)
    w0 = 0.5 * (π / 2) * weights
    kv = tan.(k0)
    wv = w0 ./ (cos.(k0) .^ 2)

    #kv, wv = gausslaguerre(Np)
    #wv = wv .* exp.(kv)

    xv, wxv = gausslegendre(Nx)

    #xv = -1 .+ 2 / Nx .* ((1:Nx) .- 0.5)
    #wxv = fill(2 / Nx, Nx)


    wd = [wignerd(qn.J / qn.Jh, acos(xv[i])) for i in 1:Nx]

    pv, wpv = gausslegendre(Nphi)
    pv = pi * (pv .+ 1.0)   # to [0, 2π]
    wpv = wpv * pi


    sp = [sin(pv[i]) for i in 1:Nphi]
    cp = [cos(pv[i]) for i in 1:Nphi]

    SYS = structSys(Sys, kv, wv, xv, wxv, wd, pv, wpv, sp, cp,
        get(cutoff, :cutoff_re_type, :Lambda),
        get(cutoff, :cutoff_ex, 0.0),
        get(cutoff, :cutoff_ex_type, :Lambda),
        get(cutoff, :FF_ex_type, 3))

    return SYS, IA, CH, IH
end

# Gauss-Legendre quadrature rule generator 
function cgqf(N::Int, a::Float64, b::Float64)

    J = zeros(N, N)
    for i in 1:N-1
        v = i / sqrt(4i^2 - 1)
        J[i, i+1] = v
        J[i+1, i] = v
    end

    vals, vecs = eigen(J)
    x_std = vals
    w_std = 2 * vecs[1, :] .^ 2

    x = 0.5 * ((b - a) * x_std .+ (b + a))
    w = 0.5 * (b - a) * w_std

    return x, w
end


end
