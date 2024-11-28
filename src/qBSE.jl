module qBSE
using FastGaussQuadrature, StaticArrays, WignerD, LinearAlgebra, Printf
using Base.Threads
using StaticArrays
using ..Xs, ..FR
#*******************************************************************************************
#store informations about system 
struct structSys
    Sys::String
    kv::Vector{Float64} #momenta of discreting
    wv::Vector{Float64} #weights of discreting
    d::Vector{Matrix{Float64}}   #Wignerd
    sp::Vector{Float64}
    cp::Vector{Float64}
end
# store the information of a interaction
struct structInterAction #IA
    Nex::Int64  #total number of exchanges
    name_ex::Vector{String}  #exchange
    key_ex::Vector{Int64}
    J_ex::Vector{Int64}
    Jh_ex::Vector{Int64}
    P_ex::Vector{Int64}
    m_ex::Vector{Float64}  #exchange
    dc::Vector{Int64}  #direct or cross
    CC::Vector{Float64}  #flavor factors
end
# store the information of a channel. The kv and wv are stored here for convenience. 
struct structChannel #CH
    p::Vector{Int64} #particles of channel
    m::Vector{Float64} #particle masses of channel
    J::Vector{Int64} #particle spins of channel
    Jh::Vector{Int64} #particle  spins of channel
    P::Vector{Int64} #particle parities of channel
    cutoff::Float64 #cutoff of channel
    name::Vector{String}  #name of channel
    name0::Vector{String} #name of channel (without charge)
    IHb::Int64
    IHe::Int64
    IHn::Int64 #total number of helicity amplitudes of a channel
end
# store the information of an independent helicity for matrix
mutable struct structIndependentHelicity #IH
    iCH::Int64 # which channel this independent helicity belongs to 
    hel::Vector{Int64} # helicities
    helh::Vector{Int64} # fermion or boson
    Dimb::Int64 #the number of begin dimension in this independent helicities
    Dime::Int64 #the number of end dimension in this independent helicities
    Dimn::Int64 #the number of dimensions in this independent helicities
    Dimo::Int64 # +1 or not for onshell W>sum m. Np+Dimo is the total dimensions of a independent helicities
end
# store the information of each dimension for matrix
mutable struct structDimension #Dim
    iIH::Int64 # which independent helicity this dimension belongs to 
    k::ComplexF64 #momenta of discreting
    w::Float64 #weights of discreting
    Dimo::Int64 # +1 or not for onshell W>sum m. Np+Dimo is the total dimensions of a independent helicities
end
#*******************************************************************************************
mutable struct structMomentum #momenta
    i1::SVector{5,ComplexF64} # initial 1
    f1::SVector{5,ComplexF64}
    i2::SVector{5,ComplexF64}
    f2::SVector{5,ComplexF64}
    q::SVector{5,ComplexF64} #exchange
    q2::Complex{Float64} # q^2
    qt::Complex{Float64} #
end
mutable struct structHelicity #l
    i1::Int64
    i1h::Int64
    f1::Int64
    f1h::Int64
    i2::Int64
    i2h::Int64
    f2::Int64
    f2h::Int64
end
#*******************************************************************************************
struct structParticle #store the information of particles
    name::String
    name0::String
    m::Float64
    J::Int64
    Jh::Int64
    P::Int64
end
function particles!(particles::Vector{structParticle}, pkey::Dict{String,Int64}, filename::String) #read the information 
    open(filename, "r") do file
        readline(file)
        i = 1
        for line in eachline(file)
            parts = split(line)
            particle = structParticle(parts[1], parts[2], parse(Float64, parts[3]),
                parse(Int64, parts[4]), parse(Int64, parts[5]), parse(Int64, parts[6]))
            push!(particles, particle)
            pkey[parts[1]] = i
            i += 1
        end
    end
end
const p = structParticle[] #store of information of particles in this global vector
const pkey = Dict{String,Int64}() #store the dictionary for key and the number of particles
#*******************************************************************************************
# Calculate the \sum|M|^2 for each channel 
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
function simpleXsection(xx, M2, CH; Ep="cm")
    Xsection = Matrix{Float64}[]
    NCH = length(CH)
    N = length(xx)
    const_factor = 0.3894 / (256.0 * pi^3)  # 预先计算不变的常数部分

    for i in 1:N
        Xs0 = zeros(size(M2[1]))  # 提前分配 `Xs0`，每次循环覆盖而非重新分配
        W = xx[i]

        for iM2 in 1:NCH
            mi1 = p[CH[iM2].p[1]].m
            mi2 = p[CH[iM2].p[2]].m

            # 检查 Ep 是否为 "L"，只在需要时重新计算 W
            if Ep == "L"
                W = sqrt((sqrt(W^2 + mi1^2) + mi2)^2 - W^2)
            end

            for fM2 in 1:NCH
                mf1 = p[CH[fM2].p[1]].m
                mf2 = p[CH[fM2].p[2]].m
                lambda_factor = lambda(W, mf1, mf2) / lambda(W, mi1, mi2) / W^2
                fac = lambda_factor * (2.0 * mi2)^2 * const_factor

                @inbounds Xs0[iM2, fM2] = M2[i][iM2, fM2] * fac
            end
        end

        push!(Xsection, Xs0)  # 将计算结果添加到 `Xsection` 中
    end

    return Xsection
end
#得到log|1-VG|以寻找极点。
function res(Range, iER, qn, SYS, IA, CH, IH, fV)
    resEct = ComplexF64[] #设置保存复的总能量Ec=ER+EI*im的数组
    reslogt = Float64[] #设置保存log|1-VG|的数组
    Dim = nothing
    TG = nothing
    NCH = length(CH)
    resM2 = zeros(Float64, NCH, NCH)
    ER = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / Range.NER #计算ER
    EI = 0.0
    if Range.Ep == "L" #here, we use the pL, so should be transfered to ER
        PL = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / Range.NER #计算ER
        ER = sqrt((sqrt(PL^2 + qBSE.p[qBSE.pkey["K_m"]].m^2) + qBSE.p[qBSE.pkey["N_p"]].m)^2 - PL^2)
    end
    NEI = Range.NEI
    for iEI in -NEI:NEI #虚部部分循环
        if NEI != 0
            EI = iEI * Range.EIt / NEI #计算虚部
        end
        Ec = ER + EI * im
        Vc, Gc, II, IH, Dim = qBSE.srAB(Ec, qn, SYS, IA, CH, IH, fV, lRm=qn.lRm) # Calculate the V, G, and unit matrix II
        VGI = II - Vc * Gc
        detVGI = det(VGI)   # Compute determinant of (1 - VGc)，调用LinearAlgebra包det函数
        logdetVGI = log(abs(detVGI)^2)
        push!(resEct, Ec)
        push!(reslogt, logdetVGI)
        if iEI == 0
            T = inv(VGI) * Vc
            TG = T * Gc

            resM2 = qBSE.M2_channel(T, CH, IH)
        end
    end

    if Range.Ep == "L"
        resEct[1] = PL
    end
    return resEct, reslogt, resM2, IH, Dim, TG
end
function interpolation(E, Et, Range, Tt, IHt, Dimt)
    ii = length(Et) - Xs.Nsij(E, Range.ERmin, Range.ERmax, length(Et) - 1)
    Emin, Emax = Et[ii+1], Et[ii]
    Tmin, Tmax = Tt[ii], Tt[ii+1]
    T, IH, Dim = nothing, nothing, nothing

    if size(Tmin) == size(Tmax)
        ww = (E - Emin) / (Emax - Emin)
        T = (1.0 - ww) * Tmin + ww * Tmax
        IH = IHt[ii]
        Dim = Dimt[ii]
    else
        mid = (Emin + Emax) / 2.0
        if E >= mid
            T = Tt[ii]
            IH = IHt[ii]
            Dim = Dimt[ii]
        else
            T = Tt[ii+1]
            IH = IHt[ii+1]
            Dim = Dimt[ii+1]
        end
    end

    return T, IH, Dim
end
function TGA(lfinal, linter, Ver, predet, para)
    TG, qn, SYS, CH, IH, Dim, k, P = para.TG, para.qn, para.SYS, para.CH, para.IH, para.Dim, para.k, para.P
    IHb = CH[lfinal].IHb
    IHe = CH[lfinal].IHe
    Dimb = IH[CH[linter].IHb].Dimb
    Dime = IH[CH[linter].IHe].Dime
    lJJ = qn.J / qn.Jh
    NJ2 = (2.0 * lJJ + 1.0) / (4.0 * pi)
    TGAdic = Dict()
    for il in IHb:IHe
        IHl = IH[il]
        l21 = IHl.hel[2] / IHl.helh[2] - IHl.hel[1] / IHl.helh[1]
        Dl21 = Int64(lJJ + l21) + 1
        TGA = 0.0
        eta=0.
        for iDim in Dimb:Dime
            p20 = Dim[iDim].k
            IHp = IH[Dim[iDim].iIH]
            l21p = IHp.hel[2] / IHp.helh[2] - IHp.hel[1] / IHp.helh[1]
            CH0 = CH[IHp.iCH]
            etap = CH0.P[1] * CH0.P[2] * qn.P * (-1)^(qn.J / qn.Jh - CH0.J[1] / CH0.Jh[1] - CH0.J[2] / CH0.Jh[2])
            Dl21pp = Int64(lJJ + l21p) + 1
            Dl21pm = Int64(lJJ - l21p) + 1
            ndx = length(SYS.d)
            ndphi = length(SYS.cp)
            dx = 2.0 / ndx
            dphi = 2.0 * pi / ndphi
            A = 0.0
            for i in 1:ndx
                ct = -1.0 + dx * (i - 0.5)
                st = sqrt(1.0 - ct^2)
                for j in 1:ndphi
                    phi = (j - 1) * dphi
                    cosphi, sinphi = SYS.cp[j], SYS.sp[j]
                    px, py, pz = p20 * st * cosphi, p20 * st * sinphi, p20 * ct
                    m1, m2 = CH[il].m[1], CH[il].m[2]
                    p1 = @SVector [-px, -py, -pz, sqrt(p20^2 + m1^2), m1]
                    p2 = @SVector [px, py, pz, sqrt(p20^2 + m2^2), m2]

                    l1, l2 = IHl.hel[1], IHl.hel[2]
                    Ap = Ver(p1, p2, l1, l2, predet, k, P)
                    Am = Ver(p1, p2, -l1, -l2, predet, k, P)

                    Dp = SYS.d[i][Dl21, Dl21pp]
                    Dm = SYS.d[i][Dl21, Dl21pm] * etap

                    A -= exp(-im * phi * l21) * (Dp * Ap + Dm * Am) * dx * dphi
                    #A -= exp( im*phi/3. ) * dphi
                end
            end

            TGA += TG[IH[il].Dime, iDim] * A
        end

        h1, h2 = IHl.hel[1], IHl.hel[2]
        TGA*=NJ2
        TGAdic[(h1, h2)] = TGA
        CHf = CH[il]
        eta = CHf.P[1] * CHf.P[2] * qn.P * (-1)^(qn.J / qn.Jh - CHf.J[1] / CHf.Jh[1] - CHf.J[2] / CHf.Jh[2])
        TGAdic[(-h1, -h2)] = TGA *eta
        if h1 == 0 && h2 == 0
            TGAdic[(h1, h2)] = TGA * sqrt(2.0)
        end

    end
    return TGAdic
end
function LorentzBoostRotation(tecm, k, p1, p2)
    kij = k[p1] + k[p2]
    s12 = kij * kij
    P = @SVector [0.0, 0.0, 0.0, tecm, tecm]
    pLB = @SVector [-kij[1], -kij[2], -kij[3], kij[4], sqrt(s12)]
    knew = Xs.LorentzBoost(k, pLB)
    Pnew = Xs.LorentzBoost(P, pLB)
    ct, st, cp, sp = FR.kph(knew[p2])
    Pnew = Xs.Rotation(Pnew, ct, st, cp, sp)
    knew = Xs.Rotation(knew, ct, st, cp, sp)
    return knew, Pnew
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
function showPoleInfo(qn, resEc, reslog, filename)
    Ampmin, Ampminx, Ampminy = 0.0, 0.0, 0.0
    for i in eachindex(resEc)
        Ec = resEc[i]
        logdetVGI = reslog[i]
        open(filename, "a") do file
            write(file, @sprintf("%.4f %.4f %.4f\n", real(Ec), imag(Ec) * 1e3, logdetVGI))
        end
        if Ampmin > logdetVGI
            Ampmin = logdetVGI
            Ampminx = real(Ec)
            Ampminy = imag(Ec) * 1e3
        end
    end


    dashline = repeat('-', 70)
    println(dashline)
    println("I(J,P)=$(qn.I)($(qn.J),$(qn.P))  pole= $Ampmin at $(Ampminx * 1e3), $Ampminy")
    println(dashline)
    return Ampmin, Ampminx, Ampminy
end
#*******************************************************************************************
function fPropFF(k, ex, L, LLi, LLf; lregu=1, lFFex=0) #form factors
    m = p[ex].m
    if lFFex >= 10
        L = m + 0.22 * L
        LLi = m + 0.22 * LLi
        LLf = m + 0.22 * LLf
        lFFex -= 10
    end

    fProp = 1.0 + 0im #propagator
    FFex = 1.0 + 0im #form factor for exchanged particles
    FFre = 0.0 + 0im #form factor for constituent particles
    if p[ex].name != "V" #not contact interaction
        fProp *= (1.0 + 0im) / (k.q2 - ComplexF64(m^2))
        if lFFex == 1
            FFex *= ((L^2 - m^2) / (L^2 - k.q2))^2  #type 1
        elseif lFFex == 2
            FFex *= (L^4 / ((m^2 - k.q2)^2 + L^4))^2 #type 2
        elseif lFFex == 3
            FFex *= exp(-(m^2 - k.q2)^2 / L^4)
        elseif lFFex == 4
            FFex *= ((L^4 + (k.qt - m^2)^2 / 4) / ((k.q2 - (qt + m^2) / 2)^2 + L^4))^2
        elseif lFFex == 5
            FFex *= exp(-(m^2 - k.q2)^2 / L^4)^2
        end
    end

    if lregu == 1
        mi1 = real(k.i1[5])
        mf1 = real(k.f1[5])
        mi2 = real(k.i2[5])
        mf2 = real(k.f2[5])
        if mi1 <= mi2
            FFre += -(mi1^2 - k.i1 * k.i1)^2 / LLi^4
        end
        if mi1 > mi2
            FFre += -(mi2^2 - k.i2 * k.i2)^2 / LLi^4
        end
        if mf1 <= mf2
            FFre += -(mf1^2 - k.f1 * k.f1)^2 / LLf^4
        end
        if mf1 > mf2
            FFre += -(mf2^2 - k.f2 * k.f2)^2 / LLf^4
        end
    end
    return fProp * FFex * exp(FFre)

end
function fKernel(kf, ki, Ec, qn, SYS, IA, CH, IHf, IHi, fV)::ComplexF64 # Calculating fKernel
    ichi = IHi.iCH # channel
    ichf = IHf.iCH

    if IA[ichi, ichf].Nex == 0 # no exchange, returen 0
        return 0 + 0im
    end

    mi1, mi2 = CH[ichi].m[1], CH[ichi].m[2]
    mf1, mf2 = CH[ichf].m[1], CH[ichf].m[2]

    Ei1, Ei2 = sqrt(ki^2 + mi1^2), sqrt(ki^2 + mi2^2)
    Ef1, Ef2 = sqrt(kf^2 + mf1^2), sqrt(kf^2 + mf2^2)

    ki10, ki20, qti = (mi1 + 1e-7 <= mi2) ? (Ec - Ei2, Ei2, mi2) : (Ei1, Ec - Ei1, Ec - mi1)
    kf10, kf20, qtf = (mf1 + 1e-7 <= mf2) ? (Ec - Ef2, Ef2, mf2) : (Ef1, Ec - Ef1, Ec - mf1)
    # @show ki10, ki20, qti 
    qt = (qti - qtf)^2 + 0im

    l = structHelicity(  #helicities
        IHi.hel[1], IHi.helh[1],
        IHf.hel[1], IHf.helh[1],
        IHi.hel[2], IHi.helh[2],
        IHf.hel[2], IHf.helh[2]
    )

    eta = CH[ichi].P[1] * CH[ichi].P[2] * qn.P *
          (-1)^(qn.J / qn.Jh - CH[ichi].J[1] / CH[ichi].Jh[1] - CH[ichi].J[2] / CH[ichi].Jh[2])

    lJJ = qn.J / qn.Jh
    l21i = -l.i2 / l.i2h + l.i1 / l.i1h
    l21f = l.f2 / l.f2h - l.f1 / l.f1h
    lf = Int64(lJJ + l21f) + 1
    ndx = length(SYS.d)
    dx = 2.0 / ndx
    Ker0 = 0 + 0im
    for i in 1:ndx
        x = -1.0 + dx * (i - 0.5)
        sqrt1_x2 = sqrt(1 - x^2)
        ki1 = @SVector [0.0 + 0im, 0.0 + 0im, -ki + 0im, ki10 + 0im, mi1 + 0im]
        kf1 = @SVector [-kf * sqrt1_x2 + 0im, 0.0 + 0im, -kf * x + 0im, kf10 + 0im, mf1 + 0im]
        ki2 = @SVector [0.0 + 0im, 0.0 + 0im, ki + 0im, ki20 + 0im, mi2 + 0im]
        kf2 = @SVector [kf * sqrt1_x2 + 0im, 0.0 + 0im, kf * x + 0im, kf20 + 0im, mf2 + 0im]
        q = kf2 - ki2
        q2 = q * q
        k = structMomentum(ki1, kf1, ki2, kf2, q, q2, qt)

        l.f1, l.i2 = -l.f1, -l.i2  #helicity to spin and the minus for fixed parity
        Kerm1 = fV(k, l, SYS, IA, CH, ichi, ichf) * SYS.d[i][Int64(lJJ + l21i)+1, lf] * eta
        l.f1, l.i2 = -l.f1, -l.i2
        l.i1, l.f1 = -l.i1, -l.f1
        Kerp1 = fV(k, l, SYS, IA, CH, ichi, ichf) * SYS.d[i][Int64(lJJ - l21i)+1, lf]
        l.i1, l.f1 = -l.i1, -l.f1
        Ker0 += (Kerp1 + Kerm1) * dx
    end

    fKernel = Ker0 * 2pi
    if l.f1 == 0 && l.f2 == 0  ## factors from fixed parity 
        fKernel /= sqrt(2.0)
    end
    if l.i1 == 0 && l.i2 == 0
        fKernel /= sqrt(2.0)
    end
    #@show fKernel
    #exit()
    return fKernel
end
#*******************************************************************************************
function fProp(k, kv, w, wv, Ec, Np, CH, IH, Dimo)
    fprop = Complex{Float64}(0, 0)
    ich = IH.iCH
    mi1, mi2 = CH[ich].m[1], CH[ich].m[2]

    # Compute mi2p2 based on helicity conditions
    mi2p2 = 1.0
    mi2p2 *= (IH.helh[1] == 2) ? 2.0 * mi1 : 1.0
    mi2p2 *= (IH.helh[2] == 2) ? 2.0 * mi2 : 1.0

    # Precompute constants
    pi_factor = (2 * pi)^3
    complex_factor = Complex{Float64}(0, 1) / (32 * pi^2 * Ec)

    # For off-momenta
    if Dimo == 0
        cE1, cE2 = sqrt(k^2 + mi1^2), sqrt(k^2 + mi2^2)
        denominator = (mi1 <= mi2) ? 2 * cE2 * ((Ec - cE2)^2 - cE1^2) : 2 * cE1 * ((Ec - cE1)^2 - cE2^2)

        if denominator != 0  # To avoid division by zero
            fprop = mi2p2 * k^2 * w / pi_factor / denominator
        end
    end

    # For on-shell
    if Dimo == 1
        delta1 = Ec^2 - (mi1 + mi2)^2
        delta2 = Ec^2 - (mi2 - mi1)^2
        konc = sqrt(delta1 * delta2) / (2 * Ec)
        konc2 = delta1 * delta2 / (4 * Ec^2)

        # Add imaginary component based on `konc` sign
        if imag(konc) > 0
            fprop += mi2p2 * konc * complex_factor
        else
            fprop -= mi2p2 * konc * complex_factor
        end

        # Loop over Np elements
        for i3 in 1:Np
            kv_i3 = kv[i3]
            denominator = 2 * Ec * (kv_i3^2 - konc2)
            fprop += mi2p2 * wv[i3] * konc2 / pi_factor / denominator
        end
    end

    return fprop
end
function srAB(Ec, qn, SYS, IA, CH, IH, fV; lRm=1)
    # Determine the dimension of the work matrix. 
    Np = length(SYS.kv)
    SYS, IH, Dim, = WORKSPACE(Ec, lRm, Np, SYS, CH, IH)
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
            Vc[i_f, i_i] = fKernel(kf, ki, Ec, qn, SYS, IA, CH, IH[Dimf.iIH], IH[Dimi.iIH], fV)::ComplexF64

            # Compute Gc and II only on the diagonal
            if i_i == i_f
                Gc[i_f, i_i] = fProp(kf, SYS.kv, Dimf.w, SYS.wv, Ec, Np, CH, IH[Dimf.iIH], Dimf.Dimo)
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

    return Vc, Gc, II, IH, Dim
end
#*******************************************************************************************
function WORKSPACE(Ec, lRm, Np, SYS, CH, IH)
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
        Rm = false
        if lRm == 1 #for first Reimann sheet
            Rm = real(Ec) > (mass1 + mass2)
        elseif lRm == 2
            Rm = true
        end
        if Rm
            kon = sqrt((Ec^2 - (mass1 + mass2)^2) * (Ec^2 - (mass1 - mass2)^2)) / (2.0 * Ec)
            IH[i1].Dimo = 1
            Nt += Np + IH[i1].Dimo
            IH[i1].Dime = Nt
            IH[i1].Dimn = Np + IH[i1].Dimo
            for idim in 1:Np
                Dim0 = structDimension(IH[i1].iCH, SYS.kv[idim], SYS.wv[idim], 0)
                push!(Dim, Dim0)
            end
            Dim0 = structDimension(IH[i1].iCH, kon, 0.0, 1)
            push!(Dim, Dim0)
        else
            Nt += Np
            IH[i1].Dime = Nt
            IH[i1].Dimn = Np

            for idim in 1:Np
                Dim0 = structDimension(IH[i1].iCH, SYS.kv[idim], SYS.wv[idim], 0)
                push!(Dim, Dim0)

            end
        end

    end


    return SYS, IH, Dim
end
#*******************************************************************************************
function Independent_amp(Sys, channels, CC, qn; Np=10, Nx=50, Nphi=100)
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
                if haskey(CC, chname) == true
                    Nex = length(CC[chname])
                    IA0 = structInterAction(
                        Nex,
                        [CC[chname][i][1][1] for i in 1:Nex],  #name
                        [pkey[CC[chname][i][1][1]] for i in 1:Nex],  #key
                        [p[pkey[CC[chname][i][1][1]]].J for i in 1:Nex],  #J
                        [p[pkey[CC[chname][i][1][1]]].Jh for i in 1:Nex],  #Jh
                        [p[pkey[CC[chname][i][1][1]]].P for i in 1:Nex],  #P
                        [p[pkey[CC[chname][i][1][1]]].m for i in 1:Nex],  #m
                        [CC[chname][i][1][2] for i in 1:Nex],  #dc
                        [CC[chname][i][2] for i in 1:Nex]     #CC
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
                        []     #CC
                    )
                end
                IA[i1, i2] = IA0
                IA[i2, i1] = IA0
            end

        end
    end
    kv, wv = gausslaguerre(Np)
    wv = wv .* exp.(kv)
    #kv = [2.0496633800321580e-002, 0.10637754206664184, 0.25725070911685544,
    #    0.47691569315652682, 0.78977094890173449, 1.2661898317497706, 2.0968065118988930,
    #    3.8872580463677919, 9.4004780129091756, 48.788435306583473]
    #wv = [5.2385549033825828e-002, 0.11870709308644416, 0.18345726134358431,
    #    0.25958276865541480, 0.37687641122072230, 0.60422211564747619, 1.1412809934307626,
    #    2.7721814583372204, 10.490025610274076, 124.69392072758504]

    wd = [wignerd(qn.J / qn.Jh, acos(-1.0 + 2.0 / Nx * (i - 0.5))) for i in 1:Nx]

    dphi = 2.0 * pi / Nphi
    sp = [sin((i - 1) * dphi) for i in 1:Nphi]
    cp = [cos((i - 1) * dphi) for i in 1:Nphi]
    SYS = structSys(Sys, kv, wv, wd, sp, cp)

    return SYS, IA, CH, IH
end
end
