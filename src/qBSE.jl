#*******************************************************************************************
#*******************************************************************************************
module qBSE
using FastGaussQuadrature, StaticArrays, WignerD
using Base.Threads
using StaticArrays
using  ..QFT
# store the information of a channel. The kv and wv are stored here for convenience. 
mutable struct ChannelBasisType{VI<:Vector{Int64},I<:Int64,VF<:Vector{Float64},F<:Float64,VS<:Vector{String}} #CB
    p::VI #particles of channel
    cutoff::F #cutoff of channel
    name::VS  #name of channel
    name0::VS #name of channel (without charge)
    Nhel::I #total number of helicity amplitudes of a channel
    kv::VF #momenta of discreting
    wv::VF #weights of discreting
end
# store the information of a interaction
mutable struct InteractionType{I<:Int64,MI<:Matrix{Int64},MF<:Matrix{Float64}} #IA
    Nex::I  #total number of exchanges
    ex::MI  #exchange
    CC::MF  #flavor factors
end
# store the information of an independent helicity for matrix
mutable struct IndependentHelicityType{I<:Int64,V<:Vector{Int64},C<:Complex{Float64}} #IH
    ich::I # which channel this independent helicity belongs to 
    hel::V # helicities
    hel_lh::V # fermion or boson
    Nmo::I # +1 or not for onshell W>sum m. Np+Nmo is the total dimensions of a independent helicities
    Nm0::I #the total number of dimensions before this independent helicities
    kon::C # onshell momentum
end
struct QuantumNumberType{I<:Int64}
    lII::I
    lII_h::I
    lJJ::I
    lJJ_h::I
    lPP::I
    lCC::I
end
mutable struct MomentaType{SV<:SVector{5,ComplexF64},C<:ComplexF64} #k
    i1::SV
    f1::SV
    i2::SV
    f2::SV
    qd::SV
    qd2::C
end


struct RangeType{F<:Float64,I<:Int64}
    ERmin::F
    ERmax::F
    NER::I
    EIt::F
    NEI::I
end
mutable struct PoleType{F<:Float64}
    Ampmin::F
    Ampminx::F
    Ampminy::F
end
mutable struct HelicitiesType{I<:Int64} #l
    i1::I
    i1_h::I
    f1::I
    f1_h::I
    i2::I
    i2_h::I
    f2::I
    f2_h::I
end
const pkey = Dict(
    "test1" => 55, "test2" => 56, "test3" => 57, "test4" => 58, "V" => 70,
    "pi" => 9, "pi_p" => 10, "pi_0" => 11, "pi_m" => 12,
    "rho" => 34, "rho_p" => 35, "rho_0" => 36, "rho_m" => 37,
    "K" => 1, "K_m" => 2, "K_p" => 3, "K_b0" => 4, "K_0" => 5,
    "B" => 78, "B_m" => 79, "B_p" => 80, "B_b0" => 81, "B_0" => 82,
    "D" => 22, "D_m" => 23, "D_p" => 24, "D_b0" => 25, "D_0" => 26,
    "KA" => 59, "KA_m" => 60, "KA_p" => 61, "KA_b0" => 62, "KA_0" => 63,
    "BA" => 83, "BA_m" => 84, "BA_p" => 85, "BA_b0" => 86, "BA_0" => 87,
    "DA" => 27, "DA_m" => 28, "DA_p" => 29, "DA_b0" => 30, "DA_0" => 31,
    "Ds" => 38, "Ds_m" => 39, "Ds_p" => 40,
    "DsA" => 100, "DsA_m" => 101, "DsA_p" => 102,
    "Ds0" => 52, "Ds0_p" => 53, "Ds0_m" => 54,
    "Ds1p" => 97, "Ds1p_m" => 98, "Ds1p_p" => 99,
    "D1p" => 50,
    "eta" => 17, "etap" => 77, "phi" => 51, "omega" => 68, "sig" => 69,
    "JPsi" => 32, "etac" => 33, "Upsilon" => 118,
    "N" => 6, "N_p" => 7, "N_n" => 8,
    "Delta" => 72, "Delta_pp" => 73, "Delta_p" => 74, "Delta_0" => 75, "Delta_m" => 76,
    "Sigma" => 13, "Sigma_p" => 14, "Sigma_0" => 15, "Sigma_m" => 16,
    "Sigmab" => 88, "Sigmab_p" => 89, "Sigmab_0" => 90, "Sigmab_m" => 91,
    "Sigmac" => 41, "Sigmac_0" => 42, "Sigmac_p" => 43, "Sigmac_pp" => 44,
    "SigmaA" => 64, "SigmaA_p" => 65, "SigmaA_0" => 66, "SigmaA_m" => 67,
    "SigmabA" => 92, "SigmabA_p" => 93, "SigmabA_0" => 94, "SigmabA_m" => 95,
    "SigmacA" => 45, "SigmacA_0" => 46, "SigmacA_p" => 47, "SigmacA_pp" => 48,
    "Xi" => 19, "Xi_m" => 20, "Xi_0" => 21,
    "XiA" => 114, "XiA_m" => 115, "XiA_0" => 116,
    "Xic123n" => 103, "Xic123n_p" => 104, "Xic123n_0" => 105,
    "Xic12p" => 106, "Xic12p_p" => 107, "Xic12p_0" => 108,
    "Xic12pA" => 109, "Xic12pA_p" => 110, "Xic12pA_0" => 111,
    "Omegac" => 112, "OmegacA" => 113, "Omegan" => 117,
    "Lambda" => 18, "Lambdab" => 96, "Lambdac" => 49,
    "Lambda_b" => 71
)
struct ParticlesType
    name::String
    name0::String
    m::Float64
    J::Int64
    Jh::Int64
    P::Int64
end
function particles!(particles::Vector{ParticlesType}, filename::String)
    open(filename, "r") do file
        Npar = parse(Int, split(readline(file))[1])
        for i in 1:Npar
            line = split(readline(file))
            name = line[2]
            name0 = line[3]
            m = parse(Float64, line[4])
            J = parse(Int64, line[5])
            Jh = parse(Int64, line[6])
            P = parse(Int64, line[7])
            push!(particles, ParticlesType(name, name0, m, J, Jh, P))
        end
    end
end
const p = ParticlesType[]
#*******************************************************************************************

function fPropFF(k, ex, L, LLi, LLf, lregu, lFFex,qt) #慢

    m = p[ex].m
    if lFFex >= 10
        L = m + 0.22 * L
        LLi = m + 0.22 * LLi
        LLf = m + 0.22 * LLf
        lFFex -= 10
    end

    fProp = 1.0 + 0im
    FFex = 1.0 + 0im
    FFre = 0.0 + 0im

    if ex != 1
        qd2=k.qd2
        fProp *= (1.0 + 0im) / (qd2 - ComplexF64(m^2))
        if lFFex == 1
            FFex *= ((L^2 - m^2) / (L^2 - k.qd2))^2
        elseif lFFex == 2
            FFex *= (L^4 / ((m^2 - k.qd2)^2 + L^4))^2
        elseif lFFex == 3
            FFex *= exp(-(m^2 - k.qd2)^2 / L^4)
        elseif lFFex == 4
            FFex *= ((L^4 + (qt - m^2)^2 / 4) / ((k.qd2 - (qt + m^2) / 2)^2 + L^4))^2
        elseif lFFex == 5
            FFex *= exp(-(m^2 - k.qd2)^2 / L^4)^2
        end
    end

    if lregu == 1
        mi1 = real(k.i1[5])
        mf1 = real(k.f1[5])
        mi2 = real(k.i2[5])
        mf2 = real(k.f2[5])

        if mi1 <= mi2
            FFre += -(mi1^2 - QFT.cdot(k.i1, k.i1))^2 / LLi^4
        end
        if mi1 > mi2
            FFre += -(mi2^2 - QFT.cdot(k.i2, k.i2))^2 / LLi^4
        end
        if mf1 <= mf2
            FFre += -(mf1^2 - QFT.cdot(k.f1, k.f1))^2 / LLf^4
        end
        if mf1 > mf2
            FFre += -(mf2^2 - QFT.cdot(k.f2, k.f2))^2 / LLf^4
        end
    end
    return fProp * FFex * exp(FFre)
end
function fKernel(kf, ki, iNih1, iNih2, Ec, qn, lregu, lFFex, CB, IA, IH, fV)::ComplexF64
    ichi = IH[iNih2].ich
    ichf = IH[iNih1].ich

    if IA[ichi, ichf].Nex == 0
        return 0 + 0im
    end

    mi1, mi2 = p[CB[ichi].p[1]].m, p[CB[ichi].p[2]].m
    mf1, mf2 = p[CB[ichf].p[1]].m, p[CB[ichf].p[2]].m

    Ei1, Ei2 = sqrt(ki^2 + mi1^2), sqrt(ki^2 + mi2^2)
    Ef1, Ef2 = sqrt(kf^2 + mf1^2), sqrt(kf^2 + mf2^2)

    ki10, ki20, qti = (mi1 + 1e-7 <= mi2) ? (Ec - Ei2, Ei2, mi2) : (Ec - Ei1, Ei1, Ec - mi1)
    kf10, kf20, qtf = (mf1 + 1e-7 <= mf2) ? (Ec - Ef2, Ef2, mf2) : (Ec - Ef1, Ef1, Ec - mf1)
    qt = (qti - qtf)^2

    l = HelicitiesType(
        IH[iNih2].hel[1], IH[iNih2].hel_lh[1],
        IH[iNih1].hel[1], IH[iNih1].hel_lh[1],
        IH[iNih2].hel[2], IH[iNih2].hel_lh[2],
        IH[iNih1].hel[2], IH[iNih1].hel_lh[2]
    )

    eta = p[CB[ichi].p[1]].P * p[CB[ichi].p[2]].P * qn.lPP *
          (-1)^(qn.lJJ / qn.lJJ_h - p[CB[ichi].p[1]].J / p[CB[ichi].p[1]].Jh - p[CB[ichi].p[2]].J / p[CB[ichi].p[2]].Jh)

    Ker0 = 0 + 0im
    d50 = 2.0 / 50.0  # 常数预计算
    for i in 1:50
        x = -1.0 + d50 * (i - 0.5)
        sqrt1_x2 = sqrt(1 - x^2)
        ki1 = @SVector [0.0 + 0im, 0.0 + 0im, -ki + 0im, ki10 + 0im, mi1 + 0im]
        kf1 = @SVector [-kf * sqrt1_x2 + 0im, 0.0 + 0im, -kf * x + 0im, kf10 + 0im, mf1 + 0im]
        ki2 = @SVector [0.0 + 0im, 0.0 + 0im, ki + 0im, ki20 + 0im, mi2 + 0im]
        kf2 = @SVector [kf * sqrt1_x2 + 0im, 0.0 + 0im, kf * x + 0im, kf20 + 0im, mf2 + 0im]
        qd = kf2 - ki2
        qd2 = QFT.cdot(qd, qd)
        k = MomentaType(ki1, kf1, ki2, kf2, qd, qd2)

        lJJ = qn.lJJ / qn.lJJ_h
        l21i = -l.i2 / l.i2_h + l.i1 / l.i1_h
        l21f = l.f2 / l.f2_h - l.f1 / l.f1_h
        acx = acos(x)
        DJ = [
            WignerD.wignerdjmn(lJJ, l21i, l21f, acx),
            WignerD.wignerdjmn(lJJ, -l21i, l21f, acx)
        ]

        ker = 0 + 0im
        for ilV in -1:2:1
            dj_idx = div(ilV + 1, 2) + 1
            ker += fV(k, l, ilV, CB[ichi].cutoff, CB[ichf].cutoff, lregu, lFFex, qt) * IA[ichi, ichf].CC[1, 1] * DJ[dj_idx]
            if ilV == -1
                ker *= eta
            end
        end
        Ker0 += ker * d50
    end

    fKernel = Ker0 / 2.0
    if l.f1 == 0 && l.f2 == 0
        fKernel /= sqrt(2.0)
    end
    if l.i1 == 0 && l.i2 == 0
        fKernel /= sqrt(2.0)
    end
    return fKernel
end
#*******************************************************************************************
function fProp(iNih, ii, rp, Ec, Np, CB, IH)
    Er = real(Ec)
    fprop = Complex{Float64}(0, 0)
    ich = IH[iNih].ich
    mi1 = p[CB[ich].p[1]].m
    mi2 = p[CB[ich].p[2]].m
    rL2 = CB[ich].cutoff
    mi2p2 = 1.0
    if IH[iNih].hel_lh[1] == 2
        mi2p2 *= 2.0 * mi1
    end
    if IH[iNih].hel_lh[2] == 2
        mi2p2 *= 2.0 * mi2
    end

    if mi1 <= mi2
        E2 = real(sqrt(rp^2 + mi2^2))
        k12 = ((Ec) - E2)^2 - rp^2
        if ii <= Np
            cE1 = sqrt(rp^2 + mi1^2)
            cE2 = sqrt(rp^2 + mi2^2)
            fprop = mi2p2 * rp^2 * CB[ich].wv[ii] / (2 * pi^2) / (2 * cE2 * ((Ec - cE2)^2 - cE1^2))
        end
    else
        E1 = real(sqrt(rp^2 + mi1^2))
        k22 = ((Ec) - E1)^2 - rp^2
        ff = exp(-(mi2^2 - k22)^2 / rL2^4)^2
        if ii <= Np
            cE2 = sqrt(rp^2 + mi2^2)
            cE1 = sqrt(rp^2 + mi1^2)
            fprop = mi2p2 * rp^2 * CB[ich].wv[ii] / (2 * pi^2) / (2 * cE1 * ((Ec - cE1)^2 - cE2^2))
        end
    end

    if ii == Np + 1
        konc = sqrt((Ec^2 - (mi1 + mi2)^2) * (Ec^2 - (mi2 - mi1)^2)) / (2 * Ec)
        konc2 = (Ec^2 - (mi1 + mi2)^2) * (Ec^2 - (mi2 - mi1)^2) / (4 * Ec^2)
        if imag(konc) > 0
            fprop += Complex{Float64}(0, 1) * mi2p2 * konc / (8 * pi * Ec)
        else
            fprop -= Complex{Float64}(0, 1) * mi2p2 * konc / (8 * pi * Ec)
        end
        for i3 in 1:Np
            fprop += mi2p2 * CB[ich].wv[i3] * konc2 / (2 * pi^2) / (2 * Ec * (CB[ich].kv[i3]^2 - konc2))
        end
    end
    return fprop
end
function srAB(Ec, qn, lregu, lFFex, lRm, Np, Nih, CB, IH, IA,fV)
    # Determine the dimension of the work matrix. At energies larger than the threshold, the dimension increases by 1
    Nt, IH = WORKSPACE(Ec, lRm, Np, Nih, CB, IH)
    Vc = zeros(Complex{Float64}, Nt, Nt)
    Gc = zeros(Complex{Float64}, Nt, Nt)
    II = zeros(Float64, Nt, Nt)
    for iNih1 in 1:Nih
        for iNih2 in 1:Nih
            for i1 in 1:(Np+IH[iNih1].Nmo)
                kf = (i1 <= Np) ? CB[IH[iNih1].ich].kv[i1] : IH[iNih1].kon
                for i2 in 1:(Np+IH[iNih2].Nmo)
                    ki = (i2 <= Np) ? CB[IH[iNih2].ich].kv[i2] : IH[iNih2].kon
                    if i1 + IH[iNih1].Nm0 <= i2 + IH[iNih2].Nm0
                        Vc[i1+IH[iNih1].Nm0, i2+IH[iNih2].Nm0] = fKernel(kf, ki, iNih1, iNih2, Ec, qn, lregu, lFFex, CB, IA, IH,fV)::ComplexF64
                    end
                    if i1 == i2 && iNih1 == iNih2
                        Gc[i1+IH[iNih1].Nm0, i2+IH[iNih2].Nm0] = fProp(iNih1, i1, kf, Ec, Np, CB, IH)
                        II[i1+IH[iNih1].Nm0, i2+IH[iNih2].Nm0] = 1.0
                    end
                end
            end
        end
    end

    for iNih1 in 1:Nih
        for iNih2 in 1:Nih
            for i1 in 1:(Np+IH[iNih1].Nmo)
                for i2 in 1:(Np+IH[iNih2].Nmo)
                    if i1 + IH[iNih1].Nm0 > i2 + IH[iNih2].Nm0
                        if IH[iNih1].ich != IH[iNih2].ich
                            Vc[i1+IH[iNih1].Nm0, i2+IH[iNih2].Nm0] = conj(Vc[i2+IH[iNih2].Nm0, i1+IH[iNih1].Nm0])
                        else
                            Vc[i1+IH[iNih1].Nm0, i2+IH[iNih2].Nm0] = Vc[i2+IH[iNih2].Nm0, i1+IH[iNih1].Nm0]
                        end
                    end
                end
            end
        end
    end

    return Vc, Gc, II
end
#*******************************************************************************************
function WORKSPACE(Ec, lRm, Np, Nih, CB, IH)
    Ec += Complex{Float64}(0.0, 1e-15)
    Nt = 0
    for i1 in 1:Nih
        IH[i1].Nmo = 0
        IH[i1].Nm0 = 0
    end
    for i1 in 1:Nih
        IH[i1].Nm0 = Nt
        mass1 = p[CB[IH[i1].ich].p[1]].m
        mass2 = p[CB[IH[i1].ich].p[2]].m
        Rm = false
        if lRm == 1
            Rm = real(Ec) > (mass1 + mass2)
        elseif lRm == 2
            Rm = true
        end
        if Rm
            IH[i1].kon = sqrt((Ec^2 - (mass1 + mass2)^2) * (Ec^2 - (mass1 - mass2)^2)) / (2.0 * Ec)
            IH[i1].Nmo = 1
            Nt += Np + IH[i1].Nmo
        else
            Nt += Np
        end

    end
    return Nt, IH
end
#*******************************************************************************************
function Independent_amp(channels, CC, lJJ, Np)
    CB = ChannelBasisType[]
    IH = IndependentHelicityType[]
    Nih, Nc = 1, 1
    for channel in channels
        p1 = pkey[channel[1]]
        p2 = pkey[channel[2]]
        Nhel = 0
        CB0 = ChannelBasisType(
            [p1, p2],
            channel[3],
            [p[p1].name, p[p2].name, "$(p[p1].name):$(p[p2].name)"],
            [p[p1].name0, p[p2].name0, "$(p[p1].name0):$(p[p2].name0)"],
            0,
            Float64[],
            Float64[]
        )



        CB0.kv, CB0.wv = gausslaguerre(Np)
        CB0.wv = CB0.wv .* exp.(CB0.kv)
        #CB0.kv = [2.0496633800321580e-002, 0.10637754206664184, 0.25725070911685544, 0.47691569315652682, 0.78977094890173449, 1.2661898317497706, 2.0968065118988930, 3.8872580463677919, 9.4004780129091756, 48.788435306583473]
        #CB0.wv = [5.2385549033825828e-002, 0.11870709308644416, 0.18345726134358431, 0.25958276865541480, 0.37687641122072230, 0.60422211564747619, 1.1412809934307626, 2.7721814583372204, 10.490025610274076, 124.69392072758504]
        Nih0 = Nih
        IH0 = IndependentHelicityType[]
        IH00 = IndependentHelicityType(0, Int64[], Int64[], 0, 0, 0.0 + 0.0 * im)
        for i1 in -p[p1].J:p[p1].Jh:p[p1].J
            for i2 in -p[p2].J:p[p2].Jh:p[p2].J
                Jh = 1.0
                if p[p1].Jh == 1 && p[p2].Jh == 2 || p[p1].Jh == 2 && p[p2].Jh == 1
                    Jh = 2.0
                end
                if abs(float(i1) / float(p[p1].Jh) - float(i2) / float(p[p2].Jh)) <= float(lJJ) + 0.01
                    lind = 1
                    for i3 in 1:Nih-Nih0 #Nih0:(Nih-1)
                        if IH0[i3].hel[1] == -i1 && IH0[i3].hel[2] == -i2
                            lind = 0
                            break
                        end
                    end
                    if lind == 1
                        IH00.hel = [i1, i2]
                        IH00.hel_lh = [p[p1].Jh, p[p2].Jh]
                        IH00.ich = Nc
                        push!(IH0, IH00)
                        Nih += 1
                        Nhel += 1
                    end
                end
            end
        end
        CB0.Nhel = Nhel
        Nc += 1

        push!(CB, CB0)
        append!(IH, IH0)
    end

    Nc -= 1
    Nih -= 1


    chnamei = Vector{String}(undef, 1)
    chnamef = Vector{String}(undef, 1)
    IA = Matrix{InteractionType}(undef, Nc, Nc)
    for i1 in 1:Nc
        for i2 in 1:Nc
            chnamei[1] = CB[i1].name[3]
            chnamef[1] = CB[i2].name[3]
            chname = "$(chnamei[1])-->$(chnamef[1])"

            IA0 = InteractionType(0, zeros(Int64, 5, 2), zeros(Float64, 5, 2))

            IA0.Nex = CC[chname][1]
            IA0.ex = Matrix{Int64}(undef, IA0.Nex, 2)
            IA0.ex[1, :] = CC[chname][2]
            IA0.CC = Matrix{Int64}(undef, IA0.Nex, 2)
            IA0.CC[1, 1] = CC[chname][3]

            IA[i1, i2] = IA0
        end
    end

    return Nih, Nc, CB, IH, IA
end
end
