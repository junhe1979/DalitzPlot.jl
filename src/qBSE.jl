#*******************************************************************************************
#*******************************************************************************************
module qBSE
using FastGaussQuadrature, StaticArrays, WignerD
using Base.Threads
using StaticArrays

# store the information of a channel. The kv and wv are stored here for convenience. 
struct ChannelBasisType #CB
    p::Vector{Int64} #particles of channel
    m::Vector{Float64} #particles of channel
    J::Vector{Int64} #particles of channel
    Jh::Vector{Int64} #particles of channel
    P::Vector{Int64} #particles of channel
    cutoff::Float64 #cutoff of channel
    name::Vector{String}  #name of channel
    name0::Vector{String} #name of channel (without charge)
    Nhel::Int64 #total number of helicity amplitudes of a channel
    kv::Vector{Float64} #momenta of discreting
    wv::Vector{Float64} #weights of discreting
end
# store the information of a interaction
struct InterActionType #IA
    Project::String
    Nex::Int64  #total number of exchanges
    name_ex::Vector{String}  #exchange
    key_ex::Vector{Int64}
    J_ex::Vector{Int64}
    Jh_ex::Vector{Int64}
    P_ex::Vector{Int64}
    m_ex::Vector{Float64}  #exchange
    dc::Vector{Int64}  #direct or cross
    CC::Vector{Float64}  #flavor factors
    D::Vector{Matrix{Float64}}   #WignerD
end
# store the information of an independent helicity for matrix
mutable struct IndependentHelicityType #IH
    ich::Int64 # which channel this independent helicity belongs to 
    hel::Vector{Int64} # helicities
    helh::Vector{Int64} # fermion or boson
    Nmo::Int64 # +1 or not for onshell W>sum m. Np+Nmo is the total dimensions of a independent helicities
    Nm0::Int64 #the total number of dimensions before this independent helicities
    kon::Complex{Float64} # onshell momentum
end
mutable struct MomentaType #momenta
    i1::SVector{5,ComplexF64} # initial 1
    f1::SVector{5,ComplexF64}
    i2::SVector{5,ComplexF64}
    f2::SVector{5,ComplexF64}
    q::SVector{5,ComplexF64} #exchange
    q2::Complex{Float64} # q^2
    qt::Complex{Float64} #
end
mutable struct HelicitiesType #l
    i1::Int64
    i1h::Int64
    f1::Int64
    f1h::Int64
    i2::Int64
    i2h::Int64
    f2::Int64
    f2h::Int64
end
struct ParticlesType #store the information of particles
    name::String
    name0::String
    m::Float64
    J::Int64
    Jh::Int64
    P::Int64
end
function particles!(particles::Vector{ParticlesType}, pkey::Dict{String,Int64}, filename::String) #read the information 
    open(filename, "r") do file
        readline(file)
        i = 1
        for line in eachline(file)
            parts = split(line)
            particle = ParticlesType(parts[1], parts[2], parse(Float64, parts[3]),
                parse(Int64, parts[4]), parse(Int64, parts[5]), parse(Int64, parts[6]))
            push!(particles, particle)
            pkey[parts[1]] = i
            i += 1
        end
    end
end
const p = ParticlesType[] #store of information of particles in this global vector
const pkey = Dict{String,Int64}() #store the dictionary for key and the number of particles
#*******************************************************************************************
# Calculate the \sum|M|^2 for each channel 
function M2_channel(T,IH,CB)
  
    Nih = length(IH)
    NCB = length(CB)
    resM2 = zeros(Float64, NCB, NCB)
    Np = length(CB[1].kv)
    resA = [IH[iNih1].Nmo == 1 && IH[iNih2].Nmo == 1 ?
            T[IH[iNih1].Nm0+Np+1, IH[iNih2].Nm0+Np+1] : 0.0
            for iNih1 in 1:Nih, iNih2 in 1:Nih]
    #@show ER, abs2(resA[2,2])
    hel1 = 1
    for iCB1 in 1:NCB
        hel2 = 1
        for iCB2 in 1:NCB
            for i in hel1:hel1+CB[iCB1].Nhel-1, j in hel2:hel2+CB[iCB2].Nhel-1
                #@show NCB, iCB1, iCB2, i, j, size(resA)
                resM2[iCB1, iCB2] += abs2(resA[i, j])
            end
            hel2 += CB[iCB2].Nhel
        end
        hel1 += CB[iCB1].Nhel
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

function simpleXsection(xx, M2, CB; Ep="cm")

    Xsection = Matrix{Float64}[]
    NCB = length(CB)
    N=length(xx)
    for i in 1:N
        Xs0 = zeros(size(M2[1]))
        for iM2 in 1:NCB
            for fM2 in 1:NCB
                mi1 = p[CB[iM2].p[1]].m
                mi2 = p[CB[iM2].p[2]].m
                mf1 = p[CB[fM2].p[1]].m
                mf2 = p[CB[fM2].p[2]].m
                W = xx[i]
                if Ep == "L"
                    W = sqrt((sqrt(W^2 + mi1^2) + mi2)^2 - W^2)
                end
                fac = lambda(W, mf1, mf2) / lambda(W, mi1, mi2) / W^2 * (2.0 * mi2 * 2.0 * mi2 * 0.3894 / 16.0 / pi)
                Xs0[iM2, fM2] = M2[i][iM2, fM2] * fac
            end
        end
        push!(Xsection,Xs0)
    end
    return Xsection
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
function fKernel(kf, ki, iNih1, iNih2, Ec, qn, CB, IA, IH, fV)::ComplexF64 # Calculating fKernel
    ichi = IH[iNih2].ich # channel
    ichf = IH[iNih1].ich

    if IA[ichi, ichf].Nex == 0 # no exchange, returen 0
        return 0 + 0im
    end

    mi1, mi2 = CB[ichi].m[1], CB[ichi].m[2]
    mf1, mf2 = CB[ichf].m[1], CB[ichf].m[2]

    Ei1, Ei2 = sqrt(ki^2 + mi1^2), sqrt(ki^2 + mi2^2)
    Ef1, Ef2 = sqrt(kf^2 + mf1^2), sqrt(kf^2 + mf2^2)

    ki10, ki20, qti = (mi1 + 1e-7 <= mi2) ? (Ec - Ei2, Ei2, mi2) : (Ei1, Ec - Ei1, Ec - mi1)
    kf10, kf20, qtf = (mf1 + 1e-7 <= mf2) ? (Ec - Ef2, Ef2, mf2) : (Ef1, Ec - Ef1, Ec - mf1)
    # @show ki10, ki20, qti 
    qt = (qti - qtf)^2 + 0im

    l = HelicitiesType(  #helicities
        IH[iNih2].hel[1], IH[iNih2].helh[1],
        IH[iNih1].hel[1], IH[iNih1].helh[1],
        IH[iNih2].hel[2], IH[iNih2].helh[2],
        IH[iNih1].hel[2], IH[iNih1].helh[2]
    )

    eta = p[CB[ichi].p[1]].P * CB[ichi].P[2] * qn.P *
          (-1)^(qn.J / qn.Jh - CB[ichi].J[1] / CB[ichi].Jh[1] - CB[ichi].J[2] / CB[ichi].Jh[2])

    lJJ = qn.J / qn.Jh
    l21i = -l.i2 / l.i2h + l.i1 / l.i1h
    l21f = l.f2 / l.f2h - l.f1 / l.f1h
    lf = Int64(lJJ + l21f) + 1
    d50 = 2.0 / 50.0
    Ker0 = 0 + 0im
    for i in 1:50
        x = -1.0 + d50 * (i - 0.5)
        sqrt1_x2 = sqrt(1 - x^2)
        ki1 = @SVector [0.0 + 0im, 0.0 + 0im, -ki + 0im, ki10 + 0im, mi1 + 0im]
        kf1 = @SVector [-kf * sqrt1_x2 + 0im, 0.0 + 0im, -kf * x + 0im, kf10 + 0im, mf1 + 0im]
        ki2 = @SVector [0.0 + 0im, 0.0 + 0im, ki + 0im, ki20 + 0im, mi2 + 0im]
        kf2 = @SVector [kf * sqrt1_x2 + 0im, 0.0 + 0im, kf * x + 0im, kf20 + 0im, mf2 + 0im]
        q = kf2 - ki2
        q2 = q * q
        k = MomentaType(ki1, kf1, ki2, kf2, q, q2, qt)

        l.f1, l.i2 = -l.f1, -l.i2  #helicity to spin and the minus for fixed parity
        Kerm1 = fV(k, l, CB, IA, ichi, ichf) * IA[ichi, ichf].D[i][Int64(lJJ + l21i)+1, lf] * eta
        l.f1, l.i2 = -l.f1, -l.i2
        l.i1, l.f1 = -l.i1, -l.f1
        Kerp1 = fV(k, l, CB, IA, ichi, ichf) * IA[ichi, ichf].D[i][Int64(lJJ - l21i)+1, lf]
        l.i1, l.f1 = -l.i1, -l.f1
        Ker0 += (Kerp1 + Kerm1) * d50
    end

    fKernel = Ker0 / 2.0 ## 2pi/4pi 
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
function fProp(iNih, ii, rp, Ec, Np, CB, IH)
    fprop = Complex{Float64}(0, 0)
    ich = IH[iNih].ich
    mi1 = CB[ich].m[1]
    mi2 = CB[ich].m[2]
    mi2p2 = 1.0
    if IH[iNih].helh[1] == 2
        mi2p2 *= 2.0 * mi1
    end
    if IH[iNih].helh[2] == 2
        mi2p2 *= 2.0 * mi2
    end
    # for off momenta
    if mi1 <= mi2
        if ii <= Np
            cE1 = sqrt(rp^2 + mi1^2)
            cE2 = sqrt(rp^2 + mi2^2)
            fprop = mi2p2 * rp^2 * CB[ich].wv[ii] / (2 * pi^2) / (2 * cE2 * ((Ec - cE2)^2 - cE1^2))
        end
    else
        if ii <= Np
            cE2 = sqrt(rp^2 + mi2^2)
            cE1 = sqrt(rp^2 + mi1^2)
            fprop = mi2p2 * rp^2 * CB[ich].wv[ii] / (2 * pi^2) / (2 * cE1 * ((Ec - cE1)^2 - cE2^2))
        end
    end
    # for onshell
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
function srAB(Ec, qn, CB, IH, IA, fV; lRm=1)
    # Determine the dimension of the work matrix. At energies larger than the threshold, the dimension increases by 1
    Np = length(CB[1].kv)
    Nih = length(IH)
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
                        Vc[i1+IH[iNih1].Nm0, i2+IH[iNih2].Nm0] = fKernel(kf, ki, iNih1, iNih2, Ec, qn, CB, IA, IH, fV)::ComplexF64
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

    return Vc, Gc, II, IH
end
#*******************************************************************************************
function WORKSPACE(Ec, lRm, Np, Nih, CB, IH)
    #produce the dimensions in a independent helicity amplitudes, especially for the onshell dimension.
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
        if lRm == 1 #for first Reimann sheet
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
function Independent_amp(Project, channels, CC, qn; Np=10, Nx=50)
    #store the channel information, interaction information in CB,IH,IA
    CB = ChannelBasisType[]
    IH = IndependentHelicityType[]
    Nih, Nc = 1, 1
    #kv, wv = gausslaguerre(Np)
    #wv = wv .* exp.(kv)
    kv = [2.0496633800321580e-002, 0.10637754206664184, 0.25725070911685544,
        0.47691569315652682, 0.78977094890173449, 1.2661898317497706, 2.0968065118988930,
        3.8872580463677919, 9.4004780129091756, 48.788435306583473]
    wv = [5.2385549033825828e-002, 0.11870709308644416, 0.18345726134358431,
        0.25958276865541480, 0.37687641122072230, 0.60422211564747619, 1.1412809934307626,
        2.7721814583372204, 10.490025610274076, 124.69392072758504]

    wD = [wignerd(qn.J / qn.Jh, acos(-1.0 + 2.0 / Nx * (i - 0.5))) for i in 1:Nx]
    for channel in channels
        p1, p2 = pkey[channel[1]], pkey[channel[2]]
        m1, m2 = p[p1].m, p[p2].m
        J1, J2 = p[p1].J, p[p2].J
        Jh1, Jh2 = p[p1].Jh, p[p2].Jh
        P1, P2 = p[p1].P, p[p2].P
        Nhel = 0
        #independent helicities
        Nih0 = Nih
        IH0 = IndependentHelicityType[]
        IH00 = IndependentHelicityType(0, Int64[], Int64[], 0, 0, 0.0 + 0.0 * im)
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
                        IH00.ich = Nc
                        push!(IH0, deepcopy(IH00))

                        Nih += 1
                        Nhel += 1
                    end
                end
            end
        end
        Nc += 1
        #channel information
        CB0 = ChannelBasisType(
            [p1, p2],
            [m1, m2],
            [J1, J2],
            [Jh1, Jh2],
            [P1, P2],
            channel[3],
            [p[p1].name, p[p2].name, "$(p[p1].name):$(p[p2].name)"],
            [p[p1].name0, p[p2].name0, "$(p[p1].name0):$(p[p2].name0)"],
            Nhel,
            kv,
            wv
        )
        push!(CB, deepcopy(CB0))
        append!(IH, deepcopy(IH0))
    end

    Nc -= 1
    Nih -= 1

    # interaction information.
    chnamei = Vector{String}(undef, 1)
    chnamef = Vector{String}(undef, 1)
    IA = Matrix{InterActionType}(undef, Nc, Nc)
    for i1 in 1:Nc
        for i2 in 1:Nc
            chnamei[1] = CB[i1].name[3]
            chnamef[1] = CB[i2].name[3]
            chname = "$(chnamei[1])-->$(chnamef[1])"

            if i1 <= i2
                if haskey(CC, chname) == true
                    Nex = length(CC[chname])
                    IA0 = InterActionType(
                        Project,
                        Nex,
                        [CC[chname][i][1][1] for i in 1:Nex],  #name
                        [pkey[CC[chname][i][1][1]] for i in 1:Nex],  #key
                        [p[pkey[CC[chname][i][1][1]]].J for i in 1:Nex],  #J
                        [p[pkey[CC[chname][i][1][1]]].Jh for i in 1:Nex],  #Jh
                        [p[pkey[CC[chname][i][1][1]]].P for i in 1:Nex],  #P
                        [p[pkey[CC[chname][i][1][1]]].m for i in 1:Nex],  #m
                        [CC[chname][i][1][2] for i in 1:Nex],  #dc
                        [CC[chname][i][2] for i in 1:Nex],     #CC
                        wD  # Vector{Matrix{Float64}} for D
                    )
                else
                    Nex = 0
                    IA0 = InterActionType(
                        Project,
                        Nex,
                        [],  #name
                        [],  #key
                        [],  #J
                        [],  #Jh
                        [],  #P
                        [],  #m
                        [],  #dc
                        [],     #CC
                        wD  # Vector{Matrix{Float64}} for D
                    )
                end
                IA[i1, i2] = IA0
                IA[i2, i1] = IA0
            end

        end
    end

    return CB, IH, IA
end
end
