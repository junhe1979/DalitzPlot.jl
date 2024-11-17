#*******************************************************************************************
#*******************************************************************************************
module qBSE
using FastGaussQuadrature, StaticArrays, WignerD
using Base.Threads
using StaticArrays

# store the information of a interaction
struct structInterAction #IA
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

# store the information of each dimension for matrix
mutable struct structDimension #Dim
    iIH::Int64 # which independent helicity this dimension belongs to 
    kv::ComplexF64 #momenta of discreting
    wv::Float64 #weights of discreting
    Dimo::Int64 # +1 or not for onshell W>sum m. Np+Dimo is the total dimensions of a independent helicities
end

# store the information of an independent helicity for matrix
mutable struct structIndependentHelicity #IH
    iCH::Int64 # which channel this independent helicity belongs to 
    hel::Vector{Int64} # helicities
    helh::Vector{Int64} # fermion or boson
    Dimi::Int64 #the number of first dimension in this independent helicities
    Dimf::Int64 #the number of last dimension in this independent helicities
    Dimn::Int64 #the number of dimensions in this independent helicities
    Dimo::Int64 # +1 or not for onshell W>sum m. Np+Dimo is the total dimensions of a independent helicities
end

# store the information of a channel. The kv and wv are stored here for convenience. 
struct structChannel #CH
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
function M2_channel(T, IH, CH)
    Nih = length(IH)
    NCH = length(CH)
    Np = length(CH[1].kv)
    resM2 = zeros(Float64, NCH, NCH)

    # 使用数组保存 resA 以减少嵌套循环中的条件判断
    resA = [IH[iNih1].Dimo == 1 && IH[iNih2].Dimo == 1 ?
            T[IH[iNih1].Dimi+Np, IH[iNih2].Dimi+Np] : 0.0
            for iNih1 in 1:Nih, iNih2 in 1:Nih]

    hel1 = 1
    for iCH1 in 1:NCH
        hel2 = 1
        nhel1 = CH[iCH1].Nhel
        for iCH2 in 1:NCH
            nhel2 = CH[iCH2].Nhel
            @inbounds for i in hel1:hel1+nhel1-1, j in hel2:hel2+nhel2-1
                resM2[iCH1, iCH2] += abs2(resA[i, j])
            end
            hel2 += nhel2
        end
        hel1 += nhel1
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
function fKernel(kf, ki, IH1, IH2, Ec, qn, CH, IA, IH, fV)::ComplexF64 # Calculating fKernel
    ichi = IH2.iCH # channel
    ichf = IH1.iCH

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
        IH2.hel[1], IH2.helh[1],
        IH1.hel[1], IH1.helh[1],
        IH2.hel[2], IH2.helh[2],
        IH1.hel[2], IH1.helh[2]
    )

    eta = p[CH[ichi].p[1]].P * CH[ichi].P[2] * qn.P *
          (-1)^(qn.J / qn.Jh - CH[ichi].J[1] / CH[ichi].Jh[1] - CH[ichi].J[2] / CH[ichi].Jh[2])

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
        k = structMomentum(ki1, kf1, ki2, kf2, q, q2, qt)

        l.f1, l.i2 = -l.f1, -l.i2  #helicity to spin and the minus for fixed parity
        Kerm1 = fV(k, l, CH, IA, ichi, ichf) * IA[ichi, ichf].D[i][Int64(lJJ + l21i)+1, lf] * eta
        l.f1, l.i2 = -l.f1, -l.i2
        l.i1, l.f1 = -l.i1, -l.f1
        Kerp1 = fV(k, l, CH, IA, ichi, ichf) * IA[ichi, ichf].D[i][Int64(lJJ - l21i)+1, lf]
        l.i1, l.f1 = -l.i1, -l.f1
        Ker0 += (Kerp1 + Kerm1) * d50
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
function fProp(IH, Dimo, rp, wv, Ec, Np, CH)
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
        cE1, cE2 = sqrt(rp^2 + mi1^2), sqrt(rp^2 + mi2^2)
        denominator = (mi1 <= mi2) ? 2 * cE2 * ((Ec - cE2)^2 - cE1^2) : 2 * cE1 * ((Ec - cE1)^2 - cE2^2)

        if denominator != 0  # To avoid division by zero
            fprop = mi2p2 * rp^2 * wv / pi_factor / denominator
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
            kv_i3 = CH[ich].kv[i3]
            denominator = 2 * Ec * (kv_i3^2 - konc2)
            fprop += mi2p2 * CH[ich].wv[i3] * konc2 / pi_factor / denominator
        end
    end

    return fprop
end
function srAB(Ec, qn, CH, IH, IA, fV; lRm=1)
    # Determine the dimension of the work matrix. 
    Np = length(CH[1].kv)
    Nih = length(IH)
    IH, Dim = WORKSPACE(Ec, lRm, Np, Nih, CH, IH)
    Nt = length(Dim)

    # Initialize matrices
    Vc = zeros(Complex{Float64}, Nt, Nt)
    Gc = zeros(Complex{Float64}, Nt, Nt)
    II = zeros(Float64, Nt, Nt)

    # Populate Vc, Gc, and II
    for i1 = 1:Nt
        Dim1 = Dim[i1]
        kf1 = Dim1.kv

        for i2 = i1:Nt
            Dim2 = Dim[i2]
            ki2 = Dim2.kv

            # Compute Vc for upper triangular part
            Vc[i1, i2] = fKernel(kf1, ki2, IH[Dim1.iIH], IH[Dim2.iIH], Ec, qn, CH, IA, IH, fV)::ComplexF64

            # Compute Gc and II only on the diagonal
            if i1 == i2
                Gc[i1, i2] = fProp(IH[Dim1.iIH], Dim1.Dimo, kf1, Dim1.wv, Ec, Np, CH)
                II[i1, i2] = 1.0
            end
        end
    end

    # Populate lower triangular part of Vc based on conjugate symmetry
    for i1 = 2:Nt
        Dim1_IH = IH[Dim[i1].iIH]

        for i2 = 1:i1-1
            Dim2_IH = IH[Dim[i2].iIH]

            if Dim1_IH.iCH != Dim2_IH.iCH
                Vc[i1, i2] = conj(Vc[i2, i1])
            else
                Vc[i1, i2] = Vc[i2, i1]
            end
        end
    end

    return Vc, Gc, II, IH
end
#*******************************************************************************************
function WORKSPACE(Ec, lRm, Np, Nih, CH, IH)
    #produce the dimensions in a independent helicity amplitudes, especially for the onshell dimension.
    Ec += Complex{Float64}(0.0, 1e-15)
    Dim = structDimension[]
    Nt = 0
    for i1 in 1:Nih
        IH[i1].Dimo = 0
        IH[i1].Dimi = 0
        IH[i1].Dimf = 0
        IH[i1].Dimn = 0
    end
    for i1 in 1:Nih
        IH[i1].Dimi = Nt+1
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
            IH[i1].Dimf = Nt
            IH[i1].Dimn = Np + IH[i1].Dimo
            for idim in 1:Np
                Dim0 = structDimension(IH[i1].iCH, CH[IH[i1].iCH].kv[idim], CH[IH[i1].iCH].wv[idim], 0)
                push!(Dim, Dim0)
            end
            Dim0 = structDimension(IH[i1].iCH, kon, 0.0, 1)
            push!(Dim, Dim0)
        else
            Nt += Np
            IH[i1].Dimf = Nt
            IH[i1].Dimn = Np
            for idim in 1:Np
                Dim0 = structDimension(IH[i1].iCH, CH[IH[i1].iCH].kv[idim], CH[IH[i1].iCH].wv[idim], 0)
                push!(Dim, Dim0)
            end
        end

    end

    return IH, Dim
end
#*******************************************************************************************
function Independent_amp(Project, channels, CC, qn; Np=10, Nx=50)
    #store the channel information, interaction information in CH,IH,IA
    CH = structChannel[]
    IH = structIndependentHelicity[]
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
        IH0 = structIndependentHelicity[]
        IH00 = structIndependentHelicity(0, Int64[], Int64[], 0,0,0,0)
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
                        Nhel += 1
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
            Nhel,
            kv,
            wv
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
                    IA0 = structInterAction(
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

    return CH, IH, IA
end
end
