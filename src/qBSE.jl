# module for quasipotential Bethe-Salpeter equation
module qBSE
using Distributed, FastGaussQuadrature, StaticArrays, WignerD, LinearAlgebra, Printf, ProgressBars, SharedArrays, ProgressMeter
using ..Xs, ..FR
#*******************************************************************************************
# Structures used
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
    channel::Dict{String,Int64}
end
# store the information of a interaction
struct structInterAction #IA[]
    Nex::Int64  #total number of exchanges
    key_ex::Vector{String}  #exchange
    J_ex::Vector{Int64}  #J for integral spin, 2J for half  integral spin
    Jh_ex::Vector{Int64} #1 for integral spin, 2 for half  integral spin
    P_ex::Vector{Int64} #parity
    m_ex::Vector{Float64} #mass
    dc::Vector{Int64}  #direct or cross
    Ff::Vector{Float64}  #flavor factors
end
# store the information of a channel.
struct structChannel #CH[] 
    p::Tuple{String,String}  #partilces of channel
    p_name0::Tuple{String,String} #partilces of channel (without charge)
    anti::Tuple{Int64,Int64} #label for antiparticle 
    m::Tuple{Float64,Float64} #particle masses of channel
    J::Tuple{Int64,Int64} #particle spins of channel
    Jh::Tuple{Int64,Int64} #particle  spins of channel
    P::Tuple{Int64,Int64} #particle parities of channel
    cutoff::Float64 #cutoff of channel
    IHb::Int64 # rank of first independent helicity
    IHe::Int64 # rank of last independent helicity
    IHn::Int64 #total number of helicity amplitudes of a channel
end
# store the information of an independent helicity for matrix
mutable struct structIndependentHelicity #IH[]
    iCH::Int64 # which channel this independent helicity belongs to 
    hel::Tuple{Int64,Int64} # helicities
    helh::Tuple{Int64,Int64} # fermion or boson
    Dimb::Int64 #the rank of begin dimension in this independent helicities
    Dime::Int64 #the rank of end dimension in this independent helicities
    k::ComplexF64 #momenta of discreting
    w::Float64 #weights of discreting
end
#------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------
#store the information of particles involved in the work. 
struct structParticle
    name0::String #without charge, used for the cases with symmetry.
    nameL::String #Latex.
    anti::Int #antiparticle
    m::Float64 #mass
    J::Int64 #spin
    Jh::Int64  #for fermion, the helicity is obtianed by J/Jh
    P::Int64 #parity
end
#* function to read the information of partilces from a file with "filename".
function particles!(particles::Dict{String,structParticle}, filename::String) #read the information 
    open(filename, "r") do file
        readline(file)
        i = 1
        for line in eachline(file)
            parts = split(line)
            particle = structParticle(parts[2], parts[3], parse(Int, parts[4]),
                parse(Float64, parts[5]), parse(Int64, parts[6]), parse(Int64, parts[7]),
                parse(Int64, parts[8]))
            particles[parts[1]] = particle
            i += 1
        end
    end
end
const p = Dict{String,structParticle}() #store of information of particles in this global vector
#*******************************************************************************************
# qBSE
#*******************************************************************************************
#  Prepare
# prepare workspace. -> res  (used in res).  The output will flow to -> res0 -> VGI ->workspace ...
function preprocessing(Sys, qn, channels, Ff, cutoff, Np, Nx, Nphi)
    #store the channel information, interaction information in CH,IH,IA
    CH = structChannel[]
    IH = structIndependentHelicity[]
    Nih, Nc = 1, 1
    channel_dict = Dict{String,Int64}()
    for channel in channels
        p1, p2 = split(channel[1], ":")
        channel_dict[channel[1]] = Nc
        anti1, anti2 = p[p1].anti, p[p2].anti
        m1, m2 = p[p1].m, p[p2].m
        J1, J2 = p[p1].J, p[p2].J
        Jh1, Jh2 = p[p1].Jh, p[p2].Jh
        P1, P2 = p[p1].P, p[p2].P
        IHn = 0
        #independent helicities
        Nih0 = Nih
        IH0 = structIndependentHelicity[]
        IH00 = structIndependentHelicity(0, (0, 0), (0, 0), 0, 0, 0., 0.)
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
                        IH00.hel = (i1, i2)

                        IH00.helh = (p[p1].Jh, p[p2].Jh)
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
            (p1, p2),
            (p[p1].name0, p[p2].name0),
            (anti1, anti2),
            (m1, m2),
            (J1, J2),
            (Jh1, Jh2),
            (P1, P2),
            channel[2],
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
    default_empty = structInterAction(0, [], [], [], [], [], [], [])
    IA = fill(default_empty, Nc, Nc)
    for i1 in 1:Nc
        for i2 in 1:Nc
            chname = "$(CH[i1].p[1]):$(CH[i1].p[2])-->$(CH[i2].p[1]):$(CH[i2].p[2])"

            #chname = ((CH[i1].p[1], CH[i1].p[2]), (CH[i2].p[1], CH[i2].p[2]))
            if haskey(Ff, chname) == true
                Nex = length(Ff[chname])
                IA0 = structInterAction(
                    Nex,
                    [Ff[chname][i][1][1] for i in 1:Nex],  #key
                    [p[Ff[chname][i][1][1]].J for i in 1:Nex],  #J
                    [p[Ff[chname][i][1][1]].Jh for i in 1:Nex],  #Jh
                    [p[Ff[chname][i][1][1]].P for i in 1:Nex],  #P
                    [p[Ff[chname][i][1][1]].m for i in 1:Nex],  #m
                    [Ff[chname][i][1][2] for i in 1:Nex],  #dc
                    [Ff[chname][i][2] for i in 1:Nex]     #Ff
                )

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
        get(cutoff, :FF_ex_type, 3),
        channel_dict)

    return SYS, IA, CH, IH
end
# generate the workspace. -> VGI, 
function workSpace(Ec, lRm, Np, SYS, CH, IH)

    Ec += Complex{Float64}(0.0, 1e-15)

    Nih = length(IH)

    # 先计算总维度 Nt
    Nt = 0
    for i1 in 1:Nih
        mass1 = p[CH[IH[i1].iCH].p[1]].m
        mass2 = p[CH[IH[i1].iCH].p[2]].m

        lRm0 = isa(lRm, Int) ? lRm : lRm[IH[i1].iCH]
        temp = real(Ec) < (mass1 + mass2)
        Rm = !((lRm0 == 0 || lRm0 == 1) && temp)

        Nt += Rm ? (Np + 1) : Np
    end

    # 预分配 3×Nt 矩阵
    Dim = zeros(Int, 3, Nt)

    # 第二遍真正填充
    Nt_count = 0

    for i1 in 1:Nih

        IH[i1].Dimb = Nt_count + 1

        mass1 = p[CH[IH[i1].iCH].p[1]].m
        mass2 = p[CH[IH[i1].iCH].p[2]].m

        lRm0 = isa(lRm, Int) ? lRm : lRm[IH[i1].iCH]
        temp = real(Ec) < (mass1 + mass2)
        Rm = !((lRm0 == 0 || lRm0 == 1) && temp)

        if Rm
            kon = sqrt((Ec^2 - (mass1 + mass2)^2) *
                       (Ec^2 - (mass1 - mass2)^2)) / (2.0 * Ec)

            # Np 个 offshell
            for idim in 1:Np
                Nt_count += 1
                Dim[1, Nt_count] = i1
                Dim[2, Nt_count] = idim
                Dim[3, Nt_count] = 0
            end

            IH[i1].k = kon
            IH[i1].w = 0.0

            # onshell
            Nt_count += 1
            Dim[1, Nt_count] = i1
            Dim[2, Nt_count] = Np + 1
            Dim[3, Nt_count] = 1

        else
            for idim in 1:Np
                Nt_count += 1
                Dim[1, Nt_count] = i1
                Dim[2, Nt_count] = idim
                Dim[3, Nt_count] = 0
            end
        end

        IH[i1].Dime = Nt_count
    end

    return SYS, IH, Dim
end
#------------------------------------------------------------------------------------------
# Potential kernel 
# Form factor. ->fV
function FFre(k, cutoffi, cutofff; cutoff_re_type=:Lambda, CHi=nothing, CHf=nothing, key_ex="V")

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
    if key_ex != "V" #not contact interaction
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
# Define interaction potentials. ->kernel
function fV(k, l, SYS, IA0, CHf, CHi, VVertex)
    """
    Arguments:
       k: momenta (4-vectors for initial/final states)
       l: helicities of particles
       SYS: system identifier (e.g., "KN" for kaon-nucleon)
       IA0: interaction amplitude object containing coupling info
       CHf, CHi: final and initial channel objects (contain cutoff parameters)
    """

    fV = 0.
    eps0 = @SVector zeros(ComplexF64, 5)
    U0 = @SVector zeros(ComplexF64, 4)
    U30 = SVector{5}(U0 for _ in 1:5)

    # 内联定义处理函数
    @inline function get_pol(J, Jh, k, l; star=false, bar=false)
        if J == 0
            (eps=eps0, U=U0, U3=U30)
        elseif J == 1 && Jh == 1
            (eps=FR.eps(k, l, star=star), U=U0, U3=U30)
        elseif J == 1 && Jh == 2
            (eps=eps0, U=FR.U(k, l, bar=bar), U3=U30)
        elseif J == 3 && Jh == 2
            (eps=eps0, U=U0, U3=FR.U3(k, l; bar=bar))
        else
            error("Spin $(J)/$(Jh) handling is pending implementation")
        end
    end
    # 一次性计算所有极化向量
    pol_f1, pol_f2, pol_i1, pol_i2 =
        get_pol(CHf.J[1], CHf.Jh[1], k.f1, l.f1; star=true, bar=true),  # pol_f1
        get_pol(CHf.J[2], CHf.Jh[2], k.f2, l.f2; star=true, bar=true),  # pol_f2
        get_pol(CHi.J[1], CHi.Jh[1], k.i1, l.i1),            # pol_i1
        get_pol(CHi.J[2], CHi.Jh[2], k.i2, l.i2)            # pol_i2



    cutoffi, cutofff = CHi.cutoff, CHf.cutoff

    KerV = 0.0  # Initialize kernel potential
    for le in 1:IA0.Nex  # Loop over all exchanged particles
        key_ex = IA0.key_ex[le]  # Name of exchanged particle
        if key_ex == "V"

            # Calculate potential contribution:
            KerV = VVertex(k, pol_f1, pol_f2, pol_i1, pol_i2, SYS, IA0, CHi, CHf)

        else
            J_ex, Jh_ex, m_ex = p[key_ex].J, p[key_ex].Jh, p[key_ex].m  # Quantum numbers and mass

            if IA0.dc[le] == 1  # Direct diagram
                vertex_if1 = (CHi.p[1], CHf.p[1]) # Vertex 1 name (e.g., (:D,:D))
                vertex_if2 = (CHi.p[2], CHf.p[2]) # Vertex 2 name
                k1i, k1f, k2i, k2f = k.i1, k.f1, k.i2, k.f2  # Momenta assignment
                m1i, m1f, m2i, m2f = CHi.m[1], CHf.m[1], CHi.m[2], CHf.m[2]  # Masses
                pol_1i, pol_1f, pol_2i, pol_2f = pol_i1, pol_f1, pol_i2, pol_f2  # Polarizations

            elseif IA0.dc[le] == 2  # Crossed diagram
                vertex_if1 = (CHi.p[1], CHf.p[2]) # Vertex 1 name (e.g., (:D,:D))
                vertex_if2 = (CHi.p[2], CHf.p[1])# Vertex 2 name
                k1i, k2f, k2i, k1f = k.i1, k.f1, k.i2, k.f2  # Crossed momenta
                m1i, m2f, m2i, m1f = CHi.m[1], CHf.m[1], CHi.m[2], CHf.m[2]  # Crossed masses
                pol_1i, pol_2f, pol_2i, pol_1f = pol_i1, pol_f1, pol_i2, pol_f2  # Crossed polarizations
            end

            k.q = k2f - k2i  # Momentum transfer
            k.q2 = -abs(k.q * k.q)  # Squared momentum transfer

            # Calculate vertices (both vector and scalar parts)
            # Here the q is from upper vertex 1 to lower vertex 2 hence a minus sign appears in Vertex 1 for q
            Ver1V, Ver1S = VVertex(k1i, m1i, pol_1i, k1f, m1f, pol_1f, -k.q, vertex_if1, key_ex, m_ex, anti=CHi.anti[1] | CHf.anti[1])
            Ver2V, Ver2S = VVertex(k2i, m2i, pol_2i, k2f, m2f, pol_2f, k.q, vertex_if2, key_ex, m_ex, anti=CHi.anti[2] | CHf.anti[2])

            # Compute propagator with form factor

            FF = FFre(k, cutoffi, cutofff, cutoff_re_type=SYS.cutoff_re_type, key_ex=key_ex) *
                 propFFex(k, key_ex, SYS.cutoff_ex, cutoff_ex_type=SYS.cutoff_ex_type, FF_ex_type=SYS.FF_ex_type)

            if J_ex == 0  # Scalar/pseudoscalar exchange
                KerV -= Ver1S * Ver2S * IA0.Ff[le] * FF  # Scalar-scalar potential
            end

            if J_ex == 1 && Jh_ex == 1  # Vector exchange
                # Full vector propagator with current conservation term
                KerV += ((Ver1V * Ver2V) - (Ver1V * k.q) * (Ver2V * k.q) / m_ex^2) * IA0.Ff[le] * FF
            end
        end

        fV = KerV  # Set final potential value

    end

    return fV  # Return computed potential
end
#potential kernel. -> VGI
function kernel(kf, ki, Ec, qn, SYS, IA, CH, IHf, IHi, VVertex)::ComplexF64 # Calculating kernel
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
        Kerm1 = fV(k, l, SYS, IA0, CHf, CHi, VVertex) * SYS.d[i][Int64(lJJ + l21i)+1, lf] * eta
        l.f1, l.i2 = -l.f1, -l.i2
        l.i1, l.f1 = -l.i1, -l.f1
        Kerp1 = fV(k, l, SYS, IA0, CHf, CHi, VVertex) * SYS.d[i][Int64(lJJ - l21i)+1, lf]
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
#------------------------------------------------------------------------------------------
#propagator. -> VGI
function propagator(k, kv, w, wv, Ec, Np, CH, IH0, Dimo, lRm, eps)
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

    Ec += eps

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
#------------------------------------------------------------------------------------------
# calculate V G I  -> res0
function fVGI(Ec, qn, SYS, IA, CH, IH, VVertex, lRm, eps)
    # Determine the dimension of the work matrix. 
    Np = length(SYS.kv)
    SYS, IH, Dim, = workSpace(Ec, lRm, Np, SYS, CH, IH)
    Nt = size(Dim, 2)

    # Initialize matrices
    Vc = zeros(Complex{Float64}, Nt, Nt)
    Gc = zeros(Complex{Float64}, Nt, Nt)
    II = zeros(Float64, Nt, Nt)

    # Populate Vc, Gc, and II
    for i_f in 1:Nt
        IHf = IH[Dim[1, i_f]]
        kf = Dim[3, i_f] == 1 ? IHf.k : real(SYS.kv[Dim[2, i_f]])
        for i_i in i_f:Nt
            IHi = IH[Dim[1, i_i]]
            ki = Dim[3, i_i] == 1 ? IHi.k : real(SYS.kv[Dim[2, i_i]])
            # Compute Vc for upper triangular part
            Vc[i_f, i_i] = kernel(kf, ki, Ec, qn, SYS, IA, CH, IHf, IHi, VVertex)::ComplexF64
            #@show Vc[i_f, i_i] 
            #exit()

            # Compute Gc and II only on the diagonal
            if i_i == i_f
                w = Dim[3, i_f] == 1 ? 0. : real(SYS.wv[Dim[2, i_f]])
                Gc[i_f, i_i] = propagator(kf, SYS.kv, w, SYS.wv, Ec, Np, CH, IHf, Dim[3, i_f], lRm, eps)
                II[i_f, i_i] = 1.0
            end
        end
    end


    # Populate lower triangular part of Vc based on conjugate symmetry
    for i_f in 2:Nt
        Dimf_IH = IH[Dim[1, i_f]]

        for i_i in 1:i_f-1
            Dimi_IH = IH[Dim[1, i_i]]

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
#*******************************************************************************************
# RESCATTERING AMPLITUDE  AND POLES 
#*******************************************************************************************
# Calculate rescattering amplitude and  serach poles.
# Calculate the \sum|M|^2 for all channels ->resc0
function M2_channel(T, CH, IH)

    NCH = length(CH)
    resM2 = zeros(Float64, NCH, NCH)
    for iCHf in 1:NCH
        for iCHi in 1:NCH
            for f in CH[iCHf].IHb:CH[iCHf].IHe, i in CH[iCHi].IHb:CH[iCHi].IHe
                #if IH[f].Dimo == 1 && IH[i].Dimo == 1
                resM2[iCHf, iCHi] += abs2(T[IH[f].Dime, IH[i].Dime])
                #end
            end
        end
    end

    return resM2
end
# Show  informations ->res
function showSYSInfo(Range, qn, IA, CH, IH)

    dashline = repeat('-', 90)
    Nc = eachindex(CH)
    println(dashline)
    @printf("%-16s%-15s\n", "channel","number of exchanges", )

    for i1 in Nc
        @printf("%-15s", String(CH[i1].p[1]) * ":" * String(CH[i1].p[2]))
        for i in Nc
            @printf("%2d ", IA[i1, i].Nex)
        end
        @printf(";  cutoff = %5.3f GeV\n", CH[i1].cutoff)
    end
    println(dashline)
         P = if isdefined(qn, :P) && !isnothing(qn.P)
           qn.P == 1 ? "+" : qn.P == -1 ? "-" : " "
       else
           " "
       end
     C = if isdefined(qn, :C) && !isnothing(qn.C)
           qn.C == 1 ? "+" : qn.C == -1 ? "-" : " "
       else
           " "
       end
    println("I(J,P,C)=$(qn.I)/$(qn.Ih)($(qn.J)/$(qn.Jh),$P,$C): independent helicities")
    for ih in eachindex(IH)
        @printf("%-15s: %2d/%1d, %2d/%1d \n", String(CH[IH[ih].iCH].p[1]) * ":" * String(CH[IH[ih].iCH].p[2]), IH[ih].hel[1],
            IH[ih].helh[1], IH[ih].hel[2], IH[ih].helh[2])
    end

    println(dashline)
    println("ER=$(Range.ERmax) to $(Range.ERmin) GeV, NER=$(Range.NER)  ;   EI=$(Range.EIt*1e3) MeV, NEI=$(Range.NEI)")
    println(dashline)
end
function showPoleInfo(qn, Ec, reslog, filename)
    open(filename, "w") do f
        # 打开文件写入模式会自动清空内容
    end
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


    dashline = repeat('-', 90)
    println(dashline)
     P = if isdefined(qn, :P) && !isnothing(qn.P)
           qn.P == 1 ? "+" : qn.P == -1 ? "-" : " "
       else
           " "
       end
     C = if isdefined(qn, :C) && !isnothing(qn.C)
           qn.C == 1 ? "+" : qn.C == -1 ? "-" : " "
       else
           " "
       end
    println("I(J,P,C)=$(qn.I)/$(qn.Ih)($(qn.J)/$(qn.Jh),$P,$C)  pole= $Ampmin at $(Ampminx * 1e3), $Ampminy")
    println(dashline)
    return Ampmin, Ampminx, Ampminy
end
# Calculate the rescattering process   ->resc
function resc0(Range, iER, qn, SYS, IA, CH, IH, VVertex, eps)
    # ⚠ 复制 IH 防止并行修改共享引用
    IH_local = deepcopy(IH)

    Ect = ComplexF64[]        # 复数能量
    reslogt = Float64[]       # log|1 - VG|
    Dim = nothing
    TG = nothing
    NCH = length(CH)
    resM2 = zeros(Float64, NCH, NCH)

    # 计算 ER
    ER = Range.ERmax - (iER - 1) * (Range.ERmax - Range.ERmin) / (Range.NER - 1)

    # 如果使用 PL 方案，则转换 ER
    if Range.Ep[1] == "L"
        PL = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / (Range.NER - 1)
        ER = sqrt((sqrt(PL^2 + p[Range.Ep[2]].m^2) + p[Range.Ep[3]].m)^2 - PL^2)
    end

    # 循环虚部 EI
    NEI = Range.NEI
    for iEI in -NEI:NEI
        EI = (NEI != 0) ? iEI * Range.EIt / NEI : 0.0
        Ec = ER + EI * im

        # 调用 fVGI 时使用本地 IH_local，避免 race condition
        Vc, Gc, II, IH_local, Dim = fVGI(Ec, qn, SYS, IA, CH, IH_local, VVertex, qn.lRm, eps)
        VGI = II - Vc * Gc
        detVGI = det(VGI)
        logdetVGI = log(abs(detVGI)^2)

        push!(Ect, Ec)
        push!(reslogt, logdetVGI)

        # 仅在 EI=0 时计算 T、TG 和 resM2
        if iEI == 0
            T = inv(VGI) * Vc
            TG = T * Gc
            resM2 = M2_channel(T, CH, IH_local)
        end
    end

    # 如果使用 PL，替换第一个能量
    if Range.Ep[1] == "L"
        Ect[1] = PL
    end

    return Ect, reslogt, resM2, IH_local, Dim, TG
end
# parallel calculation 
function worker_resc(r_range, Range, iER_start, qn, SYS, IA, CH, IH, VVertex, eps, progressbar)
    # r_range 是迭代范围，iER_start 是起始索引（用于回调）

    Ec_local = ComplexF64[]
    ER_local = Float64[]
    reslog_local = Float64[]
    resM2_local = Matrix{Float64}[]  # 应该是 Vector{Matrix{Float64}}
    TGt_local = Matrix{ComplexF64}[]
    IHt_local = Vector{structIndependentHelicity}[]
    Dimt_local = Matrix{Int64}[]

    # 进度条回调函数（仅在工作进程2且需要进度条时有效）
    function progress_callback(pb, idx)
        ProgressBars.update(pb)
    end

    # 根据工作进程ID和progressbar标志决定是否创建进度条
    if myid() == 2 && progressbar
        # 需要知道总迭代次数，这里使用r_range的长度
        pb = ProgressBar(collect(r_range))  # 创建进度条
        callback = i -> progress_callback(pb, i)  # 回调函数
    else
        callback = _ -> nothing
    end

    for (local_idx, iER) in enumerate(r_range)
        res0 = resc0(Range, iER, qn, SYS, IA, CH, deepcopy(IH), VVertex, eps)

        append!(Ec_local, res0[1])
        append!(reslog_local, res0[2])
        push!(resM2_local, res0[3])
        push!(ER_local, real(res0[1][1]))
        push!(IHt_local, res0[4])
        push!(Dimt_local, res0[5])
        push!(TGt_local, res0[6])

        # 回调更新进度条，传入全局索引
        callback(iER_start + local_idx - 1)
    end

    # 返回块起始索引 + 本块结果，确保顺序可控
    return (iER_start, Ec_local, reslog_local, resM2_local,
        ER_local, IHt_local, Dimt_local, TGt_local)
end
function paddingshare(TGt_vec::Vector{Matrix{ComplexF64}})
    # 找出最大维度
    max_dim = maximum(size(m, 1) for m in TGt_vec)
    N_total = length(TGt_vec)
    TGt_padded = zeros(ComplexF64, max_dim, max_dim, N_total)
    # 记录实际矩阵维度
    valid_dims = Vector{Int}(undef, N_total)
    # 多线程填充
    Threads.@threads for i in 1:N_total
        src = TGt_vec[i]
        r, c = size(src)
        valid_dims[i] = r
        copyto!(view(TGt_padded, 1:r, 1:c, i), src)
    end
    # 清理原向量
    TGt_vec = nothing
    TGt_shared = SharedArray(TGt_padded)
    TGt_padded = nothing
    return TGt_shared, valid_dims
end
function dim_to_3d_filled(Dimt::Vector{Matrix{Int}})
    n_vectors = length(Dimt)
    if n_vectors == 0
        return Array{Int}(undef, 3, 0, 0)
    end
    # 找最大 Nt
    max_len = maximum(size(m, 2) for m in Dimt)
    # 分配
    result = zeros(Int, 3, max_len, n_vectors)
    @inbounds for k in eachindex(Dimt)
        mat = Dimt[k]          # 3 × Nt
        Nt = size(mat, 2)
        result[:, 1:Nt, k] .= mat
    end
    Dim = SharedArray(result)
    return Dim
end
function resc(Sys, qn, Range, channels, Ff, cutoff, VVertex; Np=10, Nx=10, Nphi=5, eps=+1e-4im, progressbar=true)

    # -------------------------
    # 初始化系统
    # -------------------------
    SYS, IA, CH, IH = preprocessing(Sys, qn, channels, Ff, cutoff, Np, Nx, Nphi)

    Ec = ComplexF64[]
    ER = Float64[]
    reslog = Float64[]
    resM2 = Matrix{Float64}[]  # Vector{Matrix{Float64}}
    TGt = Matrix{ComplexF64}[]
    IHt = Vector{structIndependentHelicity}[]
    Dimt = Matrix{Int}[]

    showSYSInfo(Range, qn, IA, CH, IH)

    # -------------------------
    # 分块并行设置
    # -------------------------
    Ntot = Range.NER
    num_workers = nworkers()
    nevt_per_worker = div(Ntot, num_workers)
    # 确保覆盖所有索引
    ranges = []
    start_indices = []
    for i in 0:(num_workers-1)
        start_idx = i * nevt_per_worker + 1
        end_idx = (i == num_workers - 1) ? Ntot : (i + 1) * nevt_per_worker
        push!(ranges, start_idx:end_idx)
        push!(start_indices, start_idx)
    end

    # 使用pmap并行计算
    results = pmap(1:num_workers) do worker_id
        worker_resc(ranges[worker_id], Range, start_indices[worker_id],
            qn, SYS, IA, CH, IH, VVertex, eps,
            progressbar && worker_id == 1)  # 只在第一个工作进程显示进度条
    end

    # -------------------------
    # 主进程合并 (按块 start 索引排序)
    # -------------------------
    sort!(results, by=x -> x[1])  # 按块起始 iER 排序
    for block in results
        append!(Ec, block[2])
        append!(reslog, block[3])
        append!(resM2, block[4])
        append!(ER, block[5])
        append!(IHt, block[6])
        append!(Dimt, block[7])
        append!(TGt, block[8])
    end

    # -------------------------
    # 后处理
    # -------------------------
    # 确保输出目录存在
    mkpath("res")

    open("res/" * Sys * ".txt", "w") do f
    end
    showPoleInfo(qn, Ec, reslog, "res/output.txt")
    
    #填充成array
    TGt = paddingshare(TGt)
    Dimt = dim_to_3d_filled(Dimt)

    return Ec, reslog, resM2, ER, IHt, Dimt, TGt, SYS, IA, CH, IH
end
#------------------------------------------------------------------------------------------
# functions used to calculate the cross section.
#* lambda function 
@inline function lambda(m1, m2, m3)
    m1_sq = m1 * m1
    m2_plus_m3 = m2 + m3
    m2_minus_m3 = m2 - m3
    l = (m1_sq - m2_plus_m3 * m2_plus_m3) * (m1_sq - m2_minus_m3 * m2_minus_m3)
    return l > 0.0 ? sqrt(l) : 0.0
end
#* Calculate cross section for all channels , only for 2-2 process
function simpleXsection(ER, M2, CH, qn; Ep=("cm",))
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


            if Ep[1] == "L"
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
#*******************************************************************************************
# DECAYS
#*******************************************************************************************
#* product the complete informations about the physcial process.
function proc(pf, pin, amps)
    mf = [p[p0].m for p0 in pf]
    namef = [p[p0].nameL for p0 in pf]
    mi = [p[p0].m for p0 in pin]
    namei = [p[p0].nameL for p0 in pin]
    ranges = Tuple([-p[p0].J:p[p0].Jh:p[p0].J for p0 in pf])
    ch = (pf=pf, namef=namef, mf=mf, pin=pin, namei=namei, mi=mi, amps=amps, ranges=ranges)
end
#------------------------------------------------------------------------------------------
# Get the corresponding IH and Dim for the given case. Note: Interpolation is not performed here simultaneously, as it may cause memory leaks for unknown reasons. -> setTGA
function IHDim(E, Et, Range, TGt, IHt, Dimt)
    ii = Range.NER - Xs.Nsij(E, Range.ERmin, Range.ERmax, Range.NER - 1)
    Emin, Emax = Et[ii+1], Et[ii]
    dmin, dmax = TGt[2][ii+1], TGt[2][ii]
    if dmin == dmax
        #Tmin, Tmax = TGt[ii+1], TGt[ii]
        #if size(Tmin) == size(Tmax)

        return IHt[ii], Dimt[:, :, ii]
    else
        mid = 0.5 * (Emin + Emax)
        idx = (E < mid) + 1
        return IHt[ii+idx-1], Dimt[:, :, ii+idx-1]
    end
end
# Frame transition from the static frame of total system  to  cms of two partilce. -> setTGA
#############################################################################
function LorentzBoost(k::SVector{5,Float64}, p::SVector{5,Float64})
    kp = k[1] * p[1] + k[2] * p[2] + k[3] * p[3]  # 点积 k ⋅ p
    p5 = sqrt(p[4]^2 - p[1]^2 - p[2]^2 - p[3]^2)
    beta_factor = kp / (p[4] + p5) + k[4]       # β-related factor
    inv_p5 = 1.0 / p5


    k1 = k[1] + p[1] * beta_factor * inv_p5
    k2 = k[2] + p[2] * beta_factor * inv_p5
    k3 = k[3] + p[3] * beta_factor * inv_p5
    k4 = (p[4] * k[4] + kp) * inv_p5
    k5 = k[5]

    return @SVector [k1, k2, k3, k4, k5]
end
function LorentzBoost(momenta::Vector{SVector{5,Float64}}, p::SVector{5,Float64})
    return [LorentzBoost(i, p) for i in momenta]
end
function Rotation(k::SVector{5,Float64}, ct::Float64, st::Float64, cp::Float64, sp::Float64)
    k1 = ct * cp * k[1] + ct * sp * k[2] - st * k[3]
    k2 = -sp * k[1] + cp * k[2]
    k3 = st * cp * k[1] + st * sp * k[2] + ct * k[3]
    return @SVector [k1, k2, k3, k[4], k[5]]
end
function Rotation(momenta::Vector{SVector{5,Float64}}, ct::Float64, st::Float64, cp::Float64, sp::Float64)
    return [Rotation(k, ct, st, cp, sp) for k in momenta]
end
function LorentzBoostRotation(k, tecm, p1, p2)
    kij = k[p1] + k[p2]
    sij = kij * kij
    P = @SVector [0.0, 0.0, 0.0, tecm, tecm]
    pLB = @SVector [-kij[1], -kij[2], -kij[3], kij[4], sqrt(sij)]
    knew = LorentzBoost(k, pLB)
    Pnew = LorentzBoost(P, pLB)
    ct, st, cp, sp = FR.kph(knew[p2])
    Pnew = Rotation(Pnew, ct, st, cp, sp)
    knew = Rotation(knew, ct, st, cp, sp)
    return knew, Pnew
end
#* set the frame and other things for calculating TGA    The output will flow to -> TGA
function setTGA(para, k, tecm, i, j)
    wij = sqrt((k[i] + k[j]) * (k[i] + k[j]))
    IHt, Dimt, TGt, Et, Range = para.IHt, para.Dimt, para.TGt, para.Et, para.Range
    IH, Dim = IHDim(wij, Et, Range, TGt, IHt, Dimt)
    kn, Pn = LorentzBoostRotation(k, tecm, i, j) #reference frame thansformation
    new = (E=wij, IH=IH, Dim=Dim, k=kn, P=Pn, resc=(i, j))
    return merge(para, new)
end
# auxiliary fuction for TGA to generate the full index for all final particles. ->TGA
function build_full_idx(idx_keep::Vector{Int}, drop::Vector{Int}, drop_vals::Vector{Int}, N::Int)
    full = Vector{Int}(undef, N)

    # 直接使用向量元素，无需拆箱
    d1, d2 = drop
    v1, v2 = drop_vals

    ikeep = 1
    for i in 1:N
        if i == d1
            full[i] = v1
        elseif i == d2
            full[i] = v2
        else
            full[i] = idx_keep[ikeep]
            ikeep += 1
        end
    end

    return full  # 返回 Vector{Int} 而不是 Tuple
end
#* calculate TGA
function TGA(para, cfinal, cinter, ranges)

    ranges = collect(ranges)
    resc = collect(para.resc)

    E, IH, Dim, k, P, Et, TGt, Range = para.E, para.IH, para.Dim, para.k, para.P, para.Et, para.TGt, para.Range
    qn, CH, SYS = para.qn, para.CH, para.SYS

    # Interpolation: placing interpolation here helps prevent memory leaks
    ii = Range.NER - Xs.Nsij(E, Range.ERmin, Range.ERmax, Range.NER - 1)
    Emin, Emax = Et[ii+1], Et[ii]
    dmin, dmax = TGt[2][ii+1], TGt[2][ii]
    Tmin, Tmax = view(TGt[1], :, :, ii + 1), view(TGt[1], :, :, ii)
    TeT = dmin == dmax
    #Tmin, Tmax = TGt[ii+1], TGt[ii]
    #TeT = size(Tmin) == size(Tmax)

    if TeT
        inv_dE = 1.0 / (Emax - Emin)
        ww = (E - Emin) * inv_dE

    else
        mid = 0.5 * (Emin + Emax)
        idx = (E < mid) + 1
        TG = view(TGt[1], :, :, ii + idx - 1)
        #TG = TGt[ii+idx-1]

    end

    # Independent helicities for final states
    CHf = CH[cfinal isa String ? SYS.channel[cfinal] : cfinal]
    IHfb, IHfe = CHf.IHb, CHf.IHe
    # Rank of dimensions for intermediate state

    lJJ = qn.J / qn.Jh
    TGAdic = Dict()
    ndx = length(SYS.xv)
    ndphi = length(SYS.pv)
    #loop for all Independent helicities of final state
    keep = filter(i -> !(i in resc), eachindex(ranges))
    left = ntuple(i -> ranges[keep[i]], length(keep))
    @inbounds for iIHf in IHfb:IHfe
        IHf = IH[iIHf]
        l21 = IHf.hel[2] / IHf.helh[2] - IHf.hel[1] / IHf.helh[1]
        dl21 = Int64(lJJ + l21) + 1
        eta = 0.0
        h1, h2 = IHf.hel[1], IHf.hel[2]
        #loop for dimensions for intermediate state
        for l in Iterators.product(left...)
            l = collect(l)
            TGA = 0.0
            for cinter0 in cinter
                #CHc = CH[SYS.channel[cinter0.ch]]
                ch = cinter0.ch
                CHc = CH[ch isa String ? SYS.channel[ch] : ch]
                m1, m2 = CHc.m[1], CHc.m[2]
                @inbounds for iIHc in CHc.IHb:CHc.IHe
                    IHc = IH[iIHc]
                    l21c = IHc.hel[2] / IHc.helh[2] - IHc.hel[1] / IHc.helh[1]
                    etap = CHc.P[1] * CHc.P[2] * qn.P * (-1)^(qn.J / qn.Jh - CHc.J[1] / CHc.Jh[1] - CHc.J[2] / CHc.Jh[2])
                    dl21pp = Int64(lJJ + l21c) + 1
                    dl21pm = Int64(lJJ - l21c) + 1
                    la, lb = IHc.hel[1], IHc.hel[2]
                    lp = build_full_idx(l, resc, [la, lb], length(ranges))
                    lm = build_full_idx(l, resc, [-la, -lb], length(ranges))
                    @inbounds for iDim in IHc.Dimb:IHc.Dime

                        kc = Dim[3, iDim] == 1 ? real(IHc.k) : real(SYS.kv[Dim[2, iDim]])
                        Eon = sqrt(kc^2 + (m1 > m2 ? m1^2 : m2^2))  # 较重粒子能量
                        E1 = m1 > m2 ? Eon : E - Eon                # 粒子1能量
                        E2 = m1 > m2 ? E - Eon : Eon                # 粒子2能量
                        # loop for integration of angles
                        A = 0.0
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
                                px, py, pz = kc * st * cosphi, kc * st * sinphi, kc * ct
                                p1 = SVector{5,Float64}(-px, -py, -pz, E1, m1)
                                p2 = SVector{5,Float64}(px, py, pz, E2, m2)
                                k[resc[1]] = p1
                                k[resc[2]] = p2
                                Ap = cinter0.Vertex(k, P, lp, cinter0.cached)
                                Am = cinter0.Vertex(k, P, lm, cinter0.cached)
                                Dp = SYS.d[i][dl21, dl21pp]
                                Dm = SYS.d[i][dl21, dl21pm] * etap
                                A -= (Dp * Ap + Dm * Am) * exppl * SYS.wxv[i] * SYS.wpv[i]
                            end

                        end

                        if TeT
                            TG0 = Tmin[IHf.Dime, iDim] * (1.0 - ww) + Tmax[IHf.Dime, iDim] * ww #插值
                        else
                            TG0 = TG[IHf.Dime, iDim]
                        end

                        TGA += (A * TG0) * cinter0.weight  #note that the weights and p20^2/(2pi)^3 is in G0
                    end
                end

            end
            lp = Tuple(build_full_idx(l, resc, [h1, h2], length(ranges)))
            lm = Tuple(build_full_idx(l, resc, [-h1, -h2], length(ranges)))

            TGA *= (2.0 * lJJ + 1.0) / (8.0 * pi)
            TGAdic[lp] = TGA
            # extend from independent helicities to all helicities        
            eta = CHf.P[1] * CHf.P[2] * qn.P * (-1)^(qn.J / qn.Jh - CHf.J[1] / CHf.Jh[1] - CHf.J[2] / CHf.Jh[2])
            TGAdic[lm] = TGA * eta
            if h1 == 0 && h2 == 0
                TGAdic[lp] = TGA * sqrt(2.0)
            end

        end
    end
    return TGAdic
end

end
