
module GEN

export GENEV
using StaticArrays

function cut(a::Int32)::Int64
    binary_string = bitstring(a * 69069)
    binary_str = binary_string[end-31:end]
    if startswith(binary_str, "0")
        decimal_num = parse(Int, binary_str, base=2)
    else
        decimal_num = parse(Int, binary_str[2:end], base=2)
        decimal_num -= 2^(length(binary_str) - 1)
    end
    return decimal_num
end

global MCGN::Int32 = 12345
function RNDM()::Float32
    global MCGN
    MSK1::Int32 = 0xC000000
    MSK2::Int32 = 0x33333300
    AMAN::Float32 = 0.0f0
    MCGN = cut(MCGN::Int32)
    MANT = MCGN >>> 8
    if MANT != 0
        AMAN = MANT
        MANT = reinterpret(Int32, AMAN)
        MANT = MANT - MSK1
        AMAN = reinterpret(Float32, MANT)
        RNDM::Float32 = AMAN
    elseif MANT == 0
        MANT = MSK2
        AMAN = reinterpret(Float32, MANT)
        RNDM = AMAN
    end
    return RNDM
end

function NRAN(N::Int64)::Vector{Float64}
    return rand(N)

    VECTOR = Vector{Float64}(undef, N)
    for I in 1:N
        VECTOR[I] = RNDM()
    end
    return VECTOR
end

function RANGNR(vector_length::Int64, num_elements::Int64)
    rno = NRAN(vector_length)
    if num_elements - 2 == 0
        if rno[1] > rno[2]
            rno[1], rno[2] = rno[2], rno[1]
        end
    elseif num_elements - 2 > 0
        km1 = num_elements - 1
        @inbounds @simd for i in 1:km1
            ni = num_elements - i
            for j in 1:ni
                if rno[j] - rno[j+1] < 0
                    rno[j], rno[j+1] = rno[j+1], rno[j]
                    return
                end
            end
        end
    end
    rno
end

function ROTES2!(cos_theta::Float64, sin_theta::Float64, cos_theta2::Float64, sin_theta2::Float64, pr0::MMatrix{5,18,Float64,90}, i::Int64)
    pr = reshape(pr0, 18 * 5)
    k1 = 5 * i - 4
    k2 = k1 + 1
    sa = pr[k1]
    sb = pr[k2]
    a = sa * cos_theta - sb * sin_theta
    pr[k2] = sa * sin_theta + sb * cos_theta
    k2 += 1
    b = pr[k2]
    pr[k1] = a * cos_theta2 - b * sin_theta2
    pr[k2] = a * sin_theta2 + b * cos_theta2
    pr0 = reshape(pr, 5, 18)
end

function PDK(A::Float64, B::Float64, C::Float64)::Float64
    A_squared = A * A
    B_squared = B * B
    C_squared = C * C
    pdk = 0.5 * sqrt(abs(A_squared + (B_squared - C_squared)^2 / A_squared - 2.0 * (B_squared + C_squared)))
    return pdk
end

function GENEV(tecm::Float64, EM::Vector{Float64})
    NT = length(EM)
    EMM = @MVector zeros(Float64, 18)
    PCM = @MArray zeros(Float64, 5, 18)
    EMS = @MArray zeros(Float64, 5, 18)
    PD = @MArray zeros(Float64, 5, 18)
    SM = @MVector zeros(Float64, 18)
    FFQ = SVector{18,Float64}([0.0, 3.141592, 19.73921, 62.01255, 129.8788, 204.0131, 256.3704, 268.4705,
        240.9780, 189.2637, 132.1308, 83.0202, 47.4210, 24.8295, 12.0006, 5.3858,
        2.2560, 0.8859])
    ETC = 0.0
    INIT = 0
    TWOPI = 6.2831853073
    WTMAX = 0.0
    WTMAXQ = 0.0
    KGENEV = 2

    if tecm - ETC != 0
        if INIT <= 0
            INIT = 1
            NTM1 = NT - 1
            NTM2 = NT - 2
            NTP1 = NT + 1
            NTNM4 = 3 * NT - 4
            EMM[1] = EM[1]
            TM = 0.0
            @inbounds @simd for i in 1:NT
                EMS[i] = EM[i]^2
                TM += EM[i]
                SM[i] = TM
            end
        end
        TECMTM = tecm - TM
        ETC = tecm
        EMM[NT] = tecm

        if KGENEV <= 1
            EMMAX = TECMTM + EM[1]
            EMMIN = 0.0
            WTMAX = 1.0
            for i in 2:NT
                EMMIN += EM[i-1]
                EMMAX += EM[i]
                WTMAX *= PDK(EMMAX, EMMIN, EM[i])
            end
            WTMAXQ = 1.0 / WTMAX
        else
            WTMAXQ = TECMTM^NTM2 * FFQ[NT] / tecm
        end
    end

    RNO = RANGNR(NTNM4, NTM2)

    if NTM2 > 0
        for j in 2:NTM1
            EMM[j] = RNO[j-1] * TECMTM + SM[j]
        end
    end

    if NTM2 >= 0
        WT = WTMAXQ
        IR = NTM2
        @inbounds for i in 1:NTM1
            PD[i] = PDK(EMM[i+1], EMM[i], EM[i+1])
            WT *= PD[i]
        end

        PCM[1, 1] = 0.0
        PCM[2, 1] = PD[1]
        PCM[3, 1] = 0.0

        @inbounds @simd for i in 2:NT
            PCM[1, i] = 0.0
            PCM[2, i] = -PD[i-1]
            PCM[3, i] = 0.0
            IR += 1
            BANG = TWOPI * RNO[IR]
            CB = cos(BANG)
            SB = sin(BANG)
            IR += 1
            C = 2.0 * RNO[IR] - 1.0
            S = sqrt(1.0 - C^2)

            if i != NT
                ESYS = sqrt(PD[i]^2 + EMM[i]^2)
                BETA = PD[i] / ESYS
                GAMA = ESYS / EMM[i]
                for j in 1:i
                    AA = PCM[1, j]^2 + PCM[2, j]^2 + PCM[3, j]^2
                    PCM[5, j] = sqrt(AA)
                    PCM[4, j] = sqrt(AA + EMS[j])
                    ROTES2!(C, S, CB, SB, PCM, j)
                    PSAVE = GAMA * (PCM[2, j] + BETA * PCM[4, j])
                    PCM[2, j] = PSAVE
                end
            else
                for j in 1:i
                    AA = PCM[1, j]^2 + PCM[2, j]^2 + PCM[3, j]^2
                    PCM[5, j] = sqrt(AA)
                    PCM[4, j] = sqrt(AA + EMS[j])
                    ROTES2!(C, S, CB, SB, PCM, j)
                end
            end
        end
    end
    for i in 1:NT
    PCM[5, i] = EM[i]
    end 
    return PCM, WT
end



end
