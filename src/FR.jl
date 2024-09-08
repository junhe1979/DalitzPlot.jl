#############################################################################
# FR #
#############################################################################
module FR
#export III, GA, GS, epsilon, Uc, Ubc, LCV
using StaticArrays
const g = SMatrix{5,5,Float64}([
    -1.0 0.0 0.0 0.0 0.0;
    0.0 -1.0 0.0 0.0 0.0;
    0.0 0.0 -1.0 0.0 0.0;
    0.0 0.0 0.0 1.0 0.0;
    0.0 0.0 0.0 0.0 0.0])
const I = SMatrix{4,4,ComplexF64}([
    1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im;
    0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im;
    0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im;
    0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im])
const GA = [SMatrix{4,4,ComplexF64}([
        0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im;
        0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im;
        0.0+0.0im -1.0+0.0im 0.0+0.0im 0.0+0.0im;
        -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im]),
    SMatrix{4,4,ComplexF64}([
        0.0+0.0im 0.0+0.0im 0.0+0.0im 0.0-1.0im
        0.0+0.0im 0.0+0.0im 0.0+1.0im 0.0+0.0im
        0.0+0.0im 0.0+1.0im 0.0+0.0im 0.0+0.0im
        0.0-1.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im]),
    SMatrix{4,4,ComplexF64}([
        0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im
        -1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im]),
    SMatrix{4,4,ComplexF64}([
        1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im -1.0+0.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im 0.0+0.0im -1.0+0.0im]),
    SMatrix{4,4,ComplexF64}([
        0.0+0.0im 0.0+0.0im 1.0+0.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im 0.0+0.0im 1.0+0.0im
        1.0+0.0im 0.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im 1.0+0.0im 0.0+0.0im 0.0+0.0im])]
@inline function GS(k::SVector{5,ComplexF64})::SMatrix{4,4,ComplexF64}
    k1_im_k2 = k[1] + im * k[2]
    k1_min_im_k2 = k[1] - im * k[2]

    return @SMatrix [
        k[4] 0 -k[3] -k1_im_k2;
        0 k[4] -k1_min_im_k2 k[3];
        k[3] k1_min_im_k2 -k[4] 0;
        k1_im_k2 -k[3] 0 -k[4]
    ]
end

@inline function GS(k::SVector{5,Float64})::SMatrix{4,4,ComplexF64}
    k1_im_k2 = k[1] + im * k[2]
    k1_min_im_k2 = k[1] - im * k[2]

    return @SMatrix [
        k[4] 0 -k[3] -k1_im_k2;
        0 k[4] -k1_min_im_k2 k[3];
        k[3] k1_min_im_k2 -k[4] 0;
        k1_im_k2 -k[3] 0 -k[4]
    ]
end

@inline function kph(k::SVector{5,ComplexF64})
    zkx, zky, zkz, zk0, zm = real(k[1]), real(k[2]), real(k[3]), k[4], real(k[5])
    zkk = sqrt(zkx^2 + zky^2 + zkz^2)
    if zkk < 1e-30
        return zk0, zm, zkk, 1.0 + 0.0im, 1.0, 0.0, 1.0, 0.0
    end
    ct = zkz / zkk
    st = sqrt(1.0 - ct^2)
    cp, sp = abs(zkx) < 1e-30 && abs(zky) < 1e-30 ? (1.0, 0.0) : (zkx / zkk / st, zky / zkk / st)
    expp = cp + im * sp
    return zk0, zm, zkk, expp, ct, st, cp, sp
end
@inline function kph(k::SVector{5,Float64})
    zkx, zky, zkz, zk0, zm = k
    zkk = sqrt(zkx^2 + zky^2 + zkz^2)

    if zkk < 1e-30
        return zk0, zm, zkk, 1.0 + 0.0im, 1.0, 0.0, 1.0, 0.0
    end
    ct = zkz / zkk
    st = sqrt(1.0 - ct^2)
    cp, sp = abs(zkx) < 1e-30 && abs(zky) < 1e-30 ? (1.0, 0.0) : (zkx / zkk / st, zky / zkk / st)
    expp = cp + im * sp
    return zk0, zm, zkk, expp, ct, st, cp, sp
end
@inline function U(k, l::Int64; bar=false, V=false)::SVector{4,ComplexF64}
    zk0, zm, zkk, expp, ct, st, cp, sp = kph(k)
    zk0m = zk0 + zm
    zfac = sqrt(2.0 * zm * zk0m)
    ct2, st2 = sqrt((1.0 + ct) / 2.0) / zfac, sqrt((1.0 - ct) / 2.0) / zfac
    ct2dexp, ct2exp = ct2 / expp, ct2 * expp
    if !bar
        if !V
            return l == 1 ? SVector{4,ComplexF64}(zk0m * ct2dexp, zk0m * st2, zkk * ct2dexp, zkk * st2) :
                   SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2exp, zkk * st2, -zkk * ct2exp)
        else
            return l == 1 ? SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2) :
                   SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
        end
    else
        if !V
            return l == 1 ? SVector{4,ComplexF64}(zk0m * ct2exp, zk0m * st2, -zkk * ct2exp, -zkk * st2) :
                   SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2dexp, -zkk * st2, zkk * ct2dexp)
        else
            return l == 1 ? SVector{4,ComplexF64}(zkk * ct2exp, zkk * st2, -zk0m * ct2exp, -zk0m * st2) :
                   SVector{4,ComplexF64}(zkk * st2, -zkk * ct2dexp, zk0m * st2, -zk0m * ct2dexp)
        end
    end
end
@inline function eps(k, l::Int64; star=false)::SVector{5,ComplexF64}
    zk0, zm, zkk, expp, ct, st, cp, sp = kph(k)
    #eps = SVector{5,ComplexF64}(0.0, 0.0, 0.0, 0.0, 0.0)
    expdsqrt2 = expp * sqrt(2.0) / 2.0
    dexpdsqrt2 = 1.0 / expp * sqrt(2.0) / 2.0

    if !star
        if l == 1
            eps = SVector((im * sp - ct * cp) * dexpdsqrt2, (-im * cp - ct * sp) * dexpdsqrt2, st * expdsqrt2 + 0.0im, 0.0im, 0.0im)
        elseif l == 0
            eps = SVector(zk0 / zm * st * cp + 0.0im, zk0 / zm * st * sp + 0.0im, zk0 / zm * ct + 0.0im, zkk / zm + 0.0im, 0.0im)

        elseif l == -1
            eps = SVector((im * sp + ct * cp) * expdsqrt2, (-im * cp + ct * sp) * expdsqrt2, -st * expdsqrt2, 0.0im, 0.0im)
        end
    else
        if l == 1
            eps = SVector((-im * sp - ct * cp) * expdsqrt2, (im * cp - ct * sp) * expdsqrt2, st * expdsqrt2 + 0.0im, 0.0im, 0.0im)

        elseif l == 0
            eps = SVector(zk0 / zm * st * cp + 0.0im, zk0 / zm * st * sp + 0.0im, zk0 / zm * ct + 0.0im, zkk / zm + 0.0im, 0.0im)
        elseif l == -1
            eps = SVector((-im * sp + ct * cp) * dexpdsqrt2, (im * cp + ct * sp) * dexpdsqrt2, -st * dexpdsqrt2, 0.0im, 0.0im)
        end
    end

    return eps
end
@inline function U3(k, l::Int64; bar=false, V=false)::SVector{5,SVector{4, ComplexF64}}
    zk0, zm, zkk, expp, ct, st, cp, sp = kph(k)
    zk0m = zk0 + zm
    zfac = sqrt(2.0 * zm * zk0m)
    ct2, st2 = sqrt((1.0 + ct) / 2.0) / zfac, sqrt((1.0 - ct) / 2.0) / zfac
    ct2dexp, ct2exp = ct2 / expp, ct2 * expp
    sqrt2_3 = sqrt(2.0 / 3.0)
    sqrt1_3 = sqrt(1.0 / 3.0)
    expdsqrt2 = expp * sqrt(2.0) / 2.0
    dexpdsqrt2 = 1.0 / expp * sqrt(2.0) / 2.0
    zk0_zm = zk0 / zm
    # Calculate U3
    if !bar
        if !V
            if l == 1
                xu_m1 = SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2exp, zkk * st2, -zkk * ct2exp)
                xu_1 = SVector{4,ComplexF64}(zk0m * ct2dexp, zk0m * st2, zkk * ct2dexp, zkk * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2
                U3_1 = (sqrt2_3 * xe_01) * xu_1 + (sqrt1_3 * xe_11) * xu_m1
                U3_2 = (sqrt2_3 * xe_02) * xu_1 + (sqrt1_3 * xe_12) * xu_m1
                U3_3 = (sqrt2_3 * xe_03) * xu_1 + (sqrt1_3 * xe_13) * xu_m1
                U3_4 = (sqrt2_3 * xe_04) * xu_1 
            elseif l == -1
                xu_m1 = SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2exp, zkk * st2, -zkk * ct2exp)
                xu_1 = SVector{4,ComplexF64}(zk0m * ct2dexp, zk0m * st2, zkk * ct2dexp, zkk * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2
                U3_1 = (sqrt2_3 * xe_01) * xu_m1 + (sqrt1_3 * xe_m11) * xu_1
                U3_2 = (sqrt2_3 * xe_02) * xu_m1 + (sqrt1_3 * xe_m12) * xu_1
                U3_3 = (sqrt2_3 * xe_03) * xu_m1 + (sqrt1_3 * xe_m13) * xu_1
                U3_4 = (sqrt2_3 * xe_04) * xu_m1 
            elseif l == 3
                xu_1 = SVector{4,ComplexF64}(zk0m * ct2dexp, zk0m * st2, zkk * ct2dexp, zkk * st2)
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2

                U3_1 = xe_11 * xu_1
                U3_2 = xe_12 * xu_1
                U3_3 = xe_13 * xu_1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            elseif l == -3
                xu_m1 = SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2exp, zkk * st2, -zkk * ct2exp)
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2
                U3_1 = xe_m11 * xu_m1
                U3_2 = xe_m12 * xu_m1
                U3_3 = xe_m13 * xu_m1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            end
        else
            if l == 1
                xu_m1 = SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
                xu_1 = SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2
                U3_1 = (sqrt2_3 * xe_01) * xu_1 + (sqrt1_3 * xe_11) * xu_m1
                U3_2 = (sqrt2_3 * xe_02) * xu_1 + (sqrt1_3 * xe_12) * xu_m1
                U3_3 = (sqrt2_3 * xe_03) * xu_1 + (sqrt1_3 * xe_13) * xu_m1
                U3_4 = (sqrt2_3 * xe_04) * xu_1 
            elseif l == -1
                xu_m1 = SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
                xu_1 = SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2

                U3_1 = (sqrt2_3 * xe_01) * xu_m1 + (sqrt1_3 * xe_m11) * xu_1
                U3_2 = (sqrt2_3 * xe_02) * xu_m1 + (sqrt1_3 * xe_m12) * xu_1
                U3_3 = (sqrt2_3 * xe_03) * xu_m1 + (sqrt1_3 * xe_m13) * xu_1
                U3_4 = (sqrt2_3 * xe_04) * xu_m1 
            elseif l == 3
                xu_1 = SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2)
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2

                U3_1 = xe_11 * xu_1
                U3_2 = xe_12 * xu_1
                U3_3 = xe_13 * xu_1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            elseif l == -3
                xu_m1 = SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2
                xe_m14 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
                U3_1 = xe_m11 * xu_m1
                U3_2 = xe_m12 * xu_m1
                U3_3 = xe_m13 * xu_m1
                U3_4 = xe_m14 * xu_m1
            end
        end
        return SVector{5,SVector{4, ComplexF64}}(U3_1, U3_2, U3_3, U3_4, SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im))
        #return (U3_1, U3_2, U3_3, U3_4, SVector(0.0im, 0.0im, 0.0im, 0.0im))
    else
        if !V
            if l == 1
                xu_m1 = SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2dexp, -zkk * st2, zkk * ct2dexp)
                xu_1 = SVector{4,ComplexF64}(zk0m * ct2exp, zk0m * st2, -zkk * ct2exp, -zkk * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2

                U3_1 = (sqrt2_3 * xe_01) * xu_1 + (sqrt1_3 * xe_11) * xu_m1
                U3_2 = (sqrt2_3 * xe_02) * xu_1 + (sqrt1_3 * xe_12) * xu_m1
                U3_3 = (sqrt2_3 * xe_03) * xu_1 + (sqrt1_3 * xe_13) * xu_m1
                U3_4 = (sqrt2_3 * xe_04) * xu_1 
            elseif l == -1
                xu_m1 = SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2dexp, -zkk * st2, zkk * ct2dexp)
                xu_1 = SVector{4,ComplexF64}(zk0m * ct2exp, zk0m * st2, -zkk * ct2exp, -zkk * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2
                U3_1 = (sqrt2_3 * xe_01) * xu_m1 + (sqrt1_3 * xe_m11) * xu_1
                U3_2 = (sqrt2_3 * xe_02) * xu_m1 + (sqrt1_3 * xe_m12) * xu_1
                U3_3 = (sqrt2_3 * xe_03) * xu_m1 + (sqrt1_3 * xe_m13) * xu_1
                U3_4 = (sqrt2_3 * xe_04) * xu_m1 
            elseif l == 3
                xu_1 = SVector{4,ComplexF64}(zk0m * ct2exp, zk0m * st2+0.0im, -zkk * ct2exp, -zkk * st2+0.0im)
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2

                U3_1 = xe_11 * xu_1
                U3_2 = xe_12 * xu_1
                U3_3 = xe_13 * xu_1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            elseif l == -3
                xu_m1 = SVector{4,ComplexF64}(-zk0m * st2, zk0m * ct2dexp, -zkk * st2, zkk * ct2dexp)
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2
                U3_1 = xe_m11 * xu_m1
                U3_2 = xe_m12 * xu_m1
                U3_3 = xe_m13 * xu_m1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            end
        else
            if l == 1
                xu_m1 = SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
                xu_1 = SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2

                U3_1 = (sqrt2_3 * xe_01) * xu_1 + (sqrt1_3 * xe_11) * xu_m1
                U3_2 = (sqrt2_3 * xe_02) * xu_1 + (sqrt1_3 * xe_12) * xu_m1
                U3_3 = (sqrt2_3 * xe_03) * xu_1 + (sqrt1_3 * xe_13) * xu_m1
                U3_4 = (sqrt2_3 * xe_04) * xu_1 
            elseif l == -1
                xu_m1 = SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
                xu_1 = SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2)
                xe_01 = zk0_zm * st * cp
                xe_02 = zk0_zm * st * sp
                xe_03 = zk0_zm * ct
                xe_04 = zkk / zm
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2

                U3_1 = (sqrt2_3 * xe_01) * xu_m1 + (sqrt1_3 * xe_m11) * xu_1
                U3_2 = (sqrt2_3 * xe_02) * xu_m1 + (sqrt1_3 * xe_m12) * xu_1
                U3_3 = (sqrt2_3 * xe_03) * xu_m1 + (sqrt1_3 * xe_m13) * xu_1
                U3_4 = (sqrt2_3 * xe_04) * xu_m1 
            elseif l == 3
                xu_1 = SVector{4,ComplexF64}(zkk * ct2dexp, zkk * st2, zk0m * ct2dexp, zk0m * st2)
                xe_11 = (im * sp - ct * cp) * dexpdsqrt2
                xe_12 = (-im * cp - ct * sp) * dexpdsqrt2
                xe_13 = st * dexpdsqrt2

                U3_1 = xe_11 * xu_1
                U3_2 = xe_12 * xu_1
                U3_3 = xe_13 * xu_1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            elseif l == -3
                xu_m1 = SVector{4,ComplexF64}(zkk * st2, -zkk * ct2exp, -zk0m * st2, zk0m * ct2exp)
                xe_m11 = (im * sp + ct * cp) * expdsqrt2
                xe_m12 = (-im * cp + ct * sp) * expdsqrt2
                xe_m13 = -st * expdsqrt2
                U3_1 = xe_m11 * xu_m1
                U3_2 = xe_m12 * xu_m1
                U3_3 = xe_m13 * xu_m1
                U3_4 = SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im)
            end
        end
        return SVector{5,SVector{4, ComplexF64}}(U3_1, U3_2, U3_3, U3_4, SVector{4,ComplexF64}(0.0im, 0.0im, 0.0im, 0.0im))
        #return (U3_1', U3_2', U3_3', U3_4', SVector(0.0im, 0.0im, 0.0im, 0.0im)')
    end
end

@inline function LC(a::SVector, b::SVector, c::SVector)

    V4 = -a[1] * b[2] * c[3] + a[1] * b[3] * c[2] + a[2] * b[1] * c[3] - a[2] * b[3] * c[1] - a[3] * b[1] * c[2] + a[3] * b[2] * c[1]
    V1 = -a[4] * b[2] * c[3] + a[4] * b[3] * c[2] + a[2] * b[4] * c[3] - a[2] * b[3] * c[4] - a[3] * b[4] * c[2] + a[3] * b[2] * c[4]
    V2 = a[4] * b[1] * c[3] - a[4] * b[3] * c[1] - a[1] * b[4] * c[3] + a[1] * b[3] * c[4] + a[3] * b[4] * c[1] - a[3] * b[1] * c[4]
    V3 = -a[4] * b[1] * c[2] + a[4] * b[2] * c[1] + a[1] * b[4] * c[2] - a[1] * b[2] * c[4] - a[2] * b[4] * c[1] + a[2] * b[1] * c[4]
    return @SVector [V1, V2, V3, V4, 0]
end

@inline function LC(i0::Int64, i1::Int64, i2::Int64, i3::Int64)
    # 如果任意两个索引相等，则返回0
    if i0 == i1 || i0 == i2 || i0 == i3 || i1 == i2 || i1 == i3 || i2 == i3
        return 0.0
    end

    # 用 SVector 存储索引
    iv = @SVector [i0, i1, i2, i3]

    # 计算排列符号
    sign = 1
    for i in 1:3
        for j in i+1:4
            if iv[i] > iv[j]
                sign *= -1
            end
        end
    end
    # here we choose epsilon_{0123}=-1. Other vector should be A^\mu
    return -sign
end

@inline function LC(a::SVector, b::SVector, c::SVector, d::SVector)
    return  a[4]*b[1]*c[2]*d[3] - a[4]*b[1]*c[3]*d[2] +
    a[4]*b[2]*c[3]*d[1] - a[4]*b[2]*c[1]*d[3] +
    a[4]*b[3]*c[1]*d[2] - a[4]*b[3]*c[2]*d[1] +
    a[1]*b[4]*c[3]*d[2] - a[1]*b[4]*c[2]*d[3] +
    a[1]*b[2]*c[4]*d[3] - a[1]*b[2]*c[3]*d[4] +
    a[1]*b[3]*c[2]*d[4] - a[1]*b[3]*c[4]*d[2] +
    a[2]*b[4]*c[1]*d[3] - a[2]*b[4]*c[3]*d[1] +
    a[2]*b[1]*c[3]*d[4] - a[2]*b[1]*c[4]*d[3] +
    a[2]*b[3]*c[4]*d[1] - a[2]*b[3]*c[1]*d[4] +
    a[3]*b[4]*c[2]*d[1] - a[3]*b[4]*c[1]*d[2] +
    a[3]*b[1]*c[4]*d[2] - a[3]*b[1]*c[2]*d[4] +
    a[3]*b[2]*c[1]*d[4] - a[3]*b[2]*c[4]*d[1]
end

import Base: *

function *(Q::SVector{5,Float64}, W::SVector{5,Float64})::Float64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end
function *(Q::SVector{5,Float64}, W::SVector{5,ComplexF64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end
function *(Q::SVector{5,ComplexF64}, W::SVector{5,Float64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end
function *(Q::SVector{5,ComplexF64}, W::SVector{5,ComplexF64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

function *(A::SVector{4, ComplexF64}, M::SMatrix{4, 4, ComplexF64, 16})
    return A'*M  
end
function *(A::SVector{4, ComplexF64}, B::SVector{4, ComplexF64})
    return A'*B  
end

end
