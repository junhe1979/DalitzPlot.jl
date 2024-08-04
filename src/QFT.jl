
module QFT

#export III, GA, GS, epsilon, Uc, Ubc, LCV
using StaticArrays

#############################################################################
# QFT #
#############################################################################
const III = SMatrix{4,4,ComplexF64}([
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

function GS(k::SVector{5,ComplexF64})::SMatrix{4,4,ComplexF64}
    tmp = @MArray zeros(ComplexF64, 4, 4)
    tmp[1, 1] = k[4]
    tmp[1, 3] = -k[3]
    tmp[1, 4] = -k[1] + im * k[2]
    tmp[2, 2] = k[4]
    tmp[2, 3] = -k[1] - im * k[2]
    tmp[2, 4] = k[3]
    tmp[3, 1] = k[3]
    tmp[3, 2] = k[1] - im * k[2]
    tmp[3, 3] = -k[4]
    tmp[4, 1] = k[1] + im * k[2]
    tmp[4, 2] = -k[3]
    tmp[4, 4] = -k[4]
    return tmp
    #return GA[4] * k[4] - GA[1] * k[1] - GA[2] * k[2] - GA[3] * k[3]
end

function GS(k::Vector{ComplexF64})::Matrix{ComplexF64}
    tmp =  zeros(ComplexF64, 4, 4)
    tmp[1, 1] = k[4]
    tmp[1, 3] = -k[3]
    tmp[1, 4] = -k[1] + im * k[2]
    tmp[2, 2] = k[4]
    tmp[2, 3] = -k[1] - im * k[2]
    tmp[2, 4] = k[3]
    tmp[3, 1] = k[3]
    tmp[3, 2] = k[1] - im * k[2]
    tmp[3, 3] = -k[4]
    tmp[4, 1] = k[1] + im * k[2]
    tmp[4, 2] = -k[3]
    tmp[4, 4] = -k[4]
    return tmp
    #return GA[4] * k[4] - GA[1] * k[1] - GA[2] * k[2] - GA[3] * k[3]
end

function epsilon(k::SVector{5,Float64}, ib::Int64, l::Int64)::SVector{5,ComplexF64}
    zm = Float64(k[5])
    zk0 = Float64(k[4])
    zkx = Float64(k[1])
    zky = Float64(k[2])
    zkz = Float64(k[3])
    zkk = sqrt(zkx^2 + zky^2 + zkz^2)

    if zkk < 1e-50
        ct = 1.0
        st = 0.0
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
    elseif abs(zkx) < 1e-50 && abs(zky) < 1e-50
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
    else
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        cp = zkx / zkk / st
        sp = zky / zkk / st
        xexpp = cp + im * sp
    end

    eps = @MArray zeros(Complex{Float64}, 5)
    if ib == 0
        if l == 1
            eps[4] = 0.0
            eps[1] = (im * sp - ct * cp) / (xexpp * sqrt(2.0))
            eps[2] = (-im * cp - ct * sp) / (xexpp * sqrt(2.0))
            eps[3] = st / (xexpp * sqrt(2.0))
        elseif l == 0
            eps[4] = zkk / zm
            eps[1] = zk0 / zm * st * cp
            eps[2] = zk0 / zm * st * sp
            eps[3] = zk0 / zm * ct

        elseif l == -1
            eps[4] = 0.0
            eps[1] = (im * sp + ct * cp) * xexpp / (sqrt(2.0))
            eps[2] = (-im * cp + ct * sp) * xexpp / (sqrt(2.0))
            eps[3] = -st * xexpp / (sqrt(2.0))
        end
    elseif ib == 1
        if l == 1
            eps[4] = 0.0
            eps[1] = (-im * sp - ct * cp) * xexpp / (sqrt(2.0))
            eps[2] = (im * cp - ct * sp) * xexpp / (sqrt(2.0))
            eps[3] = st * xexpp / (sqrt(2.0))
        elseif l == 0
            eps[4] = zkk / zm
            eps[1] = zk0 / zm * st * cp
            eps[2] = zk0 / zm * st * sp
            eps[3] = zk0 / zm * ct
        elseif l == -1
            eps[4] = 0.0
            eps[1] = (-im * sp + ct * cp) / (xexpp * sqrt(2.0))
            eps[2] = (im * cp + ct * sp) / (xexpp * sqrt(2.0))
            eps[3] = -st / (xexpp * sqrt(2.0))
        end
    else
        println("wrong helicity for U")
    end
    return eps
end

function Uc(k::SVector{5,Float64}, l::Int64)::MVector{4,ComplexF64}
    zkx, zky, zkz, zk0, zm = k
    zkk::Float64 = sqrt(zkx^2 + zky^2 + zkz^2)
    zk0m::Float64 = zk0 + zm
    zfac::Float64 = sqrt(2.0 * zm * zk0m)

    if zkk < 1e-50
        ct = 1.0
        st = 0.0
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
        ct2 = 1.0 / zfac
        st2 = 0.0
    elseif abs(zkx) < 1e-50 && abs(zky) < 1e-50
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
    else
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = zkx / zkk / st
        sp = zky / zkk / st
        xexpp = cp + im * sp
    end

    U = @MVector zeros(ComplexF64, 4)

    if l == 1
        U[1] = zk0m * ct2 / xexpp
        U[2] = zk0m * st2
        U[3] = zkk * ct2 / xexpp
        U[4] = zkk * st2
    elseif l == -1
        U[1] = -zk0m * st2
        U[2] = zk0m * ct2 * xexpp
        U[3] = zkk * st2
        U[4] = -zkk * ct2 * xexpp
    else
        println("wrong helicity for U")
    end

    return U
end
function Uc(k::SVector{5,ComplexF64}, l::Int64)::MVector{4,ComplexF64}
    zkx = real(k[1])
    zky = real(k[2])
    zkz = real(k[3])
    zk0 = k[4]
    zm = real(k[5])
    zkk::Float64 = sqrt(zkx^2 + zky^2 + zkz^2)
    zk0m::ComplexF64 = zk0 + zm
    zfac::ComplexF64 = sqrt(2.0 * zm * zk0m)

    if zkk < 1e-30
        ct = 1.0
        st = 0.0
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
        ct2 = 1.0 / zfac
        st2 = 0.0
    elseif abs(zkx) < 1e-30 && abs(zky) < 1e-30
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
    else
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = zkx / zkk / st
        sp = zky / zkk / st
        xexpp = cp + im * sp
    end

    U = @MVector zeros(ComplexF64, 4)

    if l == 1
        U[1] = zk0m * ct2 / xexpp
        U[2] = zk0m * st2
        U[3] = zkk * ct2 / xexpp
        U[4] = zkk * st2
    elseif l == -1
        U[1] = -zk0m * st2
        U[2] = zk0m * ct2 * xexpp
        U[3] = zkk * st2
        U[4] = -zkk * ct2 * xexpp
    else
        println("wrong helicity for U")
    end

    return U
end

function Ubc(k::SVector{5,Float64}, l::Int64)::MVector{4,ComplexF64}
    zkx, zky, zkz, zk0, zm = k
    zkk::Float64 = sqrt(zkx^2 + zky^2 + zkz^2)
    zk0m::Float64 = zk0 + zm
    zfac::Float64 = sqrt(2.0 * zm * zk0m)

    if zkk < 1e-50
        ct = 1.0
        st = 0.0
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
        ct2 = 1.0 / zfac
        st2 = 0.0
    elseif abs(zkx) < 1e-50 && abs(zky) < 1e-50
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
    else
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = zkx / zkk / st
        sp = zky / zkk / st
        xexpp = cp + im * sp
    end

    U = @MVector zeros(ComplexF64, 4)

    if l == 1
        U[1] = zk0m * ct2 * xexpp
        U[2] = zk0m * st2
        U[3] = -zkk * ct2 * xexpp
        U[4] = -zkk * st2
    elseif l == -1
        U[1] = -zk0m * st2
        U[2] = zk0m * ct2 / xexpp
        U[3] = -zkk * st2
        U[4] = zkk * ct2 / xexpp
    else
        println("wrong helicity for U")
    end
    return U
end
function Ubc(k::SVector{5,ComplexF64}, l::Int64)::MVector{4,ComplexF64}
    zkx = real(k[1])
    zky = real(k[2])
    zkz = real(k[3])
    zk0 = k[4]
    zm = real(k[5])
    zkk::Float64 = sqrt(zkx^2 + zky^2 + zkz^2)
    zk0m::ComplexF64 = zk0 + zm
    zfac::ComplexF64 = sqrt(2.0 * zm * zk0m)

    if zkk < 1e-30
        ct = 1.0
        st = 0.0
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
        ct2 = 1.0 / zfac
        st2 = 0.0
    elseif abs(zkx) < 1e-30 && abs(zky) < 1e-30
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = 1.0
        sp = 0.0
        xexpp = 1.0
    else
        ct = zkz / zkk
        st = sqrt(1.0 - ct^2)
        ct2 = sqrt((1.0 + ct) / 2.0) / zfac
        st2 = sqrt((1.0 - ct) / 2.0) / zfac
        cp = zkx / zkk / st
        sp = zky / zkk / st
        xexpp = cp + im * sp
    end

    U = @MVector zeros(ComplexF64, 4)

    if l == 1
        U[1] = zk0m * ct2 * xexpp
        U[2] = zk0m * st2
        U[3] = -zkk * ct2 * xexpp
        U[4] = -zkk * st2
    elseif l == -1
        U[1] = -zk0m * st2
        U[2] = zk0m * ct2 / xexpp
        U[3] = -zkk * st2
        U[4] = zkk * ct2 / xexpp
    else
        println("wrong helicity for Ub")
    end
    return U
end

function LCV(a, b, c)::SVector{5,ComplexF64}
    V = @MArray zeros(ComplexF64, 5)
    V[4] = -a[1] * b[2] * c[3] + a[1] * b[3] * c[2] + a[2] * b[1] * c[3] - a[2] * b[3] * c[1] - a[3] * b[1] * c[2] + a[3] * b[2] * c[1]
    V[1] = -a[4] * b[2] * c[3] + a[4] * b[3] * c[2] + a[2] * b[4] * c[3] - a[2] * b[3] * c[4] - a[3] * b[4] * c[2] + a[3] * b[2] * c[4]
    V[2] = a[4] * b[1] * c[3] - a[4] * b[3] * c[1] - a[1] * b[4] * c[3] + a[1] * b[3] * c[4] + a[3] * b[4] * c[1] - a[3] * b[1] * c[4]
    V[3] = -a[4] * b[1] * c[2] + a[4] * b[2] * c[1] + a[1] * b[4] * c[2] - a[1] * b[2] * c[4] - a[2] * b[4] * c[1] + a[2] * b[1] * c[4]
    return V
end

function cdot(Q::SVector{5,Float64}, W::SVector{5,Float64})::Float64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

function cdot(Q::SVector{5,Float64}, W::SVector{5,ComplexF64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

function cdot(Q::SVector{5,ComplexF64}, W::SVector{5,Float64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end
function cdot(Q::SVector{5,ComplexF64}, W::SVector{5,ComplexF64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

function cdot(Q::Vector{Float64}, W::Vector{Float64})::Float64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

function cdot(Q::Vector{ComplexF64}, W::Vector{ComplexF64})::ComplexF64
    temp = Q[4] * W[4] - Q[1] * W[1] - Q[2] * W[2] - Q[3] * W[3]
    return temp
end

import Base: *

function *(A::MVector{4,ComplexF64}, B::SMatrix{4,4,ComplexF64,16})::MVector{4,ComplexF64}
    temp = @MVector zeros(Complex{Float64}, 4)  # 使用Complex{Float64}来指定元素的类型
    @inbounds for i in 1:4
        @inbounds for j in 1:4
            temp[i] += A[j] * B[j, i]
        end
    end
    return temp
end

function *(A::MVector{4,ComplexF64}, B::MVector{4,ComplexF64})::ComplexF64
    temp = ComplexF64(0.0)  # 使用complex(0.0, 0.0)来创建一个复数
    @inbounds for i in 1:4
        temp += A[i] * B[i]
    end
    return temp
end

function *(A::SMatrix{4,4,ComplexF64,16}, B::SMatrix{4,4,ComplexF64,16})::SMatrix{4,4,ComplexF64,16}
    C = @MArray zeros(Complex{Float64}, 4, 4)   # 使用complex(0.0, 0.0)来创建一个复数
    @inbounds for i in 1:4
        for j in 1:4
            for k in 1:4
                C[i, j] += A[i, k] * B[k, j]
            end
        end
    end
    return C
end


end
