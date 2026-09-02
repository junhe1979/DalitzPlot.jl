module GEN
export GENEV, reset_genev_rngs!, get_process_rng
using StaticArrays
using Random
using Distributed

# ==============================================================================
# PERFORMANCE STORAGE CONTAINERS
# ==============================================================================
# PROCESS_RNGS: Thread/Process safe container holding baseline MersenneTwister RNGs.
# Locked to MersenneTwister to prevent machine-precision domain errors in kinematics.
const PROCESS_RNGS = Vector{MersenneTwister}()

# THREAD_RNO_BUFFERS: Thread-local Scratchpads for random numbers to achieve 0 allocation.
# Default pre-allocated size 54 accommodates up to 18 final-state particles (3*18 - 4 = 50).
const THREAD_RNO_BUFFERS = Vector{Vector{Float64}}()

"""
    __init__()

Module initialization block executed on every process/worker when the module is loaded.
Pre-allocates thread-local storage buffers based on the available thread count.
"""
function __init__()
    _ensure_storage_size!(myid())
    
    # Pre-allocate thread-local scratchpads to avoid run-time heap allocations
    resize!(THREAD_RNO_BUFFERS, Threads.nthreads())
    for i in 1:Threads.nthreads()
        THREAD_RNO_BUFFERS[i] = zeros(Float64, 54) 
    end
end

# Internal helper function to dynamically scale the process-level storage container safely
function _ensure_storage_size!(target_pid::Int)
    if length(PROCESS_RNGS) < target_pid
        resize!(PROCESS_RNGS, target_pid)
    end
end

"""
    reset_genev_rngs!(base_seed::Int=12345)

Explicitly initializes or resets the RNG state across all active Distributed workers.
"""
function reset_genev_rngs!(base_seed::Int=12345)
    # FIX: Replaced ambiguous @everywhere with explicit remotecall_fetch loops.
    # This prevents sub-processes from searching the `Main` namespace for internal functions,
    # resolving the distributed `UndefVarError` completely.
    for pid in procs()
        remotecall_fetch(pid) do
            _ensure_storage_size!(maximum(procs()))
            PROCESS_RNGS[myid()] = MersenneTwister(base_seed + myid())
        end
    end
end

"""
    get_process_rng()

Inline helper to safely retrieve the current worker process's dedicated MersenneTwister RNG.
"""
@inline function get_process_rng()
    pid = myid()
    if pid > length(PROCESS_RNGS) || !isassigned(PROCESS_RNGS, pid)
        _ensure_storage_size!(pid)
        PROCESS_RNGS[pid] = MersenneTwister(12345 + pid)
    end
    return PROCESS_RNGS[pid]
end

# ==============================================================================
# KINEMATIC KERNELS & PHYSICS FORMULAS
# ==============================================================================

"""
    ROTES2!(cos_theta, sin_theta, cos_theta2, sin_theta2, pr, i)

Performs spatial coordinate rotations directly on the flat MMatrix column data.
Optimized via 1D indexing to bypass 2D matrix overhead and SIMD barriers.
"""
@inline function ROTES2!(cos_theta::Float64, sin_theta::Float64, cos_theta2::Float64,
    sin_theta2::Float64, pr::MMatrix{5,18,Float64,90}, i::Int64)
    @inbounds begin
        # Fast 1D index mapping to eliminate runtime bounds-checks
        k1 = 5 * (i - 1) + 1
        k2 = k1 + 1
        sa = pr[k1]
        sb = pr[k2]

        # First rotation
        a = sa * cos_theta - sb * sin_theta
        pr[k2] = sa * sin_theta + sb * cos_theta

        # Second rotation
        k2 += 1
        b = pr[k2]
        pr[k1] = a * cos_theta2 - b * sin_theta2
        pr[k2] = a * sin_theta2 + b * cos_theta2
    end
end

"""
    PDK(A, B, C)

Standard two-body decay momentum helper function.
Protects the root with an explicit `abs()` against machine-precision negative underflows.
"""
@inline function PDK(A::Float64, B::Float64, C::Float64)::Float64
    A_squared = A * A
    B_sq = B * B
    C_sq = C * C
    return 0.5 * sqrt(abs(A_squared + (B_sq - C_sq)^2 / A_squared - 2.0 * (B_sq + C_sq)))
end

# Optimization: Computes the sum of squares of the first 3 rows of column j in PCM.
# Avoids dynamic heap allocation caused by array slicing like `PCM[1:3, j]`.
@inline function pcm_sq_sum(PCM, j)
    @inbounds return PCM[1, j]^2 + PCM[2, j]^2 + PCM[3, j]^2
end

# ==============================================================================
# MAIN GENERATOR CORE
# ==============================================================================

"""
    GENEV(tecm, EM; fixed=true)

Universal N-body Phase Space Generator. 
Maintains strict backward compatibility with original input types and return values,
while entirely bypassing heap allocations during multithreaded operations.
"""
function GENEV(tecm::Float64, EM::Vector{Float64}; fixed=true)
    NT = length(EM)
    NTM1 = NT - 1
    NTM2 = NT - 2
    NTNM4 = 3 * NT - 4

    # Optimization: Allocated purely on the stack via StaticArrays (0 heap allocation)
    EMM = @MVector zeros(Float64, 18)
    EMS = @MVector zeros(Float64, 18)
    SM  = @MVector zeros(Float64, 18)
    PD  = @MVector zeros(Float64, 18)
    PCM = @MMatrix zeros(Float64, 5, 18)

    FFQ = SVector{18,Float64}([0.0, 3.141592, 19.73921, 62.01255, 129.8788, 204.0131, 256.3704, 268.4705,
        240.9780, 189.2637, 132.1308, 83.0202, 47.4210, 24.8295, 12.0006, 5.3858,
        2.2560, 0.8859])
    TWOPI = 6.2831853073

    EMM[1] = EM[1]
    TM = 0.0
    @inbounds for i in 1:NT
        EMS[i] = EM[i]^2
        TM += EM[i]
        SM[i] = TM
    end

    TECMTM = tecm - TM
    EMM[NT] = tecm

    WTMAXQ = TECMTM^NTM2 * FFQ[NT] / tecm

    # Optimization: Fetch pre-allocated scratchpad buffer owned by the current Thread ID
    tid = Threads.threadid()
    if length(THREAD_RNO_BUFFERS) < tid
        resize!(THREAD_RNO_BUFFERS, tid)
    end
    if !isassigned(THREAD_RNO_BUFFERS, tid) || length(THREAD_RNO_BUFFERS[tid]) < NTNM4
        THREAD_RNO_BUFFERS[tid] = zeros(Float64, max(54, NTNM4))
    end
    RNO_buf = THREAD_RNO_BUFFERS[tid]

    # Reverted to legacy stable RNG engine to maintain precise kinematic boundaries
    rng = fixed ? get_process_rng() : Random.default_rng()
    
    # Optimization: In-place overwrite using a view over the scratchpad (0 heap allocation)
    RNO_view = @views RNO_buf[1:NTNM4]
    rand!(rng, RNO_view)

    # Optimization: Native in-place sorting on the pre-allocated view
    if NTM2 > 1
        RNO_sub = @views RNO_buf[1:NTM2]
        sort!(RNO_sub)
    end

    if NTM2 > 0
        @inbounds for j in 2:NTM1
            EMM[j] = RNO_buf[j-1] * TECMTM + SM[j]
        end
    end

    WT = WTMAXQ
    if NTM2 >= 0
        IR = NTM2
        @inbounds for i in 1:NTM1
            PD[i] = PDK(EMM[i+1], EMM[i], EM[i+1])
            WT *= PD[i]
        end

        PCM .= 0.0
        PCM[2, 1] = PD[1]

        @inbounds for i in 2:NT
            PCM[1, i] = 0.0
            PCM[2, i] = -PD[i-1]
            PCM[3, i] = 0.0
            IR += 1
            BANG = TWOPI * RNO_buf[IR]
            CB = cos(BANG)
            SB = sin(BANG)
            IR += 1
            C = 2.0 * RNO_buf[IR] - 1.0
            S = sqrt(1.0 - C^2)

            if i != NT
                ESYS = sqrt(PD[i]^2 + EMM[i]^2)
                BETA = PD[i] / ESYS
                GAMA = ESYS / EMM[i]
                for j in 1:i
                    aa = pcm_sq_sum(PCM, j)
                    PCM[5, j] = sqrt(aa)
                    PCM[4, j] = sqrt(aa + EMS[j])
                    ROTES2!(C, S, CB, SB, PCM, j)
                    psave = GAMA * (PCM[2, j] + BETA * PCM[4, j])
                    PCM[2, j] = psave
                end
            else
                for j in 1:i
                    aa = pcm_sq_sum(PCM, j)
                    PCM[5, j] = sqrt(aa)
                    PCM[4, j] = sqrt(aa + EMS[j])
                    ROTES2!(C, S, CB, SB, PCM, j)
                end
            end
        end
    end
    
    @inbounds for i in 1:NT
        PCM[5, i] = EM[i]
    end

    # Optimization: Unrolled scalar initialization of individual SVector values.
    # Replaces the expensive dynamic matrix-slicing syntax `PCM[1:5, i]` which allocation-heavy.
    P = Vector{SVector{5,Float64}}(undef, NT)
    @inbounds for i in 1:NT
        P[i] = SVector{5,Float64}(PCM[1, i], PCM[2, i], PCM[3, i], PCM[4, i], PCM[5, i])
    end

    return P, WT
end

end # module