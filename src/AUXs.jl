module AUXs
using Statistics, JLD2, Distributed, Distributions
#############################################################################
# Define a wrapper function that automatically merges fixed and free parameters
# based on the mask, and extracts lower and upper bounds for free parameters
function create_fixed_obj_from_mask(original_obj, lower, upper, full_initial, mask)
    # Indices of free parameters: mask[i]==1 means the parameter is free
    free_indices = [i for i in 1:length(mask) if mask[i] == 1]
    # Indices of fixed parameters: mask[i]==0 means the parameter is fixed
    fixed_indices = [i for i in 1:length(mask) if mask[i] == 0]

    # Construct a wrapped objective function that only takes free parameters
    function fixed_obj(free_params)
        full_params = similar(full_initial)
        # Assign free parameters to their corresponding positions
        for (j, i) in enumerate(free_indices)
            full_params[i] = free_params[j]
        end
        # Use original values from full_initial for fixed parameters
        for i in fixed_indices
            full_params[i] = full_initial[i]
        end
        return original_obj(full_params)
    end

    # Initial guess for free parameters
    free_initial = [full_initial[i] for i in free_indices]
    # Lower and upper bounds for free parameters
    free_lower = [lower[i] for i in free_indices]
    free_upper = [upper[i] for i in free_indices]

    return fixed_obj, free_lower, free_upper, free_initial
end
function get_full_parameters(result, initial, mask)
    # Extract optimal values of free parameters
    free_params_opt = result.minimizer

    # Merge optimal values back into the full parameter vector
    full_params_opt = copy(initial)
    free_indices = [i for i in 1:length(mask) if mask[i] == 1]
    for (j, i) in enumerate(free_indices)
        full_params_opt[i] = free_params_opt[j]
    end

    return full_params_opt
end
#############################################################################
function uncertainty(f, x_opt; Δf=1.0, ε=1e-2, M=5, max_iter=20, tol=0.1)
    n = length(x_opt)
    fmin = mean([f(x_opt) for _ in 1:M])
    unc = zeros(n)

    for i in 1:n
        base = copy(x_opt)
        scale = max(abs(base[i]), 1.0)
        dir = zeros(n)
        dir[i] = 1.0

        # Perform binary search on both sides
        for s in (-1.0, 1.0)
            lo = 0.0
            hi = scale * 1.0
            for _ in 1:max_iter
                mid = (lo + hi) / 2
                x_try = base .+ s * mid * dir
                f_val = median([f(x_try) for _ in 1:M])
                if abs(f_val - fmin - Δf) < tol
                    break
                elseif f_val < fmin + Δf
                    lo = mid
                else
                    hi = mid
                end
            end
            unc[i] += (lo + hi) / 2
        end
        unc[i] /= 2
    end

    return unc
end
function significance(chi2_diff, ndf_diff::Int64)
    # Use high-precision BigFloat type
    chi2_diff = BigFloat(chi2_diff)

    # Create a chi-squared distribution with given degrees of freedom
    chi2_dist = Chisq(ndf_diff)

    # Compute the p-value (tail probability)
    p_value = ccdf(chi2_dist, chi2_diff)

    # Convert p-value to significance (sigma)
    sigma = quantile(Normal(), 1 - p_value)

    return sigma
end
#############################################################################
function broadcast_variable(varname::Symbol, value; filename="temp.jld2", cleanup=true, threshold=10*1024*1024)  # 默认阈值10MB
    # 估算变量大小（近似）
    approx_size = Base.summarysize(value)
    
    if approx_size < threshold
        # 小变量：直接使用 @everywhere
        println("Variable size: $(approx_size) bytes < $(threshold) bytes, using @everywhere")
        timer = time()
        @eval @everywhere const $varname = $value
        elapsed = time() - timer
        println("Broadcast time: $(elapsed)s")
    else
        # 大变量：使用文件分发
        println("Variable size: $(approx_size) bytes >= $(threshold) bytes, using file distribution")
        timer = time()
        JLD2.jldopen(filename, "w") do f
            f[string(varname)] = value
        end
        @sync for pid in workers()
            @async remotecall_wait(pid) do
                data = JLD2.load(filename, string(varname))
                if !isdefined(Main, varname)
                    Core.eval(Main, :(const $(varname) = $data))
                end
            end
        end
        load_time = time() - timer
        println("Broadcast time: $(load_time)s")
        
        cleanup && rm(filename, force=true)
    end
end


end
