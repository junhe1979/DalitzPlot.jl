module AUXs
using Statistics, JLD2, Distributed, Distributions, Printf, Dates
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
function broadcast_variable(varname::Symbol, value; filename="temp.jld2", cleanup=true, threshold=10 * 1024 * 1024)  # 默认阈值10MB
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

macro broadcast(expr)
    quote
        # 显示开始信息
        printstyled("⏳ Broadcasting... "; color=:blue)
        t_start = time_ns()

        # 执行传入的广播表达式
        $(esc(expr))

        # 计算并显示结束信息
        t_end = time_ns()
        elapsed_sec = (t_end - t_start) / 1e9
        @printf("✅ Broadcast successfully! ⏱️ Total elapsed time: %.4f seconds\n", elapsed_sec)
    end
end

macro run(ex)
    quote
        open("log_chi2.txt", "w") do io
            println(io, "#iloop chi2")
        end
        open("log_results.txt", "w") do io
        end
        # 开始信息
        println("╔════════════════════════════════════╗")
        println("║        PROGRAM BEGINNING           ║")
        println("╚════════════════════════════════════╝")
        printstyled("🚀 Started at: ", now(); color=:cyan)
        println(" \n")
        println(repeat('-', 90))

        # 记录开始时间
        t_start = time_ns()

        # 运行传入的表达式
        result = $(esc(ex))

        # 计算执行时间
        t_end = time_ns()
        elapsed_sec = (t_end - t_start) / 1e9

        # 结束信息
        printstyled("🎉 Ended at: ", now(); color=:cyan)
        println()
        @printf("⏱️ Total execution time: %.4f seconds\n", elapsed_sec)
        println(repeat('-', 90))

        # 返回表达式的执行结果
        result
    end
end

function extract_parameters(parameter)
    # 提取初始值
    initial = [p[1] for p in parameter]
    
    # 提取 mask（第4个元素）
    mask = [p[4] for p in parameter]
    
    # 提取 upper（第3个元素）
    upper = [p[3] for p in parameter]
    
    # 提取 lower（第2个元素）
    lower = [p[2] for p in parameter]
    
    return initial, upper, lower, mask
end
function bin_average(x, y, x_data; npts=100)
    # 深拷贝输入
    x_local = copy(x)
    y_local = copy(y)
    x_data_local = copy(x_data)
    
    # 排序
    p = sortperm(x_local)
    x_sorted = x_local[p]
    y_sorted = y_local[p]
    
    nbins = length(x_data_local)
    y_th = zeros(nbins)
    
    if nbins >= 2
        bin_width = x_data_local[2] - x_data_local[1]
    else
        error("至少需要两个 bin 中心")
    end
    
    # 预计算斜率（避免重复计算）
    slopes = zeros(length(x_sorted)-1)
    for i in eachindex(slopes)
        slopes[i] = (y_sorted[i+1] - y_sorted[i]) / (x_sorted[i+1] - x_sorted[i])
    end
    
    # 手动插值函数（纯函数，无状态）
    function interp(xx)
        if xx <= x_sorted[1]
            # 左边界外推
            return y_sorted[1] + slopes[1] * (xx - x_sorted[1])
        elseif xx >= x_sorted[end]
            # 右边界外推
            return y_sorted[end] + slopes[end] * (xx - x_sorted[end])
        else
            # 二分查找区间
            lo, hi = 1, length(x_sorted)
            while hi - lo > 1
                mid = (lo + hi) ÷ 2
                if x_sorted[mid] <= xx
                    lo = mid
                else
                    hi = mid
                end
            end
            # 线性插值
            return y_sorted[lo] + slopes[lo] * (xx - x_sorted[lo])
        end
    end
    
    for i in 1:nbins
        center = x_data_local[i]
        a = center - bin_width / 2
        b = center + bin_width / 2
        
        # 使用 Simpson 法则（更高精度）
        n = npts * 2
        h = (b - a) / n
        
        sum_odd = 0.0
        sum_even = 0.0
        
        for j in 1:n-1
            x_val = a + j * h
            f_val = interp(x_val)
            if j % 2 == 1
                sum_odd += f_val
            else
                sum_even += f_val
            end
        end
        
        integral = h/3 * (interp(a) + interp(b) + 4*sum_odd + 2*sum_even)
        y_th[i] = integral / bin_width
    end
    
    return y_th
end
end
