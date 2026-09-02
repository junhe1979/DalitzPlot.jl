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

    # 预计算斜率
    slopes = zeros(length(x_sorted)-1)
    for i in eachindex(slopes)
        slopes[i] = (y_sorted[i+1] - y_sorted[i]) / (x_sorted[i+1] - x_sorted[i])
    end

    # 插值函数（保持不变）
    function interp(xx)
        if xx <= x_sorted[1]
            return y_sorted[1] + slopes[1] * (xx - x_sorted[1])
        elseif xx >= x_sorted[end]
            return y_sorted[end] + slopes[end] * (xx - x_sorted[end])
        else
            lo, hi = 1, length(x_sorted)
            while hi - lo > 1
                mid = (lo + hi) ÷ 2
                if x_sorted[mid] <= xx
                    lo = mid
                else
                    hi = mid
                end
            end
            return y_sorted[lo] + slopes[lo] * (xx - x_sorted[lo])
        end
    end

    # ========== 主要修改部分：对每个 bin 积分正值 ==========
    for i in 1:nbins
        center = x_data_local[i]
        a = center - bin_width / 2
        b = center + bin_width / 2

        # Simpson 积分参数
        n = npts * 2
        h = (b - a) / n

        # 端点值：只取正值部分
        f_a = max(interp(a), 0.0)
        f_b = max(interp(b), 0.0)

        sum_odd = 0.0
        sum_even = 0.0

        for j in 1:(n-1)
            x_val = a + j * h
            # 核心改动：插值后若为负，直接置 0 再累加
            f_val = max(interp(x_val), 0.0)
            if j % 2 == 1
                sum_odd += f_val
            else
                sum_even += f_val
            end
        end

        # Simpson 积分公式
        integral = h / 3 * (f_a + f_b + 4 * sum_odd + 2 * sum_even)

        # 除以 bin 宽度得到平均值（正值部分的平均高度）
        y_th[i] = integral / bin_width
    end

    return y_th
end

#############################################################################

function estimate_errors_symmetry(p0_best, chi2f, mask, lower, upper;
                                step_fraction=0.02, verbose=true)
    """
    更健壮的版本：使用相对步长，自动适应参数尺度
    """
    n_params = length(p0_best)
    errors = zeros(n_params)
    chi2_best = chi2f(p0_best)
    total_calls = 1
    
    println("="^60)
    println("📊 健壮误差估计")
    println("="^60)
    
    for i in 1:n_params
        if mask[i] == 0
            continue
        end
        
        # 使用相对步长
        param_scale = abs(p0_best[i])
        if param_scale < 0.01
            # 非常小的参数，使用绝对步长
            step = 1e-4
        elseif param_scale < 0.1
            step = 0.001
        else
            step = step_fraction * param_scale
        end
        
        # 测试正负方向
        deltas = Vector{Float64}()
        steps_tested = Vector{Float64}()
        
        for sign in [-1.0, 1.0]
            p_test = copy(p0_best)
            p_test[i] += sign * step
            if p_test[i] >= lower[i] && p_test[i] <= upper[i]
                chi2_val = chi2f(p_test)
                total_calls += 1
                push!(deltas, chi2_val - chi2_best)
                push!(steps_tested, step)
            end
        end
        
        # 如果 Δχ² 太小，增大步长
        if !isempty(deltas) && all(d -> d < 0.1, deltas)
            step *= 3.0
            deltas = Vector{Float64}()
            steps_tested = Vector{Float64}()
            for sign in [-1.0, 1.0]
                p_test = copy(p0_best)
                p_test[i] += sign * step
                if p_test[i] >= lower[i] && p_test[i] <= upper[i]
                    chi2_val = chi2f(p_test)
                    total_calls += 1
                    push!(deltas, chi2_val - chi2_best)
                    push!(steps_tested, step)
                end
            end
        end
        
        # 计算误差
        if length(deltas) >= 2
            avg_delta = mean(deltas)
            if avg_delta > 1e-6
                # 用平均 Δχ² 估计误差
                errors[i] = step / sqrt(avg_delta)
            else
                errors[i] = step * 5.0  # 保守估计
            end
        elseif length(deltas) == 1 && deltas[1] > 0
            errors[i] = step / sqrt(deltas[1])
        else
            errors[i] = 0.1 * max(abs(p0_best[i]), 1.0)
        end
        
        if verbose
            @printf("参数 %2d: %10.6f ± %8.6f  (步长=%.4f)\n", 
                   i, p0_best[i], errors[i], step)
        end
    end
    
    println("="^60)
    println("✅ 健壮误差估计完成")
    println("  总函数调用: $total_calls 次")
    println("="^60)
    
    return errors
end

function quadratic_fit(x, y)
    """
    简单的二次拟合: y = a*x^2 + b*x + c
    返回 [a, b, c]
    """
    n = length(x)
    if n < 3
        return [0.0, 0.0, 0.0]
    end

    # 构建矩阵 A (n x 3)
    A = zeros(n, 3)
    for i in 1:n
        A[i, 1] = x[i]^2
        A[i, 2] = x[i]
        A[i, 3] = 1.0
    end

    # 最小二乘解: (A'*A) \ (A'*y)
    try
        coeffs = (A' * A) \ (A' * y)
        return coeffs
    catch
        return [0.0, 0.0, 0.0]
    end
end

function estimate_errors(p0_best, chi2f, mask, lower, upper;
                                        precision_level=3,  # 1=快速, 2=标准, 3=高精度
                                        verbose=true)
    """
    高精度误差估计：多层自适应扫描
    
    precision_level:
        1: 快速 (每个参数 2-3 个点)
        2: 标准 (每个参数 4-6 个点) 
        3: 高精度 (每个参数 8-12 个点)
    """
    n_params = length(p0_best)
    errors = Vector{Tuple{Float64, Float64}}(undef, n_params)
    
    # 设置精度参数
    precision_configs = Dict(
        1 => (n_points=3, max_step=0.05, min_step=0.005),
        2 => (n_points=6, max_step=0.08, min_step=0.002),
        3 => (n_points=12, max_step=0.12, min_step=0.001)
    )
    config = precision_configs[precision_level]
    n_points = config.n_points
    max_step = config.max_step
    min_step = config.min_step
    
    println("="^70)
    println("🎯 高精度误差估计 (精度级别: $precision_level)")
    println("="^70)
    t_start = time()
    
    chi2_best = chi2f(p0_best)
    total_calls = 1
    println("χ²_best = $chi2_best")
    println("-"^70)
    
    # ===== 第一阶段：粗扫评估重要性 =====
    println("阶段1: 评估参数重要性...")
    importance = Dict{Int, Float64}()
    
    for i in 1:n_params
        if mask[i] == 0
            continue
        end
        param_scale = max(abs(p0_best[i]), 1e-3)
        step = 0.02 * param_scale
        
        # 只测一个方向
        p_test = copy(p0_best)
        p_test[i] += step
        if p_test[i] <= upper[i]
            chi2_val = chi2f(p_test)
            total_calls += 1
            delta = max(chi2_val - chi2_best, 1e-10)
            importance[i] = delta
        else
            importance[i] = 0.0
        end
    end
    
    # 分类参数
    sorted_importance = sort([(k, v) for (k, v) in importance], by=x->x[2], rev=true)
    n_important = max(1, div(length(sorted_importance), 3))
    n_medium = max(1, div(length(sorted_importance), 3))
    
    important_params = [x[1] for x in sorted_importance[1:n_important]]
    medium_params = [x[1] for x in sorted_importance[n_important+1:n_important+n_medium]]
    weak_params = [x[1] for x in sorted_importance[n_important+n_medium+1:end]]
    
    if verbose
        println("  重要参数 ($(length(important_params))个): ", important_params)
        println("  中等参数 ($(length(medium_params))个): ", medium_params)
        println("  弱参数 ($(length(weak_params))个): ", weak_params)
    end
    println("-"^70)
    
    # ===== 第二阶段：分层精细扫描 =====
    println("阶段2: 分层精细扫描...")
    
    for i in 1:n_params
        if mask[i] == 0
            errors[i] = (0.0, 0.0)
            continue
        end
        
        param_scale = max(abs(p0_best[i]), 1e-3)
        
        # 根据重要性决定采样点数
        if i in important_params
            n_scan = n_points
            max_step_local = max_step
            min_step_local = min_step
            label = "重要"
        elseif i in medium_params
            n_scan = max(div(n_points, 2), 3)
            max_step_local = max_step * 0.8
            min_step_local = min_step * 2
            label = "中等"
        else
            n_scan = max(div(n_points, 3), 2)
            max_step_local = max_step * 0.6
            min_step_local = min_step * 3
            label = "弱"
        end
        
        # 精细扫描
        err_down, err_up, calls = scan_parameter_high_precision(
            p0_best, i, chi2f, chi2_best, lower, upper,
            param_scale, n_scan, max_step_local, min_step_local
        )
        
        total_calls += calls
        errors[i] = (err_down, err_up)
        
        if verbose
            @printf("参数 %2d [%s]: %10.6f  -%8.6f / +%8.6f  (调用 %d次)\n", 
                   i, label, p0_best[i], err_down, err_up, calls)
        end
    end
    
    elapsed = time() - t_start
    println("-"^70)
    println("✅ 高精度误差估计完成")
    println("  总函数调用: $total_calls 次")
    @printf("  总耗时: %.1f 秒 (%.1f 分钟)\n", elapsed, elapsed/60)
    @printf("  平均每次: %.1f 秒\n", elapsed/total_calls)
    println("="^70)
    
    return errors
end

function scan_parameter_high_precision(p0_best, idx, chi2f, chi2_best, lower, upper,
                                       param_scale, n_points, max_step, min_step)
    """
    高精度单参数扫描：使用自适应非均匀采样
    """
    # ===== 生成自适应采样点 =====
    # 在中心附近更密集，远离中心更稀疏
    points_down = Float64[]
    points_up = Float64[]
    
    # 对数分布：中心密集，边缘稀疏
    for j in 1:n_points
        frac = (j - 1) / (n_points - 1)
        # 使用平方分布使中心更密集
        t = frac^2
        step = min_step + t * (max_step - min_step)
        
        push!(points_down, -step * param_scale)
        push!(points_up, step * param_scale)
    end
    
    # 收集数据
    down_data = Vector{Tuple{Float64, Float64}}()  # (step, delta)
    up_data = Vector{Tuple{Float64, Float64}}()
    calls = 0
    
    # 扫描所有点
    for (step_down, step_up) in zip(points_down, points_up)
        # 向下
        if step_down != 0
            p_test = copy(p0_best)
            p_test[idx] += step_down
            if p_test[idx] >= lower[idx]
                chi2_val = chi2f(p_test)
                calls += 1
                delta = max(chi2_val - chi2_best, 1e-10)
                push!(down_data, (abs(step_down), delta))
            end
        end
        
        # 向上
        if step_up != 0
            p_test = copy(p0_best)
            p_test[idx] += step_up
            if p_test[idx] <= upper[idx]
                chi2_val = chi2f(p_test)
                calls += 1
                delta = max(chi2_val - chi2_best, 1e-10)
                push!(up_data, (abs(step_up), delta))
            end
        end
    end
    
    # ===== 高精度插值找 Δχ² = 1 =====
    err_down = interpolate_high_precision(down_data, 1.0)
    err_up = interpolate_high_precision(up_data, 1.0)
    
    # 如果插值失败，使用拟合
    if err_down <= 0 && length(down_data) >= 3
        err_down = fit_quadratic_error(down_data, 1.0)
    end
    if err_up <= 0 && length(up_data) >= 3
        err_up = fit_quadratic_error(up_data, 1.0)
    end
    
    # 保底估计
    if err_down <= 0
        err_down = 0.05 * param_scale
    end
    if err_up <= 0
        err_up = 0.05 * param_scale
    end
    
    return err_down, err_up, calls
end

function interpolate_high_precision(data, target)
    """
    高精度插值：使用三次样条或分段插值
    """
    if length(data) < 2
        return 0.0
    end
    
    # 按步长排序
    sort!(data, by=x->x[1])
    x = [d[1] for d in data]
    y = [d[2] for d in data]
    
    # 如果第一个点就超过目标
    if y[1] > target
        if length(x) >= 3
            # 用前三个点外推
            return extrapolate_quadratic(x[1:3], y[1:3], target)
        else
            return x[1] * 0.8
        end
    end
    
    # 如果最后一个点还小于目标
    if y[end] < target
        if length(x) >= 3
            # 用后三个点外推
            return extrapolate_quadratic(x[end-2:end], y[end-2:end], target)
        else
            return x[end] * 1.5
        end
    end
    
    # 正常插值
    for j in 1:length(x)-1
        if (y[j] <= target && y[j+1] >= target) || 
           (y[j] >= target && y[j+1] <= target)
            # 使用周围的4个点进行三次插值
            start = max(1, j-1)
            stop = min(length(x), j+3)
            if stop - start + 1 >= 3
                return interpolate_cubic(x[start:stop], y[start:stop], target)
            else
                # 线性插值
                if y[j+1] != y[j]
                    t = (target - y[j]) / (y[j+1] - y[j])
                    return x[j] + t * (x[j+1] - x[j])
                end
            end
        end
    end
    
    return 0.0
end

function interpolate_cubic(x, y, target)
    """
    三次插值求目标值
    """
    # 简单实现：用二次拟合 + 求根
    if length(x) >= 3
        # 二次拟合
        A = zeros(length(x), 3)
        for i in eachindex(x)
            A[i, 1] = x[i]^2
            A[i, 2] = x[i]
            A[i, 3] = 1.0
        end
        coeffs = (A' * A) \ (A' * y)
        
        # 解方程: a*x^2 + b*x + c = target
        a, b, c = coeffs[1], coeffs[2], coeffs[3] - target
        disc = b^2 - 4*a*c
        if disc >= 0 && a != 0
            roots = [(-b + sqrt(disc)) / (2*a), (-b - sqrt(disc)) / (2*a)]
            # 选择在数据范围内的根
            for root in roots
                if root >= minimum(x) && root <= maximum(x)
                    return root
                end
            end
            # 如果都不在范围内，选择接近目标的根
            return roots[argmin(abs.(roots .- mean(x)))]
        end
    end
    return 0.0
end

function extrapolate_quadratic(x, y, target)
    """
    用二次拟合外推
    """
    if length(x) >= 3
        A = zeros(length(x), 3)
        for i in eachindex(x)
            A[i, 1] = x[i]^2
            A[i, 2] = x[i]
            A[i, 3] = 1.0
        end
        coeffs = (A' * A) \ (A' * y)
        a, b, c = coeffs[1], coeffs[2], coeffs[3] - target
        disc = b^2 - 4*a*c
        if disc >= 0 && a != 0
            roots = [(-b + sqrt(disc)) / (2*a), (-b - sqrt(disc)) / (2*a)]
            # 选择靠近数据范围的根
            root_avg = mean(x)
            return roots[argmin(abs.(roots .- root_avg))]
        end
    end
    return x[end] * 1.5
end

function fit_quadratic_error(data, target)
    """
    用二次拟合估计误差
    """
    if length(data) < 3
        return 0.0
    end
    
    x = [d[1] for d in data]
    y = [d[2] for d in data]
    
    # 二次拟合
    A = zeros(length(x), 3)
    for i in eachindex(x)
        A[i, 1] = x[i]^2
        A[i, 2] = x[i]
        A[i, 3] = 1.0
    end
    coeffs = (A' * A) \ (A' * y)
    
    a, b, c = coeffs[1], coeffs[2], coeffs[3] - target
    disc = b^2 - 4*a*c
    if disc >= 0 && a != 0
        roots = [(-b + sqrt(disc)) / (2*a), (-b - sqrt(disc)) / (2*a)]
        # 选择正根
        positive_roots = filter(r -> r > 0, roots)
        if !isempty(positive_roots)
            return positive_roots[1]
        end
    end
    
    return 0.0
end


end
