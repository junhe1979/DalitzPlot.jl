module AUXusing Statistics
# 定义一个包装函数，自动根据掩码将固定参数和自由参数合并，
# 并提取自由参数的下界和上界
function create_fixed_obj_from_mask(original_obj, lower, upper, full_initial, mask)
    # 自由参数的索引：mask[i]==1表示该参数自由
    free_indices = [i for i in 1:length(mask) if mask[i] == 1]
    # 固定参数的索引：mask[i]==0表示固定参数
    fixed_indices = [i for i in 1:length(mask) if mask[i] == 0]

    # 构造包装后的目标函数，只接受自由参数
    function fixed_obj(free_params)
        full_params = similar(full_initial)
        # 将自由参数赋值到对应位置
        for (j, i) in enumerate(free_indices)
            full_params[i] = free_params[j]
        end
        # 固定参数则直接使用 full_initial 中的值
        for i in fixed_indices
            full_params[i] = full_initial[i]
        end
        return original_obj(full_params)
    end

    # 自由参数的初始猜测
    free_initial = [full_initial[i] for i in free_indices]
    # 自由参数对应的下界和上界
    free_lower = [lower[i] for i in free_indices]
    free_upper = [upper[i] for i in free_indices]

    return fixed_obj, free_lower, free_upper, free_initial
end

function get_full_parameters(result, initial, mask)
    # 获取自由参数的最优值
    free_params_opt = result.minimizer

    # 将自由参数的最优值合并回完整的参数向量
    full_params_opt = copy(initial)
    free_indices = [i for i in 1:length(mask) if mask[i] == 1]
    for (j, i) in enumerate(free_indices)
        full_params_opt[i] = free_params_opt[j]
    end

    return full_params_opt

end

function uncertainty(f, x_opt; Δf=1.0, ε=1e-2, M=5, max_iter=20, tol=0.1)
    n = length(x_opt)
    fmin = mean([f(x_opt) for _ in 1:M])
    unc = zeros(n)

    for i in 1:n
        base = copy(x_opt)
        scale = max(abs(base[i]), 1.0)
        dir = zeros(n)
        dir[i] = 1.0

        # 左右各进行一次二分搜索
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



end