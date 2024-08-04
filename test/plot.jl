using Plots,PythonCall
function plot() 
plt =pyimport("matplotlib.pyplot")
mat =pyimport("matplotlib")
gr()
ENV["GKSwstype"] = "100"
mat.use("Agg")


# 初始化三个空数组
x_values = Float64[]
y_values = Float64[]
z_values = Float64[]

# 读取文件内容并存储到数组中
filename = "outputjl.txt"
open(filename, "r") do file
    for line in eachline(file)
        nums = split(line, ",")
        push!(x_values, parse(Float64, nums[1]))
        push!(y_values, parse(Float64, nums[2]))
        push!(z_values, parse(Float64, nums[3]))
    end
end

# 打印数组内容


# 找到所有独特的 x 和 y 值
# 找到所有独特的 x 和 y 值，并排序
x_unique = sort(unique(x_values))
y_unique = sort(unique(y_values))


# 创建一个矩阵来存储 z 值
z_matrix = Matrix{Float64}(undef, length(y_unique), length(x_unique))

# 填充 z 矩阵
for i in eachindex(x_values)
    x_index = findfirst(isequal(x_values[i]), x_unique)
    y_index = findfirst(isequal(y_values[i]), y_unique)
    z_matrix[y_index, x_index] = z_values[i]
end

# 创建热力图
plt.imshow(z_matrix, extent=[minimum(x_unique), maximum(x_unique),minimum(y_unique),maximum(y_unique)], origin="lower", aspect="auto")
plt.colorbar()
plt.tick_params(axis="both", direction="in")
# 保存热力图为文件
plt.savefig("heatmap.pdf", dpi=100)
end 
plot()