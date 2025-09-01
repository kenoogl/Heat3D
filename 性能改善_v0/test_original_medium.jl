#!/usr/bin/env julia

"""
オリジナル版Heat3D中規模テスト (120x120x31)
"""

using Printf
using LinearAlgebra

println("="^50)
println("Heat3D Original Version - Medium Test")
println("="^50)

# オリジナルコードのinclude
include("../q3d/heat3d_nu.jl")

# テストパラメータ
const MODE = 3
const NXY = 120
const NZ = 31
const SOLVER = "pbicgstab"
const SMOOTHER = "gs"
const EPSILON = 1.0e-4

println("Test parameters:")
println("  Grid: $(NXY)x$(NXY)x$(NZ)")
println("  Solver: $(SOLVER)")
println("  Smoother: $(SMOOTHER)")
println("  Tolerance: $(EPSILON)")
println("  Threads: $(Threads.nthreads())")
println("")

# グローバル変数設定
global mode = MODE
global itr_tol = EPSILON

# パラメータ設定
MX = MY = NXY + 2
MZ = NZ + 2
dh = 1.2e-3 / NXY
dh = round(dh, digits=8)
SZ = (MX, MY, MZ)
Δh = (dh, dh, dh)
ox = (0.0, 0.0, 0.0)

println("Executing original version (medium scale)...")
start_time = time()

try
    # 配列初期化
    Z = zeros(Float64, SZ[3])
    ΔZ = zeros(Float64, SZ[3]-1)
    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0
    ID = zeros(UInt8, SZ[1], SZ[2], SZ[3])
    
    # 前処理
    preprocess!(SZ, λ, Z, ΔZ, ox, Δh, ID)
    
    # 初期温度
    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    θ .= Constant.θ_amb
    
    # メイン計算（プロット処理は除外）
    conv_data = main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, SOLVER, SMOOTHER)
    
    execution_time = time() - start_time
    
    # 結果の統計
    z_st = 3
    z_ed = SZ[3] - 1
    s = θ[2:SZ[1]-1, 2:SZ[2]-1, z_st:z_ed]
    min_val = minimum(s)
    max_val = maximum(s)
    l2_norm = norm(s, 2)
    
    println("="^50)
    println("ORIGINAL VERSION RESULTS (MEDIUM)")
    println("="^50)
    println("✅ Execution completed successfully!")
    println("Execution time: $(@sprintf("%.4f", execution_time)) seconds")
    
    iterations = length(conv_data.residuals)
    if iterations > 0
        final_residual = conv_data.residuals[end]
        initial_residual = conv_data.residuals[1]
        
        println("Iterations: $(iterations)")
        println("Initial residual: $(@sprintf("%.6E", initial_residual))")
        println("Final residual: $(@sprintf("%.6E", final_residual))")
        
        # 収束率の計算
        if initial_residual > 0
            reduction_factor = final_residual / initial_residual
            order_reduction = log10(initial_residual / final_residual)
            
            println("Residual reduction factor: $(@sprintf("%.6E", reduction_factor))")
            println("Order reduction: $(@sprintf("%.2f", order_reduction))")
        end
    end
    
    println("")
    println("Temperature Field Statistics:")
    println("  Min temperature: $(@sprintf("%.6f", min_val)) K")
    println("  Max temperature: $(@sprintf("%.6f", max_val)) K")
    println("  L2 norm: $(@sprintf("%.6E", l2_norm))")
    
    # 結果保存
    global original_result_medium = Dict(
        "execution_time" => execution_time,
        "iterations" => iterations,
        "initial_residual" => initial_residual,
        "final_residual" => final_residual,
        "min_temp" => min_val,
        "max_temp" => max_val,
        "l2_norm" => l2_norm,
        "grid_size" => "$(NXY)x$(NXY)x$(NZ)"
    )
    
catch e
    execution_time = time() - start_time
    println("❌ Error occurred after $(@sprintf("%.4f", execution_time)) seconds")
    println("Error: $(e)")
    println("Stack trace:")
    for (i, frame) in enumerate(stacktrace())
        println("  $i: $frame")
        if i >= 10
            break
        end
    end
end

println("\nOriginal version medium test completed.")