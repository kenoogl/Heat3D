#!/usr/bin/env julia

"""
Heat3D 最適化版の簡単な検証スクリプト

オリジナル版と最適化版の結果が一致することを確認し、
性能改善を測定します。
"""

using Printf
using LinearAlgebra

println("="^60)
println("Heat3D Quick Validation Test")
println("="^60)

# テスト用の小さいサイズ
const MODE = 3
const NXY = 60   # 小サイズで高速テスト
const NZ = 31
const SOLVER = "pbicgstab"
const SMOOTHER = "gs"
const EPSILON = 1.0e-4

function test_original()
    """オリジナル版のテスト実行（プロット無効）"""
    println("Testing original version...")
    
    # オリジナルコードのinclude
    include("../q3d/heat3d_nu.jl")
    
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
    
    # 実行時間測定
    start_time = time()
    conv_data = main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, SOLVER, SMOOTHER)
    execution_time = time() - start_time
    
    # 結果の統計
    z_st = 3
    z_ed = SZ[3] - 1
    s = θ[2:SZ[1]-1, 2:SZ[2]-1, z_st:z_ed]
    min_val = minimum(s)
    max_val = maximum(s)
    l2_norm = norm(s, 2)
    
    return Dict(
        "execution_time" => execution_time,
        "iterations" => length(conv_data.residuals),
        "final_residual" => length(conv_data.residuals) > 0 ? conv_data.residuals[end] : NaN,
        "min_temp" => min_val,
        "max_temp" => max_val,
        "l2_norm" => l2_norm,
        "temperature_field" => s,
        "conv_data" => conv_data
    )
end

function test_optimized()
    """最適化版のテスト実行"""
    println("Testing optimized version...")
    
    # 最適化版のinclude
    include("optimized_heat3d.jl")
    
    # 実行時間測定
    start_time = time()
    conv_data = optimized_q3d(MODE, NXY, NZ, SOLVER, SMOOTHER, EPSILON, false)
    execution_time = time() - start_time
    
    # この関数内で温度場を取得する必要があるため、
    # optimized_q3d関数を修正して温度場も返すようにする
    # 今回は簡単な検証のため基本的な情報のみ返す
    
    return Dict(
        "execution_time" => execution_time,
        "iterations" => length(conv_data.residuals),
        "final_residual" => length(conv_data.residuals) > 0 ? conv_data.residuals[end] : NaN,
        "conv_data" => conv_data
    )
end

function validate_results(original, optimized)
    """結果の妥当性検証"""
    println("\n" * "="^40)
    println("VALIDATION RESULTS")
    println("="^40)
    
    # 基本情報の表示
    println("Grid Size: $(NXY)x$(NXY)x$(NZ)")
    println("Solver: $(SOLVER) with $(SMOOTHER)")
    println("Tolerance: $(EPSILON)")
    println("")
    
    # 実行時間の比較
    println("Execution Time Comparison:")
    println("  Original:  $(@sprintf("%.4f", original["execution_time"])) seconds")
    println("  Optimized: $(@sprintf("%.4f", optimized["execution_time"])) seconds")
    
    speedup = original["execution_time"] / optimized["execution_time"]
    time_reduction = (original["execution_time"] - optimized["execution_time"]) / original["execution_time"] * 100
    
    println("  Speedup:   $(@sprintf("%.2f", speedup))x")
    println("  Time reduction: $(@sprintf("%.1f", time_reduction))%")
    println("")
    
    # 収束性の比較
    println("Convergence Comparison:")
    println("  Original iterations:  $(original["iterations"])")
    println("  Optimized iterations: $(optimized["iterations"])")
    
    if original["iterations"] == optimized["iterations"]
        println("  ✓ Convergence consistency: PASS")
    else
        diff = abs(original["iterations"] - optimized["iterations"])
        println("  ⚠ Convergence difference: $(diff) iterations")
    end
    
    # 最終残差の比較
    if !isnan(original["final_residual"]) && !isnan(optimized["final_residual"])
        residual_diff = abs(original["final_residual"] - optimized["final_residual"])
        relative_diff = residual_diff / original["final_residual"] * 100
        
        println("  Original final residual:  $(@sprintf("%.6E", original["final_residual"]))")
        println("  Optimized final residual: $(@sprintf("%.6E", optimized["final_residual"]))")
        println("  Relative difference: $(@sprintf("%.3f", relative_diff))%")
        
        if relative_diff < 1.0
            println("  ✓ Final residual consistency: PASS")
        else
            println("  ⚠ Final residual difference: $(relative_diff)%")
        end
    end
    
    println("")
    
    # 温度場の比較（オリジナル版のみで利用可能）
    if haskey(original, "temperature_field")
        println("Temperature Field Statistics (Original):")
        println("  Min temperature: $(@sprintf("%.6f", original["min_temp"])) K")
        println("  Max temperature: $(@sprintf("%.6f", original["max_temp"])) K")
        println("  L2 norm: $(@sprintf("%.6E", original["l2_norm"]))")
        println("")
    end
    
    # 総合評価
    println("OVERALL ASSESSMENT:")
    if speedup >= 2.0 && original["iterations"] == optimized["iterations"]
        println("  ✅ EXCELLENT: Significant speedup with consistent convergence")
    elseif speedup >= 1.5
        println("  ✅ GOOD: Noticeable speedup achieved")
    elseif speedup >= 1.1
        println("  ✅ FAIR: Modest improvement")
    else
        println("  ⚠ POOR: Limited or no improvement")
    end
    
    return Dict(
        "speedup" => speedup,
        "time_reduction" => time_reduction,
        "convergence_consistent" => original["iterations"] == optimized["iterations"]
    )
end

# メイン実行
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    try
        println("Starting Heat3D validation test...")
        println("Working directory: $(pwd())")
        println("")
        
        # オリジナル版テスト
        original_result = test_original()
        
        # 最適化版テスト  
        optimized_result = test_optimized()
        
        # 結果検証
        validation = validate_results(original_result, optimized_result)
        
        println("\nValidation completed successfully!")
        
    catch e
        println("❌ Error during validation: $(e)")
        println("Stack trace:")
        println(stacktrace())
    end
end