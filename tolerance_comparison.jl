using Printf
using LinearAlgebra

include("q3d/heat3d_nu.jl")
include("q3d/convergence_history.jl")

"""
tolerance_comparison()
異なる収束判定値(tol)での計算結果を比較
"""
function tolerance_comparison()
    println("=== Tolerance Comparison Study ===")
    
    # テスト対象のtol値（10^-4から10^-10まで）
    tol_values = [1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10]
    
    # 結果保存用
    results = []
    convergence_data_list = []
    
    # 元のtol値を保存
    original_tol = Constant.tol
    
    for (idx, tol_val) in enumerate(tol_values)
        println("\n" * "="^50)
        println("Test $idx/$(length(tol_values)): tol = $(tol_val)")
        println("="^50)
        
        # const.jlの一時的な変更
        const_content = """# common constant
module Constant

const ItrMax = 8000
const tol    = $tol_val
const FloatMin = 1.0e-37
const ω      = 1.0

# 境界条件
const θ_amb = 300.0
const θ_pcb = 300.0
const HT_top = 2.98e-6 # 5 [W/(m^2 K)] / (\\rho C)_silicon > [m/s]
const HT_side = 2.98e-6 # 5 [W/(m^2 K)] / (\\rho C)_silicon > [m/s]
const Q_src = 9.5374e4# 1.6x10^11 [W/m^3] / (\\rho C)_silicon > [K/s]
# 8x10^5 [W/m^2], 厚さ5μm > 1.6 x 10^11 [W/m^3]
# \\rho_silicon= 2330.0 [kg/m^3]
# C_silicon= 720.0 [J/(Kg K)]

end # end of module Constant
"""
        
        # const.jlを一時的に書き換え
        write("q3d/const.jl", const_content)
        
        # モジュールを再読み込み（includeを再実行）
        Base.include(Main, "q3d/const.jl")
        
        # 計算実行
        println("Running simulation with tol = $(tol_val)...")
        start_time = time()
        
        try
            # q3d関数を直接呼び出さずに、内部処理を実行
            mode = 3
            NXY = 60  # 計算時間短縮のため小さなグリッド
            NZ = 15
            solver = "pbicgstab"
            smoother = "gs"
            
            # グリッド設定
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
            θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
            θ .= Constant.θ_amb
            
            # 前処理
            preprocess!(SZ, λ, Z, ΔZ, ox, Δh, ID)
            
            # メイン計算
            conv_data = main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother)
            
            elapsed_time = time() - start_time
            
            # 結果の解析
            z_st = 3
            z_ed = SZ[3]-1
            s = θ[2:SZ[1]-1, 2:SZ[2]-1, z_st:z_ed]
            min_val = minimum(s)
            max_val = maximum(s)
            l2_norm = norm(s, 2)
            
            # log.txtから収束情報を読み取り
            info = get_convergence_info(conv_data)
            
            result = Dict(
                "tol" => tol_val,
                "iterations" => length(conv_data.residuals),
                "final_residual" => length(conv_data.residuals) > 0 ? conv_data.residuals[end] : NaN,
                "min_temp" => min_val,
                "max_temp" => max_val,
                "l2_norm" => l2_norm,
                "elapsed_time" => elapsed_time,
                "converged" => length(conv_data.residuals) > 0 ? conv_data.residuals[end] <= tol_val : false
            )
            
            push!(results, result)
            push!(convergence_data_list, conv_data)
            
            println("✓ Completed: $(result["iterations"]) iterations, final residual = $(@sprintf("%.2E", result["final_residual"]))")
            
        catch e
            println("✗ Error: $e")
            result = Dict(
                "tol" => tol_val,
                "iterations" => -1,
                "final_residual" => NaN,
                "min_temp" => NaN,
                "max_temp" => NaN,
                "l2_norm" => NaN,
                "elapsed_time" => time() - start_time,
                "converged" => false,
                "error" => string(e)
            )
            push!(results, result)
        end
    end
    
    # const.jlを元に戻す
    original_const_content = """# common constant
module Constant

const ItrMax = 8000
const tol    = $original_tol
const FloatMin = 1.0e-37
const ω      = 1.0

# 境界条件
const θ_amb = 300.0
const θ_pcb = 300.0
const HT_top = 2.98e-6 # 5 [W/(m^2 K)] / (\\rho C)_silicon > [m/s]
const HT_side = 2.98e-6 # 5 [W/(m^2 K)] / (\\rho C)_silicon > [m/s]
const Q_src = 9.5374e4# 1.6x10^11 [W/m^3] / (\\rho C)_silicon > [K/s]
# 8x10^5 [W/m^2], 厚さ5μm > 1.6 x 10^11 [W/m^3]
# \\rho_silicon= 2330.0 [kg/m^3]
# C_silicon= 720.0 [J/(Kg K)]

end # end of module Constant
"""
    write("q3d/const.jl", original_const_content)
    
    # 結果の出力
    println("\n" * "="^80)
    println("TOLERANCE COMPARISON RESULTS")
    println("="^80)
    
    println("┌─────────────┬─────────┬──────────────┬──────────────┬──────────────┬─────────────┬──────────┐")
    println("│ Tolerance   │ Iter.   │ Final Resid. │ Min Temp [K] │ Max Temp [K] │ L2 Norm     │ Time [s] │")
    println("├─────────────┼─────────┼──────────────┼──────────────┼──────────────┼─────────────┼──────────┤")
    
    for result in results
        if haskey(result, "error")
            println("│ $(@sprintf("%11.0E", result["tol"])) │ ERROR   │              │              │              │             │          │")
        else
            println("│ $(@sprintf("%11.0E", result["tol"])) │ $(@sprintf("%7d", result["iterations"])) │ $(@sprintf("%12.2E", result["final_residual"])) │ $(@sprintf("%12.6f", result["min_temp"])) │ $(@sprintf("%12.6f", result["max_temp"])) │ $(@sprintf("%11.2E", result["l2_norm"])) │ $(@sprintf("%8.1f", result["elapsed_time"])) │")
        end
    end
    println("└─────────────┴─────────┴──────────────┴──────────────┴──────────────┴─────────────┴──────────┘")
    
    # 収束比較グラフを生成
    if length(convergence_data_list) > 0
        println("\nGenerating convergence comparison plot...")
        try
            compare_convergence(convergence_data_list, "tolerance_comparison.png", show_markers=false)
            println("Convergence comparison plot saved: tolerance_comparison.png")
        catch e
            println("Error generating comparison plot: $e")
        end
    end
    
    # 詳細分析
    println("\n" * "="^80)
    println("DETAILED ANALYSIS")
    println("="^80)
    
    valid_results = filter(r -> !haskey(r, "error") && r["iterations"] > 0, results)
    
    if length(valid_results) > 1
        println("\nTemperature Solution Accuracy:")
        reference = valid_results[end]  # 最も厳しいtolの結果を参照
        for result in valid_results
            temp_diff = abs(result["max_temp"] - reference["max_temp"])
            norm_diff = abs(result["l2_norm"] - reference["l2_norm"]) / reference["l2_norm"] * 100
            println("  tol=$(result["tol"]): ΔT_max = $(@sprintf("%.6f", temp_diff)) K, ΔL2_norm = $(@sprintf("%.3f", norm_diff))%")
        end
        
        println("\nComputational Efficiency:")
        for result in valid_results
            iter_per_sec = result["iterations"] / result["elapsed_time"]
            println("  tol=$(result["tol"]): $(@sprintf("%.1f", iter_per_sec)) iter/sec, $(@sprintf("%.1f", result["elapsed_time"])) sec total")
        end
    end
    
    println("\n✓ Tolerance comparison completed!")
    return results
end

# 実行
tolerance_comparison()