#!/usr/bin/env julia

"""
Heat3D Performance Benchmark Comparison

オリジナル版と最適化版の性能を比較するベンチマークスクリプト
"""

using Printf
using BenchmarkTools
using Statistics

# 現在のディレクトリをセット
cd(@__DIR__)

println("="^60)
println("Heat3D Performance Benchmark Comparison")
println("="^60)

# ベンチマーク設定
const TEST_PARAMS = [
    (mode=3, NXY=60,  NZ=31, name="Small"),   # 開発用
    (mode=3, NXY=120, NZ=31, name="Medium"),  # 性能測定用
    (mode=3, NXY=240, NZ=31, name="Large"),   # 実用規模
]

const SOLVER = "pbicgstab"
const SMOOTHER = "gs"
const EPSILON = 1.0e-4

function benchmark_original(mode, NXY, NZ)
    """オリジナル版のベンチマーク"""
    # オリジナルコードをinclude
    include("../q3d/heat3d_nu.jl")
    
    # プロット無効化のため、元のq3d関数をオーバーライド
    function q3d_no_plot(m_mode::Int, NXY::Int, NZ::Int, solver::String="sor", smoother::String="", epsilon::Float64=1.0e-6)
        global mode = m_mode
        global itr_tol = epsilon

        MX = MY = NXY + 2
        MZ = NZ + 2

        if !(mode == 1 || mode == 2 || mode == 3 || mode == 4)
            println("mode error")
            return
        end

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

        preprocess!(SZ, λ, Z, ΔZ, ox, Δh, ID)

        θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        θ .= Constant.θ_amb

        # メイン計算のみ（プロット処理は除外）
        conv_data = main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother)
        
        return conv_data
    end
    
    return q3d_no_plot(mode, NXY, NZ, SOLVER, SMOOTHER, EPSILON)
end

function benchmark_optimized(mode, NXY, NZ)
    """最適化版のベンチマーク"""
    include("optimized_heat3d.jl")
    return optimized_q3d(mode, NXY, NZ, SOLVER, SMOOTHER, EPSILON, false)  # プロット無効
end

function run_benchmark_suite()
    """ベンチマークスイートの実行"""
    
    println("Benchmark Configuration:")
    println("  Solver: $(SOLVER)")
    println("  Smoother: $(SMOOTHER)")
    println("  Tolerance: $(EPSILON)")
    println("  Threads available: $(Threads.nthreads())")
    println("")
    
    results = []
    
    for params in TEST_PARAMS
        println("-"^60)
        println("$(params.name) Grid Test ($(params.NXY)x$(params.NXY)x$(params.NZ))")
        println("-"^60)
        
        # メモリ使用量推定
        estimated_memory = (params.NXY + 2)^2 * (params.NZ + 2) * 9 * 8 / 1024 / 1024  # 9配列 * 8bytes
        println("Estimated memory usage: $(@sprintf("%.1f", estimated_memory)) MB")
        
        try
            # オリジナル版ベンチマーク
            println("\n[Original Version]")
            original_time = @elapsed begin
                original_result = benchmark_original(params.mode, params.NXY, params.NZ)
            end
            
            original_iterations = length(original_result.residuals)
            println("  Execution time: $(@sprintf("%.4f", original_time)) seconds")
            println("  Iterations: $(original_iterations)")
            
            # 最適化版ベンチマーク
            println("\n[Optimized Version]")
            optimized_time = @elapsed begin
                optimized_result = benchmark_optimized(params.mode, params.NXY, params.NZ)
            end
            
            optimized_iterations = length(optimized_result.residuals)
            println("  Execution time: $(@sprintf("%.4f", optimized_time)) seconds")
            println("  Iterations: $(optimized_iterations)")
            
            # 性能改善の計算
            speedup = original_time / optimized_time
            time_reduction = (original_time - optimized_time) / original_time * 100
            
            println("\n[Performance Improvement]")
            println("  Speedup: $(@sprintf("%.2f", speedup))x")
            println("  Time reduction: $(@sprintf("%.1f", time_reduction))%")
            
            # 収束性の確認
            if original_iterations == optimized_iterations
                println("  ✓ Convergence consistent")
            else
                println("  ⚠ Convergence difference: $(original_iterations) vs $(optimized_iterations)")
            end
            
            # 結果を保存
            push!(results, Dict(
                "name" => params.name,
                "grid_size" => "$(params.NXY)x$(params.NXY)x$(params.NZ)",
                "original_time" => original_time,
                "optimized_time" => optimized_time,
                "speedup" => speedup,
                "time_reduction" => time_reduction,
                "original_iterations" => original_iterations,
                "optimized_iterations" => optimized_iterations,
                "memory_mb" => estimated_memory
            ))
            
        catch e
            println("  ✗ Error during benchmark: $(e)")
            push!(results, Dict(
                "name" => params.name,
                "error" => string(e)
            ))
        end
        
        println("")
    end
    
    return results
end

function generate_summary_report(results)
    """性能改善サマリーレポートの生成"""
    
    println("="^60)
    println("PERFORMANCE IMPROVEMENT SUMMARY")
    println("="^60)
    
    # 成功した結果のみ集計
    valid_results = filter(r -> !haskey(r, "error"), results)
    
    if isempty(valid_results)
        println("No successful benchmarks to summarize.")
        return
    end
    
    println("| Grid Size | Original [s] | Optimized [s] | Speedup | Time Reduction | Memory [MB] |")
    println("|-----------|--------------|---------------|---------|----------------|-------------|")
    
    for result in valid_results
        println("| $(result["grid_size"]) | $(@sprintf("%.4f", result["original_time"])) | $(@sprintf("%.4f", result["optimized_time"])) | $(@sprintf("%.2f", result["speedup"]))x | $(@sprintf("%.1f", result["time_reduction"]))% | $(@sprintf("%.1f", result["memory_mb"])) |")
    end
    
    println("")
    
    # 統計情報
    speedups = [r["speedup"] for r in valid_results]
    time_reductions = [r["time_reduction"] for r in valid_results]
    
    println("Statistical Summary:")
    println("  Average speedup: $(@sprintf("%.2f", mean(speedups)))x")
    println("  Max speedup: $(@sprintf("%.2f", maximum(speedups)))x")
    println("  Average time reduction: $(@sprintf("%.1f", mean(time_reductions)))%")
    println("  Max time reduction: $(@sprintf("%.1f", maximum(time_reductions)))%")
    
    println("")
    println("Key Optimizations Applied:")
    println("  - Plot processing disabled")
    println("  - Multi-threading enabled (@threads)")
    println("  - SIMD optimizations (@simd)")
    println("  - Bounds checking disabled (@inbounds)")
    println("  - Memory allocation optimizations")
    println("  - Loop ordering optimizations")
    
    return valid_results
end

# メイン実行
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println("Starting Heat3D performance benchmark comparison...")
    println("Working directory: $(pwd())")
    println("")
    
    # ベンチマーク実行
    benchmark_results = run_benchmark_suite()
    
    # サマリーレポート生成
    generate_summary_report(benchmark_results)
    
    println("\nBenchmark completed!")
end