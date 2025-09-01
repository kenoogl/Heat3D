#!/usr/bin/env julia

"""
最適化版Heat3D小規模テスト (60x60x31)
"""

using Printf
using LinearAlgebra

println("="^50)
println("Heat3D Optimized Version - Small Test")
println("="^50)

# 最適化版のinclude
include("optimized_heat3d.jl")

# テストパラメータ
const MODE = 3
const NXY = 60
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

println("Executing optimized version...")
start_time = time()

try
    # 最適化版実行
    conv_data = optimized_q3d(MODE, NXY, NZ, SOLVER, SMOOTHER, EPSILON, false)
    
    execution_time = time() - start_time
    
    println("="^50)
    println("OPTIMIZED VERSION RESULTS")
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
    println("Key optimizations applied:")
    println("  ✓ Plot processing disabled")
    println("  ✓ Multi-threading (@threads)")
    println("  ✓ SIMD optimizations (@simd)")
    println("  ✓ Bounds checking disabled (@inbounds)")
    println("  ✓ Memory allocation optimizations")
    
    # 結果保存
    global optimized_result = Dict(
        "execution_time" => execution_time,
        "iterations" => iterations,
        "initial_residual" => initial_residual,
        "final_residual" => final_residual,
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

println("\nOptimized version test completed.")