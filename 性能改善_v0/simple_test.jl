#!/usr/bin/env julia

"""
最適化版の基本動作確認と性能測定
"""

using Printf
using LinearAlgebra

println("="^50)
println("Heat3D Optimized Version Test")
println("="^50)

# 最適化版のみテスト
include("optimized_heat3d.jl")

# 小サイズでテスト
mode = 3
NXY = 60
NZ = 31
solver = "pbicgstab"
smoother = "gs"
epsilon = 1.0e-4

println("Test parameters:")
println("  Grid: $(NXY)x$(NXY)x$(NZ)")
println("  Solver: $(solver)")
println("  Smoother: $(smoother)")
println("  Tolerance: $(epsilon)")
println("  Threads: $(Threads.nthreads())")
println("")

# 実行時間測定
println("Executing optimized version...")
start_time = time()

try
    result = optimized_q3d(mode, NXY, NZ, solver, smoother, epsilon, false)
    
    execution_time = time() - start_time
    
    println("="^50)
    println("RESULTS")
    println("="^50)
    println("✅ Execution completed successfully!")
    println("Execution time: $(@sprintf("%.4f", execution_time)) seconds")
    
    if !isnothing(result)
        iterations = length(result.residuals)
        if iterations > 0
            final_residual = result.residuals[end]
            initial_residual = result.residuals[1]
            
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
    end
    
    println("")
    println("Key optimizations applied:")
    println("  ✓ Plot processing disabled")
    println("  ✓ Multi-threading (@threads)")
    println("  ✓ SIMD optimizations (@simd)")
    println("  ✓ Bounds checking disabled (@inbounds)")
    println("  ✓ Memory allocation optimizations")
    
catch e
    execution_time = time() - start_time
    println("❌ Error occurred after $(@sprintf("%.4f", execution_time)) seconds")
    println("Error: $(e)")
    println("Stack trace:")
    for (i, frame) in enumerate(stacktrace())
        println("  $i: $frame")
        if i >= 10  # 最初の10フレームのみ表示
            break
        end
    end
end

println("\nTest completed.")