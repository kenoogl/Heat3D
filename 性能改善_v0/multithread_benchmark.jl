#!/usr/bin/env julia

"""
Heat3D マルチスレッド性能ベンチマーク

1, 2, 4, 8スレッドでの性能比較を行う
"""

using Printf
using LinearAlgebra
using Base.Threads

println("="^60)
println("Heat3D Multithread Performance Benchmark")
println("="^60)

# テストパラメータ
const MODE = 3
const NXY = 120  # 中規模テスト
const NZ = 31
const SOLVER = "pbicgstab"
const SMOOTHER = "gs"
const EPSILON = 1.0e-4
const TEST_THREADS = [1, 2, 4, 8]

# 結果保存用
results = Dict()

function run_optimized_test(num_threads::Int)
    """指定されたスレッド数で最適化版をテスト"""
    println("\n" * "="^50)
    println("Testing with $(num_threads) thread(s)")
    println("="^50)
    
    # 新しいJuliaプロセスでスレッド数を設定して実行
    test_script = """
    using Printf
    using LinearAlgebra
    using Base.Threads
    
    # 最適化版のinclude
    include("optimized_heat3d.jl")
    
    # テストパラメータ
    const MODE = $MODE
    const NXY = $NXY
    const NZ = $NZ
    const SOLVER = "$SOLVER"
    const SMOOTHER = "$SMOOTHER"
    const EPSILON = $EPSILON
    
    println("Threads available: \$(Threads.nthreads())")
    println("Grid: \$(NXY)x\$(NXY)x\$(NZ)")
    println()
    
    # 実行時間測定
    start_time = time()
    
    try
        conv_data = optimized_q3d(MODE, NXY, NZ, SOLVER, SMOOTHER, EPSILON, false)
        execution_time = time() - start_time
        
        iterations = length(conv_data.residuals)
        final_residual = conv_data.residuals[end]
        initial_residual = conv_data.residuals[1]
        
        println("✅ Execution completed successfully!")
        println("Execution time: \$(@sprintf("%.4f", execution_time)) seconds")
        println("Iterations: \$(iterations)")
        println("Initial residual: \$(@sprintf("%.6E", initial_residual))")
        println("Final residual: \$(@sprintf("%.6E", final_residual))")
        
        # 結果をCSVに保存
        open("thread_\$(Threads.nthreads())_result.csv", "w") do f
            println(f, "threads,execution_time,iterations,initial_residual,final_residual")
            println(f, "\$(Threads.nthreads()),\$(execution_time),\$(iterations),\$(initial_residual),\$(final_residual)")
        end
        
    catch e
        execution_time = time() - start_time
        println("❌ Error occurred after \$(@sprintf("%.4f", execution_time)) seconds")
        println("Error: \$(e)")
        
        # エラー結果をCSVに保存
        open("thread_\$(Threads.nthreads())_result.csv", "w") do f
            println(f, "threads,execution_time,iterations,initial_residual,final_residual")
            println(f, "\$(Threads.nthreads()),\$(execution_time),ERROR,ERROR,ERROR")
        end
    end
    """
    
    # 一時スクリプトファイルを作成
    script_file = "temp_test_$(num_threads)t.jl"
    open(script_file, "w") do f
        write(f, test_script)
    end
    
    # Juliaコマンドを実行
    cmd = `julia -t $(num_threads) $(script_file)`
    
    println("Executing command: julia -t $(num_threads) $(script_file)")
    run(cmd)
    
    # 一時ファイルを削除
    rm(script_file, force=true)
    
    # 結果を読み込み
    result_file = "thread_$(num_threads)_result.csv"
    if isfile(result_file)
        lines = readlines(result_file)
        if length(lines) >= 2
            data = split(lines[2], ",")
            result = Dict(
                "threads" => parse(Int, data[1]),
                "execution_time" => data[2] == "ERROR" ? "ERROR" : parse(Float64, data[2]),
                "iterations" => data[3] == "ERROR" ? "ERROR" : parse(Int, data[3]),
                "initial_residual" => data[4] == "ERROR" ? "ERROR" : parse(Float64, data[4]),
                "final_residual" => data[5] == "ERROR" ? "ERROR" : parse(Float64, data[5])
            )
            rm(result_file, force=true)
            return result
        end
    end
    
    return Dict("threads" => num_threads, "execution_time" => "ERROR")
end

function run_benchmark_suite()
    """全スレッド数でのベンチマーク実行"""
    
    println("Benchmark Configuration:")
    println("  Grid Size: $(NXY)x$(NXY)x$(NZ)")
    println("  Solver: $(SOLVER)")
    println("  Smoother: $(SMOOTHER)")
    println("  Tolerance: $(EPSILON)")
    println("  Test Thread Counts: $(TEST_THREADS)")
    println()
    
    for num_threads in TEST_THREADS
        result = run_optimized_test(num_threads)
        results[num_threads] = result
        
        # 進捗表示
        if result["execution_time"] != "ERROR"
            println("  $(num_threads) threads: $(result["execution_time"]) seconds")
        else
            println("  $(num_threads) threads: ERROR")
        end
    end
    
    return results
end

function generate_summary_report(results)
    """性能改善サマリーレポートの生成"""
    
    println("\n" * "="^60)
    println("MULTITHREAD PERFORMANCE SUMMARY")
    println("="^60)
    
    # 成功した結果のみ集計
    valid_results = filter(p -> p.second["execution_time"] != "ERROR", results)
    
    if isempty(valid_results)
        println("No successful benchmarks to summarize.")
        return
    end
    
    println("| Threads | Execution Time [s] | Iterations | Speedup | Efficiency |")
    println("|---------|-------------------|------------|---------|-------------|")
    
    # 1スレッドの結果をベースラインとする
    baseline_time = haskey(valid_results, 1) ? valid_results[1]["execution_time"] : 
                   minimum([r["execution_time"] for r in values(valid_results)])
    
    for num_threads in sort(collect(keys(valid_results)))
        result = valid_results[num_threads]
        exec_time = result["execution_time"]
        iterations = result["iterations"]
        speedup = baseline_time / exec_time
        efficiency = speedup / num_threads * 100
        
        println("| $(num_threads) | $(@sprintf("%.4f", exec_time)) | $(iterations) | $(@sprintf("%.2f", speedup))x | $(@sprintf("%.1f", efficiency))% |")
    end
    
    println()
    
    # 統計情報
    speedups = [baseline_time / valid_results[t]["execution_time"] for t in sort(collect(keys(valid_results)))]
    
    println("Performance Analysis:")
    println("  Baseline (1 thread): $(@sprintf("%.4f", baseline_time)) seconds")
    println("  Best speedup: $(@sprintf("%.2f", maximum(speedups)))x")
    
    # スケーラビリティ分析
    max_threads = maximum(keys(valid_results))
    if max_threads in keys(valid_results)
        max_speedup = baseline_time / valid_results[max_threads]["execution_time"]
        max_efficiency = max_speedup / max_threads * 100
        println("  $(max_threads)-thread efficiency: $(@sprintf("%.1f", max_efficiency))%")
    end
    
    # 結果をCSVファイルに保存
    open("multithread_results.csv", "w") do f
        println(f, "threads,execution_time,iterations,speedup,efficiency")
        for num_threads in sort(collect(keys(valid_results)))
            result = valid_results[num_threads]
            exec_time = result["execution_time"]
            iterations = result["iterations"]
            speedup = baseline_time / exec_time
            efficiency = speedup / num_threads * 100
            println(f, "$(num_threads),$(exec_time),$(iterations),$(speedup),$(efficiency)")
        end
    end
    
    println("\nResults saved to: multithread_results.csv")
    
    return valid_results
end

# メイン実行
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println("Starting Heat3D multithread performance benchmark...")
    println("Working directory: $(pwd())")
    
    # ベンチマーク実行
    benchmark_results = run_benchmark_suite()
    
    # サマリーレポート生成
    generate_summary_report(benchmark_results)
    
    println("\nMultithread benchmark completed!")
end