include("../Vinokur_Julia/stretch.jl")
using Printf
using BenchmarkTools

"""
C++とJuliaの実行時間比較用ベンチマーク
"""

function benchmark_stretch()
  println("=== C++ vs Julia Performance Benchmark ===\n")
  
  # ベンチマークケース定義（収束保証パラメータ）
  test_cases = [
    (name="Small Uniform (51 nodes)", size=51, x1=0.0, x2=1.0, sp1=0.02, sp2=0.02),
    (name="Medium Uniform (101 nodes)", size=101, x1=0.0, x2=1.0, sp1=0.01, sp2=0.01),
    (name="Large Uniform (201 nodes)", size=201, x1=0.0, x2=1.0, sp1=0.005, sp2=0.005),
    (name="Huge Uniform (501 nodes)", size=501, x1=0.0, x2=1.0, sp1=0.002, sp2=0.002)
  ]
  
  println("C++実行結果:")
  println("-" * repeat("-", 60))
  
  for case in test_cases
    println("$(case.name):")
    @printf("パラメータ: size=%d, x1=%.1f, x2=%.1f, sp1=%.3f, sp2=%.3f\n", 
            case.size, case.x1, case.x2, case.sp1, case.sp2)
    
    # C++実行
    cpp_cmd = `../vinokur_code/test $(case.size) $(case.x1) $(case.x2) $(case.sp1) $(case.sp2)`
    cpp_time = @elapsed run(cpp_cmd)
    @printf("C++実行時間: %.6f 秒\n\n", cpp_time)
  end
  
  println("Julia実行結果:")
  println("-" * repeat("-", 60))
  
  # JIT precompilation
  println("JITウォームアップ...")
  stretch(0.0, 1.0, 0.1, 0.1, 11, debug=false)
  println("ウォームアップ完了\n")
  
  for case in test_cases
    println("$(case.name):")
    @printf("パラメータ: size=%d, x1=%.1f, x2=%.1f, sp1=%.3f, sp2=%.3f\n", 
            case.size, case.x1, case.x2, case.sp1, case.sp2)
    
    # Julia実行（単回）
    julia_time = @elapsed begin
      coord, success = stretch(case.x1, case.x2, case.sp1, case.sp2, case.size, debug=false)
    end
    @printf("Julia実行時間: %.6f 秒 (収束: %s)\n", julia_time, success ? "成功" : "失敗")
    
    # Julia実行（平均）
    julia_bench = @benchmark stretch($case.x1, $case.x2, $case.sp1, $case.sp2, $case.size, debug=false) samples=10 seconds=30
    @printf("Julia平均時間: %.6f 秒 (min=%.6f, max=%.6f)\n\n", 
            mean(julia_bench.times) * 1e-9, 
            minimum(julia_bench.times) * 1e-9, 
            maximum(julia_bench.times) * 1e-9)
  end
end

function detailed_benchmark()
  println("=== 詳細ベンチマーク（中サイズグリッド: 101ノード）===\n")
  
  # パラメータ
  size = 101
  x1, x2 = 0.0, 1.0
  sp1, sp2 = 0.01, 0.01
  
  println("C++詳細測定:")
  cpp_times = Float64[]
  
  for i in 1:10
    cpp_cmd = `../vinokur_code/test $size $x1 $x2 $sp1 $sp2`
    cpp_time = @elapsed run(cpp_cmd)
    push!(cpp_times, cpp_time)
    @printf("実行 %2d: %.6f 秒\n", i, cpp_time)
  end
  
  @printf("C++平均: %.6f 秒 (std=%.6f)\n", mean(cpp_times), std(cpp_times))
  @printf("C++最小: %.6f 秒\n", minimum(cpp_times))
  @printf("C++最大: %.6f 秒\n\n", maximum(cpp_times))
  
  println("Julia詳細測定:")
  
  # ウォームアップ
  stretch(x1, x2, sp1, sp2, size, debug=false)
  
  julia_bench = @benchmark stretch($x1, $x2, $sp1, $sp2, $size, debug=false) samples=100 seconds=60
  
  @printf("Julia平均: %.6f 秒\n", mean(julia_bench.times) * 1e-9)
  @printf("Julia最小: %.6f 秒\n", minimum(julia_bench.times) * 1e-9)
  @printf("Julia最大: %.6f 秒\n", maximum(julia_bench.times) * 1e-9)
  @printf("Julia標準偏差: %.6f 秒\n", std(julia_bench.times) * 1e-9)
  
  # 速度比較
  cpp_avg = mean(cpp_times)
  julia_avg = mean(julia_bench.times) * 1e-9
  
  println("\n=== 速度比較結果 ===")
  if julia_avg < cpp_avg
    speedup = cpp_avg / julia_avg
    @printf("Juliaが %.2fx 高速\n", speedup)
  else
    slowdown = julia_avg / cpp_avg
    @printf("C++が %.2fx 高速\n", slowdown)
  end
end

function memory_benchmark()
  println("\n=== メモリ使用量比較 ===")
  
  size = 201
  x1, x2 = 0.0, 1.0
  sp1, sp2 = 0.005, 0.005
  
  println("Julia メモリ使用量:")
  # ウォームアップ
  stretch(x1, x2, sp1, sp2, 11, debug=false)
  
  # メモリベンチマーク
  memory_bench = @benchmark stretch($x1, $x2, $sp1, $sp2, $size, debug=false) samples=10
  
  @printf("平均メモリ使用量: %.2f KB\n", memory_bench.memory / 1024)
  @printf("アロケーション回数: %d\n", memory_bench.allocs)
end

# メイン実行
function main()
  benchmark_stretch()
  detailed_benchmark()
  memory_benchmark()
  
  println("\nベンチマーク完了")
end

main()