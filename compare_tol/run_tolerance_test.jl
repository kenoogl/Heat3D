#!/usr/bin/env julia

"""
許容誤差比較実験 - 手動実行スクリプト

Usage:
  julia run_tolerance_test.jl <tolerance_value>
  
Example:
  julia run_tolerance_test.jl 1.0e-4
  julia run_tolerance_test.jl 1.0e-5
"""

using Printf

# 元のコードをinclude
include("../q3d/heat3d_nu.jl")

function run_tolerance_test(epsilon::Float64)
    println("="^60)
    println("Tolerance Test: ε = $(epsilon)")
    println("="^60)
    
    # パラメータ設定
    mode = 3
    NXY = 240
    NZ = 31
    solver = "pbicgstab"
    smoother = "gs"
    
    println("Parameters:")
    println("  Mode: $(mode)")
    println("  Grid: $(NXY)x$(NXY)x$(NZ)")
    println("  Solver: $(solver)")
    println("  Smoother: $(smoother)")
    println("  Tolerance: $(epsilon)")
    println("")
    
    # 実行時間計測開始
    start_time = time()
    
    try
        # メイン計算実行
        println("Starting calculation...")
        result = q3d(mode, NXY, NZ, solver, smoother, epsilon)
        
        elapsed_time = time() - start_time
        
        println("")
        println("="^60)
        println("Calculation Completed")
        println("="^60)
        println("Elapsed time: $(@sprintf("%.2f", elapsed_time)) seconds")
        
        # ファイル名に許容誤差値を含める
        epsilon_str = replace(string(epsilon), "." => "_", "-" => "m")
        
        # Z方向ライン出力ファイル名を変更
        original_files = [
            "temp3Z_ctr.csv",
            "temp3Z_tsv.csv"
        ]
        
        new_files = [
            "temp3Z_ctr_eps$(epsilon_str).csv",
            "temp3Z_tsv_eps$(epsilon_str).csv"
        ]
        
        # ファイル名変更
        for (old_file, new_file) in zip(original_files, new_files)
            if isfile(old_file)
                mv(old_file, new_file)
                println("Renamed: $(old_file) -> $(new_file)")
            else
                println("Warning: $(old_file) not found - skipping rename")
            end
        end
        
        # 収束履歴ファイルもリネーム
        conv_files = filter(f -> startswith(f, "convergence_") && endswith(f, ".csv"), readdir("."))
        for conv_file in conv_files
            base_name = replace(conv_file, ".csv" => "")
            new_conv_file = "$(base_name)_eps$(epsilon_str).csv"
            mv(conv_file, new_conv_file)
            println("Renamed: $(conv_file) -> $(new_conv_file)")
        end
        
        # PNG収束履歴ファイルもリネーム
        conv_png_files = filter(f -> startswith(f, "convergence_") && endswith(f, ".png"), readdir("."))
        for conv_png_file in conv_png_files
            base_name = replace(conv_png_file, ".png" => "")
            new_conv_png_file = "$(base_name)_eps$(epsilon_str).png"
            mv(conv_png_file, new_conv_png_file)
            println("Renamed: $(conv_png_file) -> $(new_conv_png_file)")
        end
        
        println("")
        println("Key output files:")
        println("  - $(new_files[1]) (center line temperature)")
        println("  - $(new_files[2]) (TSV line temperature)")
        if !isempty(conv_files)
            conv_renamed = replace(conv_files[1], ".csv" => "_eps$(epsilon_str).csv")
            println("  - $(conv_renamed) (convergence history)")
        end
        
        return true
        
    catch e
        elapsed_time = time() - start_time
        println("")
        println("ERROR: Calculation failed after $(@sprintf("%.2f", elapsed_time)) seconds")
        println("Error message: $(e)")
        return false
    end
end

# コマンドライン引数の処理
if length(ARGS) >= 1
    epsilon = parse(Float64, ARGS[1])
    success = run_tolerance_test(epsilon)
    exit(success ? 0 : 1)
else
    println("Usage: julia run_tolerance_test.jl <tolerance_value>")
    println("")
    println("Example:")
    println("  julia run_tolerance_test.jl 1.0e-3")
    println("  julia run_tolerance_test.jl 1.0e-4")
    println("  julia run_tolerance_test.jl 1.0e-5")
    println("  julia run_tolerance_test.jl 1.0e-6")
    println("  julia run_tolerance_test.jl 1.0e-7")
    println("  julia run_tolerance_test.jl 1.0e-8")
    exit(1)
end