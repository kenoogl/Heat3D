#!/usr/bin/env julia

"""
heat3d_nu.jl プロファイリングスクリプト

このスクリプトはheat3d_nu.jlの性能を分析し、計算時間のホットスポットを特定します。
"""

using Profile
using Printf
using InteractiveUtils

# 元のコードをinclude
include("../q3d/heat3d_nu.jl")

function profile_heat3d_numerical_only()
    println("="^60)
    println("Heat3D Numerical Computation Profiling")
    println("="^60)
    
    # プロファイル用パラメータ
    mode = 3
    NXY = 120  # 実用的なサイズでテスト
    NZ = 31
    solver = "pbicgstab"
    smoother = "gs"
    epsilon = 1.0e-3
    
    println("Profile parameters:")
    println("  Mode: $(mode)")
    println("  Grid: $(NXY)x$(NXY)x$(NZ)")
    println("  Solver: $(solver)")
    println("  Smoother: $(smoother)")
    println("  Tolerance: $(epsilon)")
    println("")
    
    # 数値計算のみのプロファイリング用関数を定義
    function numerical_computation_only()
        global mode = 3
        global itr_tol = epsilon
        
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
        
        # メイン計算のみ（プロット処理は除外）
        conv_data = main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother)
        
        # 基本統計のみ
        z_st = 3
        z_ed = SZ[3] - 1
        s = θ[2:SZ[1]-1, 2:SZ[2]-1, z_st:z_ed]
        min_val = minimum(s)
        max_val = maximum(s)
        l2_norm = norm(s, 2)
        
        return conv_data, min_val, max_val, l2_norm
    end
    
    # プロファイル開始
    println("Starting numerical computation profiling...")
    Profile.clear()
    
    # 最初に一回実行してJITコンパイルを済ませる
    println("JIT compilation run...")
    numerical_computation_only()
    
    # 実際のプロファイル実行
    println("Profiling run...")
    Profile.clear()
    @profile numerical_computation_only()
    
    println("\n" * "="^60)
    println("Profiling Results")
    println("="^60)
    
    # プロファイル結果をファイルに出力
    open("profile_results.txt", "w") do f
        Profile.print(f, format=:flat, sortedby=:count)
    end
    
    # プロファイル結果を表示（上位20件）
    println("Top functions by sample count:")
    println("-"^40)
    Profile.print(format=:flat, sortedby=:count, maxdepth=20)
    
    println("\n" * "-"^40)
    println("Top functions by cumulative time:")
    println("-"^40)
    Profile.print(format=:tree, maxdepth=10)
    
    # 統計情報
    profile_data = Profile.fetch()
    total_samples = length(profile_data[1])
    println("\n" * "="^40)
    println("Profile Statistics:")
    println("  Total samples: $(total_samples)")
    println("  Profile saved to: profile_results.txt")
    println("="^40)
    
    return total_samples
end

function detailed_timing_analysis()
    println("\n" * "="^60)
    println("Detailed Timing Analysis")
    println("="^60)
    
    # パラメータ設定
    mode = 3
    NXY = 60
    NZ = 31
    solver = "pbicgstab"
    smoother = "gs"
    epsilon = 1.0e-4
    
    # 各段階の時間計測
    MX = MY = NXY + 2
    MZ = NZ + 2
    dh = 1.2e-3 / NXY
    dh = round(dh, digits=8)
    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0)
    
    println("Timing individual components:")
    println("-"^40)
    
    # 配列初期化の時間
    t_start = time()
    Z = zeros(Float64, SZ[3])
    ΔZ = zeros(Float64, SZ[3]-1)
    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0
    ID = zeros(UInt8, SZ[1], SZ[2], SZ[3])
    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    θ .= Constant.θ_amb
    t_init = time() - t_start
    println("Array initialization: $(@sprintf("%.4f", t_init)) seconds")
    
    # 前処理時間
    t_start = time()
    preprocess!(SZ, λ, Z, ΔZ, ox, Δh, ID)
    t_preprocess = time() - t_start
    println("Preprocessing: $(@sprintf("%.4f", t_preprocess)) seconds")
    
    # メイン計算時間
    t_start = time()
    conv_data = main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother)
    t_main = time() - t_start
    println("Main computation: $(@sprintf("%.4f", t_main)) seconds")
    
    # 後処理時間（プロット等）
    t_start = time()
    # プロット処理（実際のq3d関数の最後の部分）
    try
        plot_slice_xz_nu(1, λ, 0.3e-3, SZ, ox, Δh, Z, "alpha3_profile.png", "α")
        plot_slice_xz_nu(2, θ, 0.3e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.3_profile.png")
        plot_line_z_nu(θ, SZ, ox, Δh, Z, 0.6e-3, 0.6e-3,"temp3Z_ctr_profile", "Center")
    catch e
        println("Plot error (ignored): $e")
    end
    t_postprocess = time() - t_start
    println("Post-processing: $(@sprintf("%.4f", t_postprocess)) seconds")
    
    total_time = t_init + t_preprocess + t_main + t_postprocess
    println("-"^40)
    println("Total time: $(@sprintf("%.4f", total_time)) seconds")
    println("")
    
    # 相対的な時間配分
    println("Time distribution:")
    println("  Initialization: $(@sprintf("%.1f", t_init/total_time*100))%")
    println("  Preprocessing:  $(@sprintf("%.1f", t_preprocess/total_time*100))%")
    println("  Main computation: $(@sprintf("%.1f", t_main/total_time*100))%")
    println("  Post-processing: $(@sprintf("%.1f", t_postprocess/total_time*100))%")
    
    return Dict(
        "init" => t_init,
        "preprocess" => t_preprocess,
        "main" => t_main,
        "postprocess" => t_postprocess,
        "total" => total_time
    )
end

function memory_analysis()
    println("\n" * "="^60)
    println("Memory Usage Analysis")
    println("="^60)
    
    # パラメータ設定
    mode = 3
    NXY = 240  # フルサイズで確認
    NZ = 31
    
    MX = MY = NXY + 2
    MZ = NZ + 2
    SZ = (MX, MY, MZ)
    
    # 各配列のメモリ使用量計算
    float64_size = sizeof(Float64)
    uint8_size = sizeof(UInt8)
    
    # 主要配列のメモリ使用量
    θ_memory = prod(SZ) * float64_size
    λ_memory = prod(SZ) * float64_size
    ID_memory = prod(SZ) * uint8_size
    Z_memory = SZ[3] * float64_size
    ΔZ_memory = (SZ[3]-1) * float64_size
    
    # PBiCGSTAB用の作業配列
    pbicgstab_arrays = 8  # pcg_p, pcg_p_, pcg_r, pcg_r0, pcg_q, pcg_s, pcg_s_, pcg_t_
    pbicgstab_memory = pbicgstab_arrays * prod(SZ) * float64_size
    
    total_memory = θ_memory + λ_memory + ID_memory + Z_memory + ΔZ_memory + pbicgstab_memory
    
    println("Memory usage for $(NXY)x$(NXY)x$(NZ) grid:")
    println("  Temperature field (θ): $(@sprintf("%.1f", θ_memory/1024/1024)) MB")
    println("  Thermal conductivity (λ): $(@sprintf("%.1f", λ_memory/1024/1024)) MB")
    println("  Material ID: $(@sprintf("%.1f", ID_memory/1024/1024)) MB")
    println("  Z coordinates: $(@sprintf("%.3f", Z_memory/1024)) KB")
    println("  ΔZ array: $(@sprintf("%.3f", ΔZ_memory/1024)) KB")
    println("  PBiCGSTAB work arrays: $(@sprintf("%.1f", pbicgstab_memory/1024/1024)) MB")
    println("  Total estimated: $(@sprintf("%.1f", total_memory/1024/1024)) MB")
    
    return total_memory
end

# メイン実行
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println("Starting Heat3D numerical computation analysis...")
    
    # 数値計算のみのプロファイル実行
    total_samples = profile_heat3d_numerical_only()
    
    # 詳細タイミング解析
    timing_results = detailed_timing_analysis()
    
    # メモリ解析
    memory_usage = memory_analysis()
    
    println("\n" * "="^60)
    println("Analysis Summary")
    println("="^60)
    println("Profile samples collected: $(total_samples)")
    println("Main computation time: $(@sprintf("%.4f", timing_results["main"])) seconds")
    println("Estimated memory usage: $(@sprintf("%.1f", memory_usage/1024/1024)) MB")
    println("\nProfile results saved to: profile_results.txt")
    println("Use ProfileView.jl for graphical analysis")
    
    println("\n" * "="^60)
    println("Key Hotspots (Numerical Computation Only)")
    println("="^60)
    println("This analysis excludes plotting functions and focuses on:")
    println("  - PBiCGSTAB solver iterations")
    println("  - Matrix-vector operations")
    println("  - Boundary condition setup")
    println("  - Material property assignment")
    println("  - Convergence calculations")
end