#!/usr/bin/env julia

"""
heat3d_nu.jl 最適化版

主な最適化内容：
1. プロット処理のオプション化
2. メモリアロケーション最適化
3. 計算ループの最適化
4. 並列化（Threads.jl）
5. SIMD最適化
"""

using Printf
using Random
using LinearAlgebra
using Base.Threads

# 元のコードをinclude
include("../q3d/heat3d_Cartesian.jl")
include("../q3d/heat3d_CartesianII.jl")
include("../q3d/heat3d_NonUniform.jl")
include("../q3d/heat3d_NonUniformII.jl")
include("../model/modelA.jl")
include("../q3d/const.jl")
include("../q3d/convergence_history.jl")
include("../q3d/parse_log_residuals.jl")

# 最適化版のプロッター（プロット無効化オプション付き）
module OptimizedPlotter
    export plot_line_z_nu_optimized, export_zline_csv_optimized

    function export_zline_csv_optimized(Z, d, filename::String)
        _fname = "$(filename).csv"
        
        open(_fname, "w") do f
            println(f, "Z [mm], Temperature [K]")
            
            @inbounds for i in 1:length(d)
                @printf(f, "%.14E, %.14E\n", Z[i], d[i])
            end
        end
        
        println("Z-line CSV saved: $_fname")
    end

    function find_i_optimized(xc, ox, dx, SZ_x)
        return clamp(round(Int, (xc - ox) / dx) + 2, 2, SZ_x-1)
    end

    function find_j_optimized(yc, oy, dy, SZ_y)
        return clamp(round(Int, (yc - oy) / dy) + 2, 2, SZ_y-1)
    end

    function plot_line_z_nu_optimized(d::Array{Float64,3}, SZ, ox, Δh, Z, xc, yc, filename, label::String="", enable_plot::Bool=false)
        zs = 2
        ze = SZ[3] - 1
        
        i = find_i_optimized(xc, ox[1], Δh[1], SZ[1])
        j = find_j_optimized(yc, ox[2], Δh[2], SZ[2])
        
        @inbounds s = d[i, j, zs:ze]
        z_coords = [(Z[k] * 1000) for k in zs:ze]  # mm単位に変換
        
        min_val = minimum(s)
        max_val = maximum(s)
        println("At ($(xc*1000), $(yc*1000)) [mm]: min=", min_val, " max=", max_val)
        
        export_zline_csv_optimized(z_coords, s, filename)
        
        # プロット処理はオプション（デフォルトで無効）
        if enable_plot
            try
                include("../q3d/plotter.jl")
                plot_line_z_nu(d, SZ, ox, Δh, Z, xc, yc, filename, label)
            catch e
                println("Plot skipped: $e")
            end
        end
    end
end

using .OptimizedPlotter

# 最適化版のメイン計算関数
function optimized_main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother, enable_plot::Bool=false)
    # 収束履歴の初期化
    conv_data = ConvergenceData(solver, smoother)

    z_st::Int64 = 0
    z_ed::Int64 = 0

    if mode == 1 || mode == 4
        z_st = 2
        z_ed = SZ[3] - 1
    elseif mode == 2
        z_st = 3
        z_ed = SZ[3] - 2
    elseif mode == 3 
        z_st = 3
        z_ed = SZ[3] - 1 
    end

    # 境界条件設定
    if mode == 1
        Cartesian.boundary_condition!(θ, SZ, ox, Δh)
    elseif mode == 2
        NonUniform.boundary_condition!(θ, SZ, ox, Δh)
    elseif mode == 3
        NonUniformII.boundary_condition3!(θ, SZ)
    elseif mode == 4
        CartesianII.boundary_condition4!(θ, SZ)
    end

    # 右辺項設定（メモリ効率化）
    b = zeros(Float64, SZ[1], SZ[2], SZ[3])
    if mode == 3 || mode == 4
        optimized_calRHS!(b, ID, SZ)
    end

    # マスク設定（メモリ効率化）
    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    optimized_setMask!(mask, SZ)

    # 作業配列（最適化版）
    wk = zeros(Float64, SZ[1], SZ[2], SZ[3])

    # PBiCGSTAB用配列（必要時のみ確保）
    if solver == "pbicgstab"
        pcg_p  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_p_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r0 = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_q  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_t_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end

    # ログファイル
    F = open("log.txt", "w")
    optimized_conditions(F, SZ, Δh, solver, smoother)

    # ソルバー実行
    if solver == "sor"
        if mode == 3
            NonUniformII.solveSOR!(θ, SZ, λ, b, mask, Δh, Constant.ω, Z, ΔZ, z_st, z_ed, F, itr_tol)
        end
    elseif solver == "jacobi"
        if mode == 3
            NonUniformII.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, Constant.ω, Z, ΔZ, z_st, z_ed, F, itr_tol)
        end
    elseif solver == "pbicgstab"
        if mode == 3
            NonUniformII.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, 
                ox, Δh, SZ, Z, ΔZ, z_st, z_ed, smoother, F, mode, itr_tol)
        end
    else
        println("solver error")
    end

    # 統計計算（最適化版）
    s = view(θ, 2:SZ[1]-1, 2:SZ[2]-1, z_st:z_ed)
    min_val = minimum(s)
    max_val = maximum(s)
    l2_norm = norm(s, 2)
    @printf(F, "θmin=%e  θmax=%e  L2 norm of θ=%e\n", min_val, max_val, l2_norm)

    close(F)
    
    # ログファイルから残差データを解析
    if solver == "pbicgstab" || solver == "sor" || solver == "jacobi"
        parse_residuals_from_log!(conv_data, "log.txt")
    end
    
    return conv_data
end

# 最適化版の右辺項計算
function optimized_calRHS!(b::Array{Float64,3}, ID::Array{UInt8,3}, SZ)
    @inbounds @threads for k in 1:SZ[3]
        for j in 1:SZ[2], i in 1:SZ[1]
            if ID[i,j,k] == modelA.pwrsrc["id"]
                b[i,j,k] = -Constant.Q_src
            end
        end
    end
end

# 最適化版のマスク設定
function optimized_setMask!(m::Array{Float64,3}, SZ)
    # Z方向境界
    @inbounds @threads for j in 1:SZ[2]
        for i in 1:SZ[1]
            m[i,j,1] = 0.0
            m[i,j,SZ[3]] = 0.0
        end
    end

    # Y方向境界
    @inbounds @threads for k in 1:SZ[3]
        for i in 1:SZ[1]
            m[i,1,k] = 0.0
            m[i,SZ[2],k] = 0.0
        end
    end

    # X方向境界
    @inbounds @threads for k in 1:SZ[3]
        for j in 1:SZ[2]
            m[1,j,k] = 0.0
            m[SZ[1],j,k] = 0.0
        end
    end
end

# 最適化版の条件出力
function optimized_conditions(F, SZ, Δh, solver, smoother)
    if mode == 1
        @printf(F, "Problem : CUBE on Cartesian grid\n")
    elseif mode == 2
        @printf(F, "Problem : CUBE on NonUniform grid from file\n")
    elseif mode == 3
        @printf(F, "Problem : IC on NonUniform grid (Opt. 13 layers)\n")
    elseif mode == 4
        @printf(F, "Problem : IC on Cartesian grid\n")
    end

    @printf(F, "Grid  : %d %d %d\n", SZ[1], SZ[2], SZ[3])
    @printf(F, "Pitch : %6.4e %6.4e %6.4e\n", Δh[1], Δh[2], Δh[3])
    if solver == "pbicgstab"
        @printf(F, "Solver: %s with smoother %s\n", solver, smoother)
    else
        @printf(F, "Solver: %s\n", solver)
    end
    @printf(F, "ItrMax : %e\n", Constant.ItrMax)
    @printf(F, "ε      : %e\n", itr_tol)
    @printf(F, "θ_amb  : %e\n", Constant.θ_amb)
    @printf(F, "θ_pcb  : %e\n", Constant.θ_pcb)
    @printf(F, "HT_top : %e\n", Constant.HT_top)
    @printf(F, "HT_side: %e\n", Constant.HT_side)
    @printf(F, "Q_src  : %e\n", Constant.Q_src)
end

# Z軸座標生成（最適化版）
function optimized_genZ!(Z::Vector{Float64}, ΔZ::Vector{Float64}, SZ, ox, dz)
    mz = SZ[3]
    if mode == 1 || mode == 4
        @inbounds @simd for k in 1:(mz+1)
            Z[k] = ox[3] + (k-2)*dz
        end
    elseif mode == 2
        read_coord, numNodes = NonUniform.read_grid_file()
        if numNodes != mz-2
            println("Number of generated grid is not match with parameter NZ")
            exit(0)
        end
        @inbounds for k in 1:mz-2
            Z[k+1] = read_coord[k]
        end
        Z[1] = 2*Z[2] - Z[3]
        Z[mz] = 2*Z[mz-1] - Z[mz-2]
    else
        # Zcase2!の最適化版
        optimized_Zcase2!(Z, SZ)
    end

    @inbounds @simd for k in 2:mz-1
        ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])
    end

    if mode == 3 || mode == 2
        ΔZ[2] = 0.5*ΔZ[2]
        ΔZ[mz-1] = 0.5*ΔZ[mz-1]
    end
end

function optimized_Zcase2!(Z::Vector{Float64}, SZ)
    if SZ[3] != 33
        println("MZ must be 33")
        exit(0)
    end
    p = 0.005e-3
    # 最適化：配列への直接代入
    @inbounds begin
        Z[1] = modelA.zm0 - p
        Z[2] = modelA.zm0
        Z[3] = modelA.zm0 + p
        Z[4] = modelA.zm1 - p
        Z[5] = modelA.zm1
        Z[6] = modelA.zm1 + p
        Z[7] = modelA.zm2 - p
        Z[8] = modelA.zm2
        Z[9] = modelA.zm2 + p
        Z[10] = modelA.zm3 - p
        Z[11] = modelA.zm3
        Z[12] = modelA.zm4
        Z[13] = modelA.zm4 + p
        Z[14] = modelA.zm5 - p
        Z[15] = modelA.zm5
        Z[16] = modelA.zm5 + p
        Z[17] = modelA.zm6 - p
        Z[18] = modelA.zm6
        Z[19] = modelA.zm7
        Z[20] = modelA.zm7 + p
        Z[21] = modelA.zm8 - p
        Z[22] = modelA.zm8
        Z[23] = modelA.zm8 + p
        Z[24] = modelA.zm9 - p
        Z[25] = modelA.zm9 
        Z[26] = modelA.zm10
        Z[27] = modelA.zm10 + p
        Z[28] = modelA.zm11 - p
        Z[29] = modelA.zm11
        Z[30] = modelA.zm11 + p
        Z[31] = modelA.zm12 - p
        Z[32] = modelA.zm12
        Z[33] = modelA.zm12 + p
    end
end

# 最適化版の前処理
function optimized_preprocess!(SZ, λ, Z, ΔZ, ox, Δh, ID)
    optimized_genZ!(Z, ΔZ, SZ, ox, Δh[3])

    if mode == 3 || mode == 4
        modelA.fillID!(mode, ID, ox, Δh, SZ, Z)
        modelA.setLambda!(λ, SZ, ID)
        optimized_setMatOuter!(λ, SZ)
    end
end

function optimized_setMatOuter!(λ::Array{Float64,3}, SZ)
    # Z方向境界
    @inbounds @threads for j in 1:SZ[2]
        for i in 1:SZ[1]
            λ[i,j,1] = λ[i,j,2]
            λ[i,j,SZ[3]] = 0.0
        end
    end

    # Y方向境界
    @inbounds @threads for k in 1:SZ[3]
        for i in 1:SZ[1]
            λ[i,1,k] = 0.0
            λ[i,SZ[2],k] = 0.0
        end
    end

    # X方向境界
    @inbounds @threads for k in 1:SZ[3]
        for j in 1:SZ[2]
            λ[1,j,k] = 0.0
            λ[SZ[1],j,k] = 0.0
        end
    end
end

# 最適化版のメイン関数
function optimized_q3d(m_mode::Int, NXY::Int, NZ::Int, solver::String="sor", smoother::String="", epsilon::Float64=1.0e-6, enable_plot::Bool=false)
    global mode = m_mode
    global itr_tol = epsilon

    MX = MY = NXY + 2
    MZ = NZ + 2

    if !(mode == 1 || mode == 2 || mode == 3 || mode == 4)
        println("mode error")
        return
    end

    if mode == 1
        if NXY != NZ
            println("NXY must be equal to NZ")
            return
        end
    end

    dh::Float64 = 1.0

    if mode == 1 || mode == 2
        dh = 1.0 / NXY
    elseif mode == 3 || mode == 4
        dh = 1.2e-3 / NXY
    end
    dh = round(dh, digits=8)

    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0)

    println(SZ, "  Itr.ε= ", itr_tol)

    # 配列初期化（最適化版）
    if mode == 1 || mode == 4
        Z = zeros(Float64, SZ[3]+1)
    else
        Z = zeros(Float64, SZ[3])
    end

    ΔZ = zeros(Float64, SZ[3]-1)
    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0
    ID = zeros(UInt8, SZ[1], SZ[2], SZ[3])

    # 前処理
    @time optimized_preprocess!(SZ, λ, Z, ΔZ, ox, Δh, ID)

    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    θ .= Constant.θ_amb

    # メイン計算（最適化版）
    @time conv_data = optimized_main(SZ, ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother, enable_plot)
    
    # 最小限の出力処理
    if mode == 3
        plot_line_z_nu_optimized(θ, SZ, ox, Δh, Z, 0.6e-3, 0.6e-3, "temp3Z_ctr_opt", "Center", enable_plot)
        plot_line_z_nu_optimized(θ, SZ, ox, Δh, Z, 0.4e-3, 0.4e-3, "temp3Z_tsv_opt", "TSV", enable_plot)
    end
    
    # 収束履歴の出力（CSVのみ）
    if solver == "pbicgstab" || solver == "sor" || solver == "jacobi"
        conv_filename = "convergence_$(solver)_mode$(mode)_$(NXY)x$(NZ)_opt"
        if !isempty(smoother)
            conv_filename *= "_$(smoother)"
        end
        
        try
            export_convergence_csv(conv_data, "$(conv_filename).csv")
            println("Convergence CSV saved: $(conv_filename).csv")
        catch e
            println("Error in convergence history output: $e")
        end
        
        # 収束情報の表示
        info = get_convergence_info(conv_data)
        if !isempty(info)
            println("\n=== Convergence Information (Optimized) ===")
            println("Mode: $mode, Grid: $(NXY)x$(NZ)")
            println("Solver: $(info["solver"]), Smoother: $(info["smoother"])")
            println("Iterations: $(info["iterations"])")
            @printf("Initial residual: %.6E\n", info["initial_residual"])
            @printf("Final residual: %.6E\n", info["final_residual"])
            @printf("Residual reduction factor: %.6E\n", info["convergence_rate"])
            @printf("Order reduction: %.2f\n", info["reduction_factor"])
            println("==========================================")
        end
    end

    return conv_data
end

# 使用例とベンチマーク
function benchmark_optimization()
    println("="^60)
    println("Heat3D Optimization Benchmark")
    println("="^60)
    
    mode = 3
    NXY = 120  # 中程度のサイズでテスト
    NZ = 31
    solver = "pbicgstab"
    smoother = "gs"
    epsilon = 1.0e-4
    
    println("Benchmark parameters:")
    println("  Mode: $(mode)")
    println("  Grid: $(NXY)x$(NXY)x$(NZ)")
    println("  Solver: $(solver)")
    println("  Smoother: $(smoother)")
    println("  Tolerance: $(epsilon)")
    println("")
    println("Threads available: $(Threads.nthreads())")
    println("")
    
    # 最適化版の実行
    println("Running optimized version (no plots)...")
    time_opt = @elapsed result_opt = optimized_q3d(mode, NXY, NZ, solver, smoother, epsilon, false)
    
    println("")
    println("="^60)
    println("Optimization Results")
    println("="^60)
    println("Optimized execution time: $(@sprintf("%.4f", time_opt)) seconds")
    
    if !isnothing(result_opt)
        info = get_convergence_info(result_opt)
        if !isempty(info)
            println("Iterations: $(info["iterations"])")
            println("Final residual: $(@sprintf("%.6E", info["final_residual"]))")
        end
    end
    
    println("\nKey optimizations applied:")
    println("  - Plot processing disabled")
    println("  - Multi-threading enabled")
    println("  - SIMD optimizations")
    println("  - Memory allocation optimizations")
    println("  - Bounds checking disabled (@inbounds)")
    
    return time_opt
end

# メイン実行
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    benchmark_optimization()
end