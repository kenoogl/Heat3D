using Printf
using Random
using LinearAlgebra

include("heat3d_Cartesian.jl")
include("heat3d_NonUniform.jl")

const ω      = 1.0

# 元のプロッタと新しいプロッタの両方をインクルード
include("plotter.jl")         # 元の関数（plot_slice, plot_slice2）
include("plotter_nu.jl")      # NonUniform対応関数

#=
@brief マスク指定
@param [in]     m      マスク
@param [in]     SZ     配列長
@note  非計算セルは0.0、計算セルには1.0
        ディリクレ型　温度指定 : mask_BC=0.0, θ=θ_BC, λ_BC=λ_inner
        ノイマン型　断熱境界 : λ_BC=0.0
=#
function setMask!(m::Array{Float64,3}, SZ)
    for j in 1:SZ[2], i in 1:SZ[1]
        m[i,j,1    ] = 0.0
        m[i,j,SZ[3]] = 0.0
    end

    for k in 1:SZ[3], i in 1:SZ[1]
        m[i,1    ,k] = 0.0
        m[i,SZ[2],k] = 0.0
    end

    for k in 1:SZ[3], j in 1:SZ[2]
        m[1    ,j,k] = 0.0
        m[SZ[1],j,k] = 0.0
    end
end

#=
@brief 熱伝導率の設定 例題1 一様分布
@param [in]     λ      熱伝導率
@param [in]     SZ     配列長
=#
function setMat_1!(λ::Array{Float64,3}, SZ)
    λ .= 1.0
end

#=
@brief RHSの設定
@param [in]     b    RHSベクトル
@param [in]     SZ   配列長
=#
function calRHS!(b::Array{Float64,3}, SZ)
    b .= 0.0
end

#=
@param [in] SZ       内部セル数
@param [in] Δh       セル幅
@param [in] exact    厳密解
@param [in] θ        解ベクトル
@param [in] solver   ["jacobi", "sor", "pbicgstab"]
@param [in] smoother ["jacobi", "gs", ""]
@param [in] mode     動作モード   
=#
function main(SZ, Δh, exact, θ, solver, smoother, mode)
    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0 # default

    Cartesian.boundary_condition!(θ, SZ, Δh) # NonUniform 

    b = zeros(Float64, SZ[1], SZ[2], SZ[3])
    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    setMask!(mask, SZ)
    wk = zeros(Float64, SZ[1], SZ[2], SZ[3])

    if solver=="pbicgstab"
        pcg_p  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_p_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r0 = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_q  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_t_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end

    F = open("res.txt", "w")

    if solver=="sor"
        Cartesian.solveSOR!(θ, SZ, λ, b, mask, Δh, ω, F)
    elseif solver=="jacobi"
        Cartesian.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, F)
    elseif solver=="pbicgstab"
        Cartesian.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, Δh, SZ, ItrMax, smoother, F)
    else
        println("solver error")
    end
    nrm = sqrt(norm(θ-exact,2))
    @printf(F, "Sum of norm = %24.14E\n", nrm)
    close(F)
end

#=
@param [in] mode (
                    1--Uniform in Z-dir, (Cartesian)
                    2--Uniform for X&Y and Non-uniform in Z, but data is uniform
                    3--Uniform for X&Y and Non-uniform in Z
@param NXY  Number of inner cells for X&Y dir.
@param NZ   Number of inner cells for Z dir.
@param [in] solver    ["jacobi", "sor", "pbicgstab"]
@param [in] smoother  ["jacobi", "gs", ""]
=#
function q3d_fixed(mode::Int64, NXY::Int64, NZ::Int64, solver::String="sor", smoother::String="")
    MX = MY = NXY + 2  # Number of CVs including boundaries
    MZ = NZ + 2

    if !(mode==1 ||mode==2 || mode==3)
        println("mode error")
        return
    end

    if NXY != NZ
        println("NXY must be equal to NZ")
        return
    end
    
    dh::Float64 = 1.0 / NXY
    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0) #原点を仮定
    
    println("=== 計算条件 ===")
    println("SZ = ", SZ)
    println("Δh = ", Δh)
    println("Mode = ", mode)

    Z = zeros(Float64, SZ[3])
    ΔZ = zeros(Float64, SZ[3]-1)

    if mode >= 2
        # NonUniform格子の読み込み
        read_coord, numNodes = NonUniform.read_grid_file()
        if numNodes != Int(NZ)
            println("Number of generated grid is not match with parameter NZ")
            return
        end
        
        for k in 1:NZ
            Z[k+1] = read_coord[k]
        end
        Z[1] = 2*Z[2] - Z[3]
        Z[MZ] = 2*Z[MZ-1] - Z[MZ-2]
        
        for k in 2:MZ-1
            ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])
        end
        ΔZ[2] = 0.5*ΔZ[2]
        ΔZ[MZ-1] = 0.5*ΔZ[MZ-1]
        
        println("=== NonUniform格子情報 ===")
        @printf("Z座標範囲: [%.6f, %.6f]\n", Z[1], Z[MZ])
        @printf("最小格子間隔: %.6f\n", minimum(diff(Z)))
        @printf("最大格子間隔: %.6f\n", maximum(diff(Z)))
        
        # 格子分布の可視化
        plot_grid_distribution_nu(SZ, ox, Δh, Z, "grid_distribution_nu.png")
    end

    # 厳密解の計算
    exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    if mode == 1
        Cartesian.exact_solution!(exact, SZ, ox, Δh)
        println("=== Cartesian厳密解 ===")
        plot_slice(exact, SZ, "exact_cartesian.png")
    else
        NonUniform.exact_solution!(exact, SZ, ox, Δh, Z)
        println("=== NonUniform厳密解 ===")
        
        # 従来の方法（問題のあるプロット）
        plot_slice(exact, SZ, "exact_nu_org.png")
        println("従来プロット保存: exact_nu_org.png（格子座標系）")
        
        # 修正された方法（物理座標系プロット）
        plot_slice_nu(exact, SZ, ox, Δh, Z, "exact_nu_fixed.png")
        println("修正プロット保存: exact_nu_fixed.png（物理座標系）")
        
        # Cartesian版との比較用（参考）
        exact_cart = zeros(Float64, SZ[1], SZ[2], SZ[3])
        Cartesian.exact_solution!(exact_cart, SZ, ox, Δh)
        plot_slice(exact_cart, SZ, "exact_cartesian_ref.png")
        println("Cartesian参考プロット保存: exact_cartesian_ref.png")
        
        # 比較プロット作成
        plot_comparison_nu_cart(exact, exact_cart, SZ, ox, Δh, Z, "comparison_nu_cart.png")
        println("比較プロット保存: comparison_nu_cart.png")
    end
    
    return  # 計算部分は省略して可視化のみ
    
    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    @time main(SZ, Δh, exact, θ, solver, smoother, mode)

    if mode >= 2
        plot_slice2_nu(θ, SZ, ox, Δh, Z, "p_nu_fixed.png")
        plot_slice2_nu(θ-exact, SZ, ox, Δh, Z, "diff_nu_fixed.png")
    else
        plot_slice2(θ, SZ, "p.png")
        plot_slice2(θ-exact, SZ, "diff.png")
    end
end

# テスト実行
println("=== NonUniform格子可視化修正テスト ===")
q3d_fixed(2, 25, 25, "sor")  # NonUniform格子での実行