using Printf
using Random
using LinearAlgebra

include("heat3d_Cartesian.jl")
include("heat3d_NonUniform.jl")

const ω      = 1.0


#include("../base/pbicgstab.jl")
#
include("plotter.jl")

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
@brief 熱伝導率の設定 例題2 5分割
@param [in]     λ      熱伝導率
@param [in]     SZ     配列長

function setMat_2!(λ::Array{Float64,3}, SZ)

    d = div(SZ[1]-2, 5)
    i1 = 2 + d
    i2 = 2 + d*2
    i3 = 2 + d*3
    i4 = 2 + d*4

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        for i in 2:i1-1
            λ[i,j,k] = 1.0
        end
        for i in i1:i2-1
            λ[i,j,k] = 5.0
        end
        for i in i2:i3-1
            λ[i,j,k] = 10.0
        end
        for i in i3:i4-1
            λ[i,j,k] = 20.0
        end
        for i in i4:SZ[1]-1
            λ[i,j,k] = 30.0
        end
    end
end


#=
@brief 熱伝導率の設定 例題3 ランダム分布
@param [in]     λ      熱伝導率
@param [in]     SZ     配列長
=#
function setMat_3!(λ::Array{Float64,3}, SZ)
    copy!(λ, rand(SZ[1], SZ[2], SZ[3]) .* 50.0)
end
=#


#=
@brief RHSの設定
@param [in]     b    RHSベクトル
@param [in]     SZ   配列長
=#
function calRHS!(b::Array{Float64,3}, SZ)
    b .= 0.0
end



#=
@param [in] mode     動作モード   
@param [in] SZ       内部セル数
@param [in] Δh       セル幅
@param [in] exact    厳密解
@param [in] θ        解ベクトル
@param [in] solver   ["jacobi", "sor", "pbicgstab"]
@param [in] smoother ["jacobi", "gs", ""]
=#
function main(mode, SZ, ox, Δh, exact, θ, Z, ΔZ, solver, smoother)

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0 # default
    #=
    if example==1
        setMat_1!(λ, SZ)
    elseif example==2
        setMat_2!(λ, SZ)
    elseif example==3
        setMat_3!(λ, SZ)
    else
        println("example error")
        return
    end
    plot_slice2(λ, SZ, "lambda.png")
    =#
    if mode==1
        Cartesian.boundary_condition!(θ, SZ, ox, Δh)
    else
        NonUniform.boundary_condition!(θ, SZ, ox, Δh)
    end
    #plot_slice2(θ, SZ, "p0.png")

    b = zeros(Float64, SZ[1], SZ[2], SZ[3])
    #calRHS!(b, SZ)

    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    setMask!(mask, SZ)
    #plot_slice2(mask, SZ, "mask.png")

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
        if mode==1
            Cartesian.solveSOR!(θ, SZ, λ, b, mask, Δh, ω, F)
        else
            NonUniform.solveSOR!(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, 3, SZ[3]-2, F)
        end
    elseif solver=="jacobi"
        if mode==1
            Cartesian.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, F)
        else
            NonUniform.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, Z, ΔZ, 3, SZ[3]-2, F)
        end
    elseif solver=="pbicgstab"
        if mode==1
            Cartesian.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, Δh, SZ, ItrMax, smoother, F)
        else
        end
    else
        println("solver error")
    end


    nrm = sqrt(norm(θ-exact,2))
    @printf(F, "Sum of norm = %24.14E\n", nrm)
    close(F)
end

#=
# Visio用のsphフォーマット出力
function writeSPH(size, org, pch, step, time, var)
    ccall((:write_sph_d, "./iosph.so"), Nothing,
    (Ref{Int32},    # size
     Ref{Float64},  # org
     Ref{Float64},  # pch
     Ref{Int32},    # step
     Ref{Float64},  # time
     Ref{Float64}), # var
     size, org, pch, step, time, var)
end
=#

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
function q3d(mode::Int, NXY::Int, NZ::Int, solver::String="sor", smoother::String="")

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
    @show typeof(mode)
    println(SZ)
    println(Δh)

    Z = zeros(Float64, SZ[3])
    ΔZ= zeros(Float64, SZ[3]-1)
    #@show typeof(Z)

    if mode>=2
        read_coord, numNodes = NonUniform.read_grid_file()
        if numNodes != NZ
            println("Number of genereted grid is not match with parameter NZ")
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
    end

    exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    if mode==1
        Cartesian.exact_solution!(exact, SZ, ox, Δh)
        plot_slice(exact, SZ, ox, Δh,"exact.png")
    else
        NonUniform.exact_solution!(exact, SZ, ox, Δh, Z)
        plot_slice_nu(exact, 0.5, SZ, ox, Δh, Z, "exact_nu.png")
    end
    
    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])

    @time main(mode, SZ, ox, Δh, exact, θ, Z, ΔZ, solver, smoother)

    if mode==1
        plot_slice2(θ, SZ, ox, Δh, "p.png", "solution")
        plot_slice2(θ-exact, SZ, ox, Δh, "diff.png", "diff")
    else
        plot_slice2_nu(θ, 0.5, SZ, ox, Δh, Z, "p_nu.png", "solution")
        plot_slice2_nu(θ-exact, 0.5, SZ, ox, Δh, Z, "diff_nu.png", "diff")
    end
end

q3d(2, 1, 1) # just compile　JITコンパイルを行うためにパラメータは1
q3d(2, 25, 25, "jacobi") # ここで本実行し、計測