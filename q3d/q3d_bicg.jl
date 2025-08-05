using Printf
using Plots

const ItrMax = 2000
const ω      = 1.0
const tol    = 1.0e-8

include("pbicgstab.jl")
include("modelA.jl")


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
@brief XZ断面（内部セル）
@param [in]     ｄ      解ベクトル
@param [in]     SZ     配列長
@param [in]     fname  ファイル名
=#
function plot_slice(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[2],2)
    s = d[2:SZ[1]-1,j,2:SZ[3]-1]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end


#=
@brief XZ断面(全セル)
@param [in]     ｄ      解ベクトル
@param [in]     SZ     配列長
@param [in]     fname  ファイル名
=#
function plot_slice2(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[2],2)
    s = d[1:SZ[1],j,1:SZ[3]]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
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
@brief SOR法による求解
@param [in/out] θ    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    RHSベクトル
@param [in]     mask マスク配列
@param [in]     Δh   セル幅
@param [in]     F    ファイルディスクリプタ
=#
function solveSOR!(θ, SZ, λ, b, mask, Δh, F)

    res0 = resSOR(θ, SZ, λ, b, mask, Δh, ω)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        res = sor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        #res = rbsor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        #println(n, " ", res)
        @printf(F, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end


#=
@brief 緩和ヤコビ法による求解
@param [in/out] θ    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    RHSベクトル
@param [in]     mask マスク配列
@param [in]     wk   ワーク配列
@param [in]     Δh   セル幅
@param [in]     F    ファイルディスクリプタ
=#
function solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, F)

    res0 = resJCB(θ, SZ, λ, b, mask, Δh, ω)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        res = jacobi!(θ, SZ, λ, b, mask, Δh, ω, wk) / res0
        #res = rbsor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        
        @printf(F, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end


#=
@param [in] SZ       内部セル数
@param [in] Δh       セル幅
@param [in] exact    厳密解
@param [in] θ        解ベクトル
@param [in] solver   ["jacobi", "sor", "pbicgstab"]
@param [in] smoother ["jacobi", "gs", ""]
=#
function main(SZ, Δh, exact, θ, solver, smoother)

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0 # default
 
    plot_slice2(λ, SZ, "lambda.png")

    boundary_condition!(θ, SZ, Δh)
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
        solveSOR!(θ, SZ, λ, b, mask, Δh, F)
    elseif solver=="jacobi"
        solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, F)
    elseif solver=="pbicgstab"
        PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, Δh, SZ, ItrMax, smoother, F)
    else
        println("solver error")
    end
    close(F)
end

#=
# Visio用のsphフォーマット出力
function writeSPH(size, org, pch, step, time, var)
    ccall((:write_sph_d, "./iosph.so"), Nothing,
    (Ref{Int64},    # size
     Ref{Float64},  # org
     Ref{Float64},  # pch
     Ref{Int64},    # step
     Ref{Float64},  # time
     Ref{Float64}), # var
     size, org, pch, step, time, var)
end
=#

#=
@param mode (1--Uniform in Z-dir, 2--Non-uniform in Z), while Uniform for X&Y
@param NXY  Number of inner cells for X&Y dir.
@param NZ   Number of inner cells for Z dir.
@param [in] solver    ["jacobi", "sor", "pbicgstab"]
@param [in] smoother  ["jacobi", "gs", ""]
=#
function q3d(mode::Int64, NXY::Int64, NZ::Int64=14)
    solver = "sor"
    smoother=""
    
    MX = MY = NXY + 2  # Number of CVs including boundaries
    if mode==2
        NZ = 14
    end
    MZ = NZ + 2

    dx::Float64 = 1.2 / NXY
    dy::Float64 = 1.2 / NXY
    dz::Float64 = zm5 / NZ
    Δh = (dx, dy, dz) 
    ox = (0.0, 0.0, 0.0)
    SZ = (MX, MY, MZ)
    println(SZ)
    println(Δh)
    
    ID = Array{Int64}(undef, SZ[1], SZ[2], SZ[3])
    Z = zeros(Float64, SZ[3]+1)

    coordZ!(Z, SZ, mode, ox, Δh[3])
    modelA!(ID, Z, SZ, ox, Δh)

    λ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    setLambda!(λ, SZ, ID)

    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])

    @time main(SZ, Δh, exact, θ, solver, smoother)

    plot_slice2(θ, SZ, "p.png")
    plot_slice2(θ-exact, SZ, "diff.png")
end

q3d(1, 1, 1) # just compile　JITコンパイルを行うためにパラメータは1
q3d(1, 60, 20) # ここで本実行し、計測