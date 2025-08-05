using Printf
using Plots
using Random

const ItrMax = 8000
const ω      = 1.0
const tol    = 1.0e-8

include("pbicgstab.jl")


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
=#
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


#=
@brief 厳密解
@param [in]     e      解ベクトル
@param [in]     SZ     配列長
@param [in]     Δh     セル幅
=#
function exact_solution!(e::Array{Float64,3}, SZ, Δh)
    r2 = sqrt(2.0)
    ox = (0.0, 0.0, 0.0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        z = ox[3] + Δh[3]*(k-1.5)
        e[i,j,k] = sin(π*x)*sin(π*y) / sinh(π*r2) * ( sinh(r2*π*z)-sinh(π*r2*(z-1.0)) )
    end
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
@param [in] example   例題番号
=#
function main(SZ, Δh, exact, θ, solver, smoother, example)

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0 # default
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
@param [in] NN_inner  内部セル数
@param [in] example   例題番号
@param [in] solver    ["jacobi", "sor", "pbicgstab"]
@param [in] smoother  ["jacobi", "gs", ""]
=#
function diff3d(NN_inner::Int64, example=1, solver::String="sor", smoother::String="")
    NN = NN_inner + 2  # Number of CVs including boundaries

    dh::Float64 = 1.0 / NN_inner
    Δh = (dh, dh, dh)
    # ox = (0.0, 0.0, 0.0) 原点を仮定
    SZ = (NN, NN, NN)
    

    exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    exact_solution!(exact, SZ, Δh)
    plot_slice(exact, SZ, "exact.png")
    #writeSPH(SZ, ox, Δh, 0, 0.0, exact)

    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])

    @time main(SZ, Δh, exact, θ, solver, smoother, example)

    plot_slice2(θ, SZ, "p.png")
    plot_slice2(θ-exact, SZ, "diff.png")
end

diff3d(1) # just compile　JITコンパイルを行うためにパラメータは1
diff3d(25, 2, "jacobi", "") # ここで本実行し、計測
# solver : ["jacobi", "sor", "pbicgstab"]
# smoother : ["jacobi", "gs", ""]