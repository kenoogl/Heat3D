module Cartesian
export exact_solution!, boundary_condition!, solveSOR!, solveJACOBI!, PBiCGSTAB!

using Printf

const ItrMax = 8000
const tol    = 1.0e-8
const FloatMin = 1.0e-37

# Harmonic mean
# @param a left value
# @param b right value
# @param ma mask for left
# @param mb mask for right
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))

#=
@brief 厳密解
@param [in]     e      解ベクトル
@param [in]     SZ     配列長
@param [in]     ox     原点座標
@param [in]     Δh     セル幅
=#
function exact_solution!(e::Array{Float64,3}, SZ, ox, Δh)
    r2 = sqrt(2.0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        z = ox[3] + Δh[3]*(k-1.5)
        e[i,j,k] = sin(π*x)*sin(π*y) / sinh(π*r2) * ( sinh(r2*π*z)-sinh(π*r2*(z-1.0)) )
    end
end


#=
@brief 境界条件
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     Δh   セル幅
=#
function boundary_condition!(p::Array{Float64,3}, SZ, ox, Δh)

    for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        a = sin(π*x)*sin(π*y)
        p[i,j,1    ] = a
        p[i,j,SZ[3]] = a
    end

    for k in 2:SZ[3]-1, i in 2:SZ[1]-1
        p[i,1    ,k] = 0.0
        p[i,SZ[2],k] = 0.0
    end

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        p[1    ,j,k] = 0.0
        p[SZ[1],j,k] = 0.0
    end
end

#=
@brief SOR法による求解
@param [in/out] θ    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    RHSベクトル
@param [in]     mask マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     F    ファイルディスクリプタ
=#
function solveSOR!(θ, SZ, λ, b, mask, Δh, ω, F)

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
@param [in]     ω    緩和係数
@param [in]     F    ファイルディスクリプタ
=#
function solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, F)

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
@brief SOR法の残差
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    加速係数
@ret                 1セルあたりの残差RMS
=#
function resSOR(p::Array{Float64,3}, 
                SZ, 
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64)

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = (dd + ω*(aw+as+ab))*dp / ω
        res += r*r
    end
    
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end

#=
@brief 緩和Jacobi法の残差
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    緩和係数
@ret                 1セルあたりの残差RMS
=#
function resJCB(p::Array{Float64,3},
                SZ,
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64)

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = dd*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end

#=
@brief 緩和Jacobi法
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    緩和係数
@param [out]    wk   ワーク用配列
@ret                 1セルあたりの残差RMS
=#
function jacobi!(p::Array{Float64,3},
                 SZ,
                 λ::Array{Float64,3}, 
                 b::Array{Float64,3},
                 m::Array{Float64,3}, 
                 Δh, 
                 ω::Float64,
                wk::Array{Float64,3})

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        pn = pp + ω * dp
        wk[i,j,k] = pn
        r = dd*dp / ω
        res += r*r
    end

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        p[i,j,k] = wk[i,j,k]
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


#=
@brief SOR法
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    加速係数
@ret                 1セルあたりの残差RMS
=#
function sor!(p::Array{Float64,3}, 
              SZ, 
              λ::Array{Float64,3}, 
              b::Array{Float64,3},
              m::Array{Float64,3}, 
              Δh, 
              ω::Float64)

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
        + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
        + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        pn = pp + ω * dp
        p[i,j,k] = pn
        r = (dd + ω*(aw+as+ab))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


#=
@brief RB-SOR法のカーネル
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    加速係数
@param [in] 　　　　　　　　color R or B
@ret                 残差2乗和
=#
function rbsor_core!(p::Array{Float64,3}, 
                     SZ, 
                     λ::Array{Float64,3}, 
                     b::Array{Float64,3},
                     m::Array{Float64,3}, 
                     Δh, 
                     ω::Float64, 
                     color::Int)

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            pp = p[i,j,k]
            bb = b[i,j,k]
            λ0 = λ[i,j,k]
            m0 = m[i,j,k]
            ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
            aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
            an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
            as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
            at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
            ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
            dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
            ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
               + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
               + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
            dp = (((ss-bb)/dd - pp)) * m0
            pn = pp + ω * dp
            p[i,j,k] = pn
            r = (dd + ω*(aw+as+ab))*dp / ω
            res += r*r
        end
    end

    return res
end


#=
@brief RB-SOR法
@param [in,out] θ    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     mask マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    加速係数
@ret                 セルあたりの残差RMS
=#
function rbsor!(θ::Array{Float64,3}, 
                SZ, 
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
             mask::Array{Float64,3}, 
                Δh, 
                ω::Float64)

    res::Float64 = 0.0

    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        r = rbsor_core!(θ, SZ, λ, b, mask, Δh, ω, c)
        res += r
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


#=
@brief PBiCGSTAB反復
@param [in,out] X      解ベクトル
@param [in]     B      RHSベクトル
@param [in]     pcg*   ワークベクトル
@param [in]     λ      熱伝導率
@param [in]     mask   マスク配列
@param [in]     Δh     セル幅
@param [in]     SZ     配列長
@param [in] 　　　　　　　　ω      加速係数
@param [in] 　　　　　　　　ItrMax 前処理反復数
@param [in] 　　　　　　　　smoother ["gs", "jacobi", ""]
@param [in]     F      ファイルディスクリプタ
@ret                   セルあたりの残差RMS
=#
function PBiCGSTAB!(X::Array{Float64,3}, 
                    B::Array{Float64,3},
                pcg_q::Array{Float64,3},
                pcg_r::Array{Float64,3},
               pcg_r0::Array{Float64,3},
                pcg_p::Array{Float64,3},
               pcg_p_::Array{Float64,3},
                pcg_s::Array{Float64,3},
               pcg_s_::Array{Float64,3},
               pcg_t_::Array{Float64,3},
                    λ::Array{Float64,3},
                 mask::Array{Float64,3},
                   wk::Array{Float64,3},
                    Δh,
                    SZ,
               ItrMax::Int,
             smoother::String,
                    F)
   
    fill!(pcg_q, 0.0)
    res0 = CalcRK!(SZ, pcg_r, X, B, λ, mask, Δh)
    println("Inital residual = ", res0)
    copy!(pcg_r0, pcg_r)

    rho_old::Float64 = 1.0
    alpha::Float64 = 0.0
    omega::Float64  = 1.0
    r_omega::Float64 = -omega
    beta::Float64 = 0.0

    for itr in 1:ItrMax
        rho = Fdot2(pcg_r, pcg_r0, SZ) # 非計算部分はゼロのこと

        if abs(rho) < FloatMin
            itr = 0
            break
        end

        if itr == 1
            copy!(pcg_p, pcg_r)
        else
            beta = rho / rho_old * alpha / omega
            BiCG1!(pcg_p, pcg_r, pcg_q, beta, omega, SZ)
        end

        fill!(pcg_p_, 0.0)
        Preconditioner!(pcg_p_, pcg_p, λ, mask, wk, SZ, Δh, smoother)

        CalcAX!(pcg_q, pcg_p_, SZ, Δh, λ, mask)
        alpha = rho / Fdot2(pcg_q, pcg_r0, SZ)
        r_alpha = -alpha
        Triad!(pcg_s, pcg_q, pcg_r, r_alpha, SZ)

        fill!(pcg_s_, 0.0)
        Preconditioner!(pcg_s_, pcg_s, λ, mask, wk, SZ, Δh, smoother);

        CalcAX!(pcg_t_, pcg_s_, SZ, Δh, λ, mask)
        omega = Fdot2(pcg_t_, pcg_s, SZ) / Fdot1(pcg_t_, SZ)
        r_omega = -omega

        BICG2!(X, pcg_p_, pcg_s_, alpha , omega, SZ)
        boundary_condition!(X, SZ, Δh)

        Triad!(pcg_r, pcg_t_, pcg_s, r_omega, SZ)
        res = sqrt(Fdot1(pcg_r, SZ))/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
        res /= res0
        #println(itr, " ", res)
        @printf(F, "%10d %24.14E\n", itr, res) # 時間計測の場合にはコメントアウト

        if res<tol
            println("Converged at ", itr)
            break
        end

        rho_old = rho
    end # itr
end


#=
@brief 残差ベクトルの計算
@param [in]     SZ   配列寸法
@param [out]    r    残差ベクトル
@param [in]     p    解ベクトル
@param [in]     b    右辺ベクトル
@param [in]     λ    係数
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@ret                 セルあたりの残差RMS
=#
function CalcRK!(SZ, 
                r::Array{Float64,3}, 
                p::Array{Float64,3}, 
                b::Array{Float64,3},
                λ::Array{Float64,3}, 
                m::Array{Float64,3}, 
                Δh)

    res::Float64 = 0.0
    dx0::Float64 = Δh[1]
    dy0::Float64 = Δh[2]
    dz0::Float64 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        rs = (b[i,j,k] - (ss - dd * p[i,j,k]))* m0
        r[i,j,k] = rs
        res += rs*rs
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


#=
@brief ベクトルの内積
@param [in]     x    ベクトル
@param [in]     SZ   配列長
@ret            　    内積
=#
function Fdot1(x::Array{Float64,3}, SZ)

    y::Float64 = 0.0
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        y += x[i,j,k] * x[i,j,k]
    end
    return y
end



#=
@brief 2ベクトルの内積
@param [in]     x    ベクトル
@param [in]     y    ベクトル
@param [in]     SZ   配列長
@ret            　　   内積
=#
function Fdot2(x::Array{Float64,3}, y::Array{Float64,3}, SZ)

    xy::Float64 = 0.0
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        xy += x[i,j,k] * y[i,j,k]
    end
    return xy
end



#=
@brief BiCGstabの部分演算1
@param [in,out] p    ベクトル
@param [in]     r    ベクトル
@param [in]     q    ベクトル
@param [in]     beta 係数
@param [in]     omg  係数
@param [in]     SZ   配列長
=#
function BiCG1!(p::Array{Float64,3}, 
                r::Array{Float64,3}, 
                q::Array{Float64,3}, 
             beta::Float64, 
              omg::Float64, 
                SZ)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        p[i,j,k] = r[i,j,k] + beta * (p[i,j,k] - omg * q[i,j,k])
    end
end


#=
@brief 前処理
@param [in,out] xx   解ベクトル
@param [in]     bb   RHSベクトル
@param [in]     λ    熱伝導率
@param [in]     mask マスク配列
@param [in]     wk   作業ベクトル
@param [in]     SZ   配列長
@param [in]     Δh   セル幅
@param [in]     smoother  ["jacobi", "dor", ""]
=#
function Preconditioner!(xx::Array{Float64,3},
                         bb::Array{Float64,3}, 
                          λ::Array{Float64,3}, 
                       mask::Array{Float64,3}, 
                         wk::Array{Float64,3}, 
                         SZ,
                         Δh,
                   smoother::String)

    res::Float64 = 0.0
    LCmax::Int = 5

    if smoother=="jacobi"
        #P_Jacobi!(xx, bb, LCmax, SZ, λ, mask, Δh, wk)
        for _ in 1:LCmax
            res = jacobi!(xx, SZ, λ, bb, mask, Δh, 0.8, wk)
        end
    elseif smoother=="gs"
        #P_SOR!(xx, bb, LCmax, SZ, λ, mask, Δh)
        for _ in 1:LCmax
            res = sor!(xx, SZ, λ, bb, mask, Δh, 1.0)
        end
    else
        copy!(xx, bb)
    end
end



#=
@brief AXの計算
@param [out] ap   AX
@param [in]  p    解ベクトル
@param [in]  SZ   配列長
@param [in]  Δh   セル幅
@param [in]  λ    熱伝導率
@param [in]  mask マスク配列
=#
function CalcAX!(ap::Array{Float64,3}, 
                  p::Array{Float64,3}, 
                  SZ,
                  Δh,
                  λ::Array{Float64,3}, 
                  m::Array{Float64,3})

    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        ap[i,j,k] = (ss - dd*p[i,j,k]) * m0
    end
end


#=
@brief AXPYZ
@param [out]    z    ベクトル
@param [in]     y    ベクトル
@param [in]     x    ベクトル
@param [in]     a    係数
@param [in]     SZ   配列長
=#
function Triad!(z::Array{Float64,3}, 
                x::Array{Float64,3}, 
                y::Array{Float64,3}, 
                a::Float64, 
                SZ)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] = a * x[i,j,k] + y[i,j,k]
    end
end


#=
!! @brief BiCGstab 2
!! @param [in,out] z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     b    係数
!! @param [in]     SZ   配列長
=#
function BICG2!(z::Array{Float64,3}, 
                x::Array{Float64,3}, 
                y::Array{Float64,3}, 
                a::Float64, 
                b::Float64, 
                SZ)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] += a * x[i,j,k] + b * y[i,j,k]
    end
end

end # end of module Cartesian