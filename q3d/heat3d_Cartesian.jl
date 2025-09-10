module Cartesian
export exact_solution!, solveSOR!, PBiCGSTAB!, CG!

using Printf

include("const.jl")

"""
Harmonic mean
@param a left value
@param b right value
@param ma mask for left
@param mb mask for right
"""
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))

"""
@brief 厳密解
@param [in]     e      解ベクトル
@param [in]     ox     原点座標
@param [in]     Δh     セル幅
"""
function exact_solution!(e::Array{Float64,3}, ox, Δh)
    SZ = size(e)
    r2 = sqrt(2.0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        z = ox[3] + Δh[3]*(k-1.5)
        e[i,j,k] = sin(π*x)*sin(π*y) / sinh(π*r2) * ( sinh(r2*π*z)-sinh(π*r2*(z-1.0)) )
    end
end


"""
@brief SOR法による求解
@param [in/out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    RHSベクトル
@param [in]     mask マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     F    ファイルディスクリプタ
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
"""
function solveSOR!(θ, λ, b, mask, Δh, ω, F, tol, HF::Vector{Float64}, HT::Vector{Float64})

    res0 = resSOR(θ, λ, b, mask, Δh, ω, HF, HT)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:Constant.ItrMax
        #res = sor!(θ, λ, b, mask, Δh, ω, HF, HT) / res0
        res = rbsor!(θ, λ, b, mask, Δh, ω, HF, HT) / res0
        #println(n, " ", res)
        @printf(F, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        @printf(stdout, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end


"""
@brief SOR法の残差
@param [in,out] p    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                 1セルあたりの残差RMS
"""
function resSOR(p::Array{Float64,3}, 
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64,
                HF::Vector{Float64},
                HT::Vector{Float64})
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    dz1 = 1.0 / dz0

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        me = 1.0-m[i+1,j  ,k  ]
        mw = 1.0-m[i-1,j  ,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mt = 1.0-m[i  ,j  ,k+1]
        mb = 1.0-m[i  ,j  ,k-1]
        bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +(mb*HF[5] - mt*HF[6])*dz1 ) # 熱流束境界
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2 + mb*dz1*HT[5]
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2 + mt*dz1*HT[6]
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = (dd + ω*(axm+aym+azm))*dp / ω
        res += r*r
    end
    
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


"""
@brief SOR法
@param [in,out] p    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                 1セルあたりの残差RMS
"""
function sor!(p::Array{Float64,3}, 
              λ::Array{Float64,3}, 
              b::Array{Float64,3},
              m::Array{Float64,3}, 
              Δh, 
              ω::Float64,
              HF::Vector{Float64},
              HT::Vector{Float64})
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    dz1 = 1.0 / dz0

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        me = 1.0-m[i+1,j  ,k  ]
        mw = 1.0-m[i-1,j  ,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mt = 1.0-m[i  ,j  ,k+1]
        mb = 1.0-m[i  ,j  ,k-1]
        bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +(mb*HF[5] - mt*HF[6])*dz1 ) # 熱流束境界
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2 + mb*dz1*HT[5]
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2 + mt*dz1*HT[6]
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
        + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
        + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        pn = pp + ω * dp
        p[i,j,k] = pn
        r = (dd + ω*(axm+aym+azm))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


"""
@brief RB-SOR法のカーネル
@param [in,out] p    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@param [in]     color R or B
@ret                 残差2乗和
"""
function rbsor_core!(p::Array{Float64,3}, 
                     λ::Array{Float64,3}, 
                     b::Array{Float64,3},
                     m::Array{Float64,3}, 
                     Δh, 
                     ω::Float64,
                     HF::Vector{Float64},
                     HT::Vector{Float64},
                     color::Int)
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    dz1 = 1.0 / dz0

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            pp = p[i,j,k]
            λ0 = λ[i,j,k]
            m0 = m[i,j,k]
            me = 1.0-m[i+1,j  ,k  ]
            mw = 1.0-m[i-1,j  ,k  ]
            mn = 1.0-m[i  ,j+1,k  ]
            ms = 1.0-m[i  ,j-1,k  ]
            mt = 1.0-m[i  ,j  ,k+1]
            mb = 1.0-m[i  ,j  ,k-1]
            bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +(mb*HF[5] - mt*HF[6])*dz1 ) # 熱流束境界
            axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
            axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
            aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
            ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
            azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2 + mb*dz1*HT[5]
            azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2 + mt*dz1*HT[6]
            dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
            ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
               + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
               + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
            dp = (((ss-bb)/dd - pp)) * m0
            pn = pp + ω * dp
            p[i,j,k] = pn
            r = (dd + ω*(axm+aym+azm))*dp / ω
            res += r*r
        end
    end

    return res
end


"""
@brief RB-SOR法
@param [in,out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     mask マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                 セルあたりの残差RMS
"""
function rbsor!(θ::Array{Float64,3}, 
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
             mask::Array{Float64,3}, 
                Δh, 
                ω::Float64,
                HF::Vector{Float64},
                HT::Vector{Float64})
    SZ = size(b)
    res::Float64 = 0.0

    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        r = rbsor_core!(θ, λ, b, mask, Δh, ω, HF, HT, c)
        res += r
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


"""
@brief PBiCGSTAB反復
@param [in,out] X      解ベクトル
@param [in]     B      RHSベクトル
@param [in]     pcg*   ワークベクトル
@param [in]     λ      熱伝導率
@param [in]     mask   マスク配列
@param [in]     Δh     セル幅
@param [in]     ω      加速係数
@param [in]     ItrMax 前処理反復数
@param [in]            smoother ["gs", "jacobi", ""]
@param [in]     F      ファイルディスクリプタ
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                   セルあたりの残差RMS
"""
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
                   ox,
                    Δh,
             smoother::String,
                    F,
                    mode,
                    tol,
                    HF::Vector{Float64},
                    HT::Vector{Float64})
    SZ = size(X)
    pcg_q .= 0.0  #fill!(pcg_q, 0.0)
    res0 = CalcRK!(pcg_r, X, B, λ, mask, Δh, HF, HT)
    println("Inital residual = ", res0)
    pcg_r0 .= pcg_r  #copy!(pcg_r0, pcg_r)

    rho_old::Float64 = 1.0
    alpha::Float64 = 0.0
    omega::Float64  = 1.0
    r_omega::Float64 = -omega
    beta::Float64 = 0.0

    for itr in 1:Constant.ItrMax
        rho = Fdot2(pcg_r, pcg_r0) # 非計算部分はゼロのこと

        if abs(rho) < Constant.FloatMin
            itr = 0
            break
        end

        if itr == 1
            pcg_p .= pcg_r  #copy!(pcg_p, pcg_r)
        else
            beta = rho / rho_old * alpha / omega
            BiCG1!(pcg_p, pcg_r, pcg_q, beta, omega)
        end

        pcg_p_ .= 0.0  #fill!(pcg_p_, 0.0)
        Preconditioner!(pcg_p_, pcg_p, λ, mask, Δh, smoother, HF, HT)

        CalcAX!(pcg_q, pcg_p_, Δh, λ, mask, HF, HT)
        alpha = rho / Fdot2(pcg_q, pcg_r0)
        r_alpha = -alpha
        Triad!(pcg_s, pcg_q, pcg_r, r_alpha)

        pcg_s_ .= 0.0  #fill!(pcg_s_, 0.0)
        Preconditioner!(pcg_s_, pcg_s, λ, mask, Δh, smoother, HF, HT);

        CalcAX!(pcg_t_, pcg_s_, Δh, λ, mask, HF, HT)
        omega = Fdot2(pcg_t_, pcg_s) / Fdot1(pcg_t_)
        r_omega = -omega

        BICG2!(X, pcg_p_, pcg_s_, alpha , omega)

        #=
        if mode==4
            # boundary_condition4!(X, SZ)
        else
            boundary_condition!(X, ox, Δh)
        end
        =#

        Triad!(pcg_r, pcg_t_, pcg_s, r_omega)
        res = sqrt(Fdot1(pcg_r))/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
        res /= res0
        #println(itr, " ", res)
        @printf(F, "%10d %24.14E\n", itr, res) # 時間計測の場合にはコメントアウト
        @printf(stdout, "%10d %24.14E\n", itr, res) # 時間計測の場合にはコメントアウト

        if res<tol
            println("Converged at ", itr)
            break
        end

        rho_old = rho
    end # itr
    @printf(stdout, "\n")
end


"""
@brief 残差ベクトルの計算
@param [out]    r    残差ベクトル
@param [in]     p    解ベクトル
@param [in]     b    RHS
@param [in]     λ    係数
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                 セルあたりの残差RMS
"""
function CalcRK!(
                r::Array{Float64,3}, 
                p::Array{Float64,3}, 
                b::Array{Float64,3},
                λ::Array{Float64,3}, 
                m::Array{Float64,3}, 
                Δh,
                HF::Vector{Float64},
                HT::Vector{Float64})
    SZ = size(r)
    res::Float64 = 0.0
    dx0::Float64 = Δh[1]
    dy0::Float64 = Δh[2]
    dz0::Float64 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    dz1 = 1.0 / dz0

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        me = 1.0-m[i+1,j  ,k  ]
        mw = 1.0-m[i-1,j  ,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mt = 1.0-m[i  ,j  ,k+1]
        mb = 1.0-m[i  ,j  ,k-1]
        bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +(mb*HF[5] - mt*HF[6])*dz1 ) # 熱流束境界
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2 + mb*dz1*HT[5]
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2 + mt*dz1*HT[6]
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        rs = (bb - (ss - dd * p[i,j,k]))* m0
        r[i,j,k] = rs
        res += rs*rs
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


"""
@brief ベクトルの内積
@param [in]     x    ベクトル
@ret                 内積
"""
function Fdot1(x::Array{Float64,3})
    SZ = size(x)
    y::Float64 = 0.0
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        y += x[i,j,k] * x[i,j,k]
    end
    return y
end


"""
@brief 2ベクトルの内積
@param [in]     x    ベクトル
@param [in]     y    ベクトル
@ret                 内積
"""
function Fdot2(x::Array{Float64,3}, y::Array{Float64,3})
    SZ = size(x)
    xy::Float64 = 0.0
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        xy += x[i,j,k] * y[i,j,k]
    end
    return xy
end


"""
@brief BiCGstabの部分演算1
@param [in,out] p    ベクトル
@param [in]     r    ベクトル
@param [in]     q    ベクトル
@param [in]     beta 係数
@param [in]     omg  係数
"""
function BiCG1!(p::Array{Float64,3}, 
                r::Array{Float64,3}, 
                q::Array{Float64,3}, 
             beta::Float64, 
              omg::Float64)
    SZ = size(p)
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        p[i,j,k] = r[i,j,k] + beta * (p[i,j,k] - omg * q[i,j,k])
    end
end


"""
@brief 前処理
@param [in,out] xx   解ベクトル
@param [in]     bb   RHSベクトル
@param [in]     λ    熱伝導率
@param [in]     mask マスク配列
@param [in]     Δh   セル幅
@param [in]     smoother  ["jacobi", "gs", ""]
"""
function Preconditioner!(xx::Array{Float64,3},
                         bb::Array{Float64,3}, 
                          λ::Array{Float64,3}, 
                       mask::Array{Float64,3}, 
                         Δh,
                   smoother::String,
                   HF::Vector{Float64},
                   HT::Vector{Float64})

    res::Float64 = 0.0
    LCmax::Int = 5
    #=
    if smoother=="jacobi"
        #P_Jacobi!(xx, bb, LCmax, SZ, λ, mask, Δh, wk)
        for _ in 1:LCmax
            res = jacobi!(xx, λ, bb, mask, Δh, 0.8, wk)
        end
    else
    =#
    if smoother=="gs"
        #P_SOR!(xx, bb, LCmax, SZ, λ, mask, Δh)
        for _ in 1:LCmax
            res = rbsor!(xx, λ, bb, mask, Δh, 1.0, HF, HT)
        end
    else
        xx .= bb  #copy!(xx, bb)
    end
end



"""
@brief AXの計算
@param [out] ap   AX
@param [in]  p    解ベクトル
@param [in]  Δh   セル幅
@param [in]  λ    熱伝導率
@param [in]  mask マスク配列
@param [in]  HF   熱流束境界の値
@param [in]  HT   熱伝達境界の値
"""
function CalcAX!(ap::Array{Float64,3}, 
                  p::Array{Float64,3}, 
                  Δh,
                  λ::Array{Float64,3}, 
                  m::Array{Float64,3},
                  HF::Vector{Float64},
                  HT::Vector{Float64})
    SZ = size(p)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    dz1 = 1.0 / dz0

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        me = 1.0-m[i+1,j  ,k  ]
        mw = 1.0-m[i-1,j  ,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mt = 1.0-m[i  ,j  ,k+1]
        mb = 1.0-m[i  ,j  ,k-1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2 + mb*dz1*HT[5]
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2 + mt*dz1*HT[6]
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        ap[i,j,k] = (ss - dd*p[i,j,k]) * m0
    end
end


"""
@brief AXPYZ
@param [out]    z    ベクトル
@param [in]     y    ベクトル
@param [in]     x    ベクトル
@param [in]     a    係数
"""
function Triad!(z::Array{Float64,3}, 
                x::Array{Float64,3}, 
                y::Array{Float64,3}, 
                a::Float64)
    SZ = size(z)
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] = a * x[i,j,k] + y[i,j,k]
    end
end


"""
!! @brief BiCGstab 2
!! @param [in,out] z    ベクトル
!! @param [in]     y    ベクトル
!! @param [in]     x    ベクトル
!! @param [in]     a    係数
!! @param [in]     b    係数
"""
function BICG2!(z::Array{Float64,3}, 
                x::Array{Float64,3}, 
                y::Array{Float64,3}, 
                a::Float64, 
                b::Float64)
    SZ = size(z)
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] += a * x[i,j,k] + b * y[i,j,k]
    end
end

"""
@brief CG反復
@param [in,out] X      解ベクトル
@param [in]     B      RHSベクトル
@param [in]     p,r,ax,ap   ワークベクトル
@param [in]     λ      熱伝導率
@param [in]     mask   マスク配列
@param [in]     F      ファイルディスクリプタ
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                   セルあたりの残差RMS
"""
function CG!(X::Array{Float64,3}, 
             B::Array{Float64,3},
             p::Array{Float64,3},
             r::Array{Float64,3},
            ax::Array{Float64,3},
            ap::Array{Float64,3},
             λ::Array{Float64,3},
          mask::Array{Float64,3},
            Δh,
            F,
            tol,
            HF::Vector{Float64},
            HT::Vector{Float64})
    sd=false
    SZ = size(X)
    res0 = CalcRK!(r, X, B, λ, mask, Δh, HF, HT)
    p .= r

    println("Inital residual = ", res0)

    for itr in 1:Constant.ItrMax
        CalcAX!(ap, p, Δh, λ, mask, HF, HT)
        alpha = Fdot2(p, r) / Fdot2(p, ap)

        for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
            X[i,j,k] += alpha *  p[i,j,k]
            r[i,j,k] -= alpha * ap[i,j,k]
        end

        res = sqrt(Fdot1(r))/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
        res /= res0

        @printf(F, "%10d %24.14E\n", itr, res)
        @printf(stdout, "%10d %24.14E\n", itr, res) # 時間計測の場合にはコメントアウト

        if res<tol
            println("Converged at ", itr)
            break
        end

        if sd==false # CG method
            beta = -Fdot2(r, ap) / Fdot2(p, ap)
            p .= r .+ beta * p
        else
            p .= r # SD method
        end

    end # itr
    @printf(stdout, "\n")
end

end # end of module Cartesian

#=
"""
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
"""
function solveJACOBI!(θ, λ, b, mask, wk, Δh, ω, F, tol)

    res0 = resJCB(θ, λ, b, mask, Δh, ω)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:Constant.ItrMax
        res = jacobi!(θ, λ, b, mask, Δh, ω, wk) / res0
        #res = rbsor!(θ, λ, b, mask, Δh, ω) / res0
        
        @printf(F, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end

"""
@brief 緩和Jacobi法の残差
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    緩和係数
@ret                 1セルあたりの残差RMS
"""
function resJCB(p::Array{Float64,3},
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64)
    SZ = size(p)
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

"""
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
"""
function jacobi!(p::Array{Float64,3},
                 λ::Array{Float64,3}, 
                 b::Array{Float64,3},
                 m::Array{Float64,3}, 
                 Δh, 
                 ω::Float64,
                wk::Array{Float64,3})
    SZ = size(p)
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
=#