module NonUniform
export exact_solution!, solveSOR!, PBiCGSTAB!

using Printf

include("common.jl")


"""
@brief 厳密解
@param [out]    e      解ベクトル
@param [in]     ox     領域基点
@param [in]     Δh     セル幅
@param [in]     Z      Z座標
"""
function exact_solution!(e::Array{Float64,3}, ox, Δh, Z::Vector{Float64})
    SZ = size(e)
    r2 = sqrt(2.0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        z = Z[k]
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     F    ファイルディスクリプタ
@param [in]     tol  収束判定基準
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
"""
function solveSOR!(θ, λ, b, mask, Δh, ω, Z, ΔZ, 
                    z_range::Vector{Int64}, 
                    F, tol,
                    HF::Vector{Float64}, 
                    HT::Vector{Float64})

    res0 = resSOR(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range, HF, HT)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:Constant.ItrMax
        #res = sor!(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range, HF, HT) / res0
        res = rbsor!(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range, HF, HT) / res0
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
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
                Z::Vector{Float64},
                ΔZ::Vector{Float64},
                z_range::Vector{Int64},
                HF::Vector{Float64}, 
                HT::Vector{Float64})
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +((1.0-mb)*HF[5] - (1.0-m0)*HF[6])/ΔZ[k] ) # 熱流束境界
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]
        
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = (dd + ω*(axm+aym+azm))*dp / ω
        res += r*r
    end
    
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end

"""
@brief SOR法
@param [in,out] p    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
@ret                 セルあたりの残差RMS
"""
function sor!(p::Array{Float64,3}, 
              λ::Array{Float64,3}, 
              b::Array{Float64,3},
              m::Array{Float64,3}, 
              Δh, 
              ω::Float64,
              Z::Vector{Float64},
             ΔZ::Vector{Float64},
             z_range::Vector{Int64},
             HF::Vector{Float64}, 
             HT::Vector{Float64})
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +((1.0-mb)*HF[5] - (1.0-m0)*HF[6])/ΔZ[k] ) # 熱流束境界
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        p[i,j,k] = pp + ω * dp
        r = (dd + ω*(axm+aym+azm))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end


"""
@brief RB-SOR法のカーネル
@param [in,out] p    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
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
                     Z::Vector{Float64},
                    ΔZ::Vector{Float64},
                    z_range::Vector{Int64},
                    HF::Vector{Float64},
                    HT::Vector{Float64},
                    color::Int)
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            pp = p[i,j,k]
            λ0 = λ[i,j,k]
            m0 = m[i,j,k]
            me = 1.0-m[i+1,j  ,k  ]
            mw = 1.0-m[i-1,j  ,k  ]
            mn = 1.0-m[i  ,j+1,k  ]
            ms = 1.0-m[i  ,j-1,k  ]
            mt = m[i  ,j  ,k+1]
            mb = m[i  ,j  ,k-1]
            bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                           +(ms*HF[3] - mn*HF[4])*dy1
                           +((1.0-mb)*HF[5] - (1.0-m0)*HF[6])/ΔZ[k] ) # 熱流束境界
            axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
            axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
            aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
            ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
            zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
            zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
            azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]
            azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]
            dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
            ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
                 + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
                 + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
            dp = (((ss-bb)/dd - pp)) * m0
            p[i,j,k] = pp + ω * dp
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
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
                Z::Vector{Float64},
               ΔZ::Vector{Float64},
               z_range::Vector{Int64},
               HF::Vector{Float64}, 
               HT::Vector{Float64})
    SZ = size(b)
    res::Float64 = 0.0

    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        res += rbsor_core!(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range, HF, HT, c)
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
@param [in]     Z      Z座標
@param [in]     ΔZ     格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     ItrMax 前処理反復数
@param [in]     smoother ["gs", ""]
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
                    Z::Vector{Float64}, 
                   ΔZ::Vector{Float64}, 
                 z_range::Vector{Int64},
             smoother::String,
                    F,
                    mode,
                    tol,
                    HF::Vector{Float64},
                    HT::Vector{Float64})
    SZ = size(X)
    pcg_q .= 0.0  #fill!(pcg_q, 0.0)
    res0 = CalcRK!(pcg_r, X, B, λ, mask, Δh, Z, ΔZ, z_range, HF, HT)
    println("Inital residual = ", res0)
    pcg_r0 .= pcg_r  #copy!(pcg_r0, pcg_r)

    rho_old::Float64 = 1.0
    alpha::Float64 = 0.0
    omega::Float64  = 1.0
    r_omega::Float64 = -omega
    beta::Float64 = 0.0

    for itr in 1:Constant.ItrMax
        rho = Fdot2(pcg_r, pcg_r0, z_range) # 非計算部分はゼロのこと

        if abs(rho) < Constant.FloatMin
            itr = 0
            break
        end

        if itr == 1
            pcg_p .= pcg_r  #copy!(pcg_p, pcg_r)
        else
            beta = rho / rho_old * alpha / omega
            BiCG1!(pcg_p, pcg_r, pcg_q, beta, omega, z_range)
        end

        pcg_p_ .= 0.0  #fill!(pcg_p_, 0.0)
        Preconditioner!(pcg_p_, pcg_p, λ, mask, Δh, smoother, Z, ΔZ, z_range, HF, HT)

        CalcAX!(pcg_q, pcg_p_, Δh, λ, mask, Z, ΔZ, z_range, HF, HT)
        alpha = rho / Fdot2(pcg_q, pcg_r0, z_range)
        r_alpha = -alpha
        Triad!(pcg_s, pcg_q, pcg_r, r_alpha, z_range)

        pcg_s_ .= 0.0  #fill!(pcg_s_, 0.0)
        Preconditioner!(pcg_s_, pcg_s, λ, mask, Δh, smoother, Z, ΔZ, z_range, HF, HT);

        CalcAX!(pcg_t_, pcg_s_, Δh, λ, mask, Z, ΔZ, z_range, HF, HT)
        omega = Fdot2(pcg_t_, pcg_s, z_range) / Fdot1(pcg_t_, z_range)
        r_omega = -omega

        BICG2!(X, pcg_p_, pcg_s_, alpha , omega, z_range)

        Triad!(pcg_r, pcg_t_, pcg_s, r_omega, z_range)
        res = sqrt(Fdot1(pcg_r, z_range))/((SZ[1]-2)*(SZ[2]-2)*(z_range[2]-z_range[1]+1))
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
@param [in]     b    右辺ベクトル
@param [in]     λ    係数
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
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
                Z::Vector{Float64}, 
               ΔZ::Vector{Float64}, 
             z_range::Vector{Int64},
             HF::Vector{Float64},
             HT::Vector{Float64})
    SZ = size(p)
    res::Float64 = 0.0
    dx0::Float64 = Δh[1]
    dy0::Float64 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        bb = b[i,j,k]+( (mw*HF[1] - me*HF[2])*dx1
                       +(ms*HF[3] - mn*HF[4])*dy1
                       +((1.0-mb)*HF[5] - (1.0-m0)*HF[6])/ΔZ[k] ) # 熱流束境界
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )
        rs = (b[i,j,k] - (ss - dd * p[i,j,k]))* m0
        r[i,j,k] = rs
        res += rs*rs
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end


"""
@brief ベクトルの内積
@param [in]     x    ベクトル
@param [in]     z_range Zループ開始/終了インデクス
@ret            内積
"""
function Fdot1(x::Array{Float64,3}, z_range::Vector{Int64})
    SZ = size(x)
    z_st = z_range[1]
    z_ed = z_range[2]
    y::Float64 = 0.0
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        y += x[i,j,k] * x[i,j,k]
    end
    return y
end



"""
@brief 2ベクトルの内積
@param [in]     x    ベクトル
@param [in]     y    ベクトル
@param [in]     z_range Zループ開始/終了インデクス
@ret            内積
"""
function Fdot2(x::Array{Float64,3}, y::Array{Float64,3}, z_range::Vector{Int64})
    SZ = size(x)
    z_st = z_range[1]
    z_ed = z_range[2]
    xy::Float64 = 0.0
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
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
@param [in]     z_range Zループ開始/終了インデクス
"""
function BiCG1!(p::Array{Float64,3}, 
                r::Array{Float64,3}, 
                q::Array{Float64,3}, 
             beta::Float64, 
              omg::Float64, 
             z_range::Vector{Int64})
    SZ = size(p)
    z_st = z_range[1]
    z_ed = z_range[2]
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        p[i,j,k] = r[i,j,k] + beta * (p[i,j,k] - omg * q[i,j,k])
    end
end


"""
@brief 前処理
@param [in,out] xx   解ベクトル
@param [in]     bb   RHSベクトル
@param [in]     λ    熱伝導率
@param [in]     mask マスク配列
@param [in]     wk   作業ベクトル
@param [in]     Δh   セル幅
@param [in]     smoother  ["jacobi", "dor", ""]
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HF   熱流束境界の値
@param [in]     HT   熱伝達境界の値
"""
function Preconditioner!(xx::Array{Float64,3},
                         bb::Array{Float64,3}, 
                          λ::Array{Float64,3}, 
                       mask::Array{Float64,3},  
                         Δh,
                   smoother::String,
                   Z, ΔZ, 
                   z_range::Vector{Int64},
                   HF::Vector{Float64},
                   HT::Vector{Float64})

    res::Float64 = 0.0
    LCmax::Int = 5

    #=
    if smoother=="jacobi"
        #P_Jacobi!(xx, bb, LCmax, SZ, λ, mask, Δh, wk)
        for _ in 1:LCmax
            res = jacobi!(xx, λ, bb, mask, Δh, 0.8, wk, Z, ΔZ, z_range)
        end
    else
    =#
    if smoother=="gs"
        #P_SOR!(xx, bb, LCmax, SZ, λ, mask, Δh)
        for _ in 1:LCmax
            res = sor!(xx, λ, bb, mask, Δh, 1.0, Z, ΔZ, z_range, HF, HT)
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
@param [in]  z_range Zループ開始/終了インデクス
@param [in]  HF   熱流束境界の値
@param [in]  HT   熱伝達境界の値
"""
function CalcAX!(ap::Array{Float64,3}, 
                  p::Array{Float64,3}, 
                  Δh,
                  λ::Array{Float64,3}, 
                  m::Array{Float64,3},
                  Z::Vector{Float64}, 
                 ΔZ::Vector{Float64}, 
               z_range::Vector{Int64},
               HF::Vector{Float64}, 
               HT::Vector{Float64})
    SZ = size(p)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]
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
@param [in]     z_range Zループ開始/終了インデクス
"""
function Triad!(z::Array{Float64,3}, 
                x::Array{Float64,3}, 
                y::Array{Float64,3}, 
                a::Float64, 
             z_range::Vector{Int64})
    SZ = size(z)
    z_st = z_range[1]
    z_ed = z_range[2]
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] = a * x[i,j,k] + y[i,j,k]
    end
end


"""
@brief BiCGstab 2
@param [in,out] z    ベクトル
@param [in]     y    ベクトル
@param [in]     x    ベクトル
@param [in]     a    係数
@param [in]     b    係数
@param [in]     z_range Zループ開始/終了インデクス
"""
function BICG2!(z::Array{Float64,3}, 
                x::Array{Float64,3}, 
                y::Array{Float64,3}, 
                a::Float64, 
                b::Float64, 
             z_range::Vector{Int64})
    SZ = size(z)
    z_st = z_range[1]
    z_ed = z_range[2]
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] += a * x[i,j,k] + b * y[i,j,k]
    end
end


end # end of module

#=

"""
@brief 緩和ヤコビ法による求解
@param [in/out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    RHSベクトル
@param [in]     mask マスク配列
@param [in]     wk   ワーク配列
@param [in]     Δh   セル幅
@param [in]     ω    緩和係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     F    ファイルディスクリプタ
"""
function solveJACOBI!(θ, λ, b, mask, wk, Δh, ω, Z, ΔZ, z_range::Array{Int64,2}, F, tol)

    res0 = resJCB(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:Constant.ItrMax
        res = jacobi!(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range, wk) / res0
        #res = rbsor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        
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
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    緩和係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@ret                 1セルあたりの残差RMS
"""
function resJCB(p::Array{Float64,3},
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64, 
                Z::Vector{Float64}, 
                ΔZ::Vector{Float64}, 
                z_range::Array{Int64,2})
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mb = m[i,j,k-1]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k]
        at = λ0         / (ΔZ[k]*zt)
        ab = λ[i,j,k-1] / (ΔZ[k]*zb)
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = dd*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end


"""
@brief 緩和Jacobi法
@param [in,out] p    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    緩和係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [out]    wk   ワーク用配列
@ret                 1セルあたりの残差RMS
"""
function jacobi!(p::Array{Float64,3},
                 λ::Array{Float64,3}, 
                 b::Array{Float64,3},
                 m::Array{Float64,3}, 
                 Δh, 
                 ω::Float64,
                 Z::Vector{Float64}, 
                ΔZ::Vector{Float64}, 
              z_range::Array{Int64,2},
                wk::Array{Float64,3})
    SZ = size(p)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    z_st = z_range[1]
    z_ed = z_range[2]

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mb = m[i,j,k-1]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k]
        at = λ0         / (ΔZ[k]*zt)
        ab = λ[i,j,k-1] / (ΔZ[k]*zb)
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

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        p[i,j,k] = wk[i,j,k]
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end
=#