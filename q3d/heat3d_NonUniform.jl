module NonUniform
export exact_solution!, read_grid_file, solveSOR!

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

# Adiabatic Mask
# @param a value
# @param ma mask for left
# @param mb mask for right
λa(a, ma, mb) = a*div(ma+mb,2)

"""
    read_grid_file(filename)

ASCIIファイルから格子点座標を読み込む

# Arguments
- `filename::String`: 入力ファイル名

# Returns
- `coord::Vector{Float64}`: 格子点座標配列
- `numNodes::Int`: 格子点数
"""
function read_grid_file(filename::String="grid.dat")
  if !isfile(filename)
    error("ファイルが見つかりません: $filename")
  end
  
  coord = Float64[]
  numNodes = 0
  
  open(filename, "r") do f
    # 1行目: 格子点数を読み込み
    line = readline(f)
    numNodes = parse(Int, strip(line))
    
    # 座標配列を初期化
    coord = zeros(Float64, numNodes)
    
    # 格子点データを読み込み
    for i in 1:numNodes
      line = readline(f)
      parts = split(strip(line))
      
      if length(parts) != 2
        error("ファイル形式エラー（行$(i+1)）: $line")
      end
      
      grid_index = parse(Int, parts[1])
      grid_coord = parse(Float64, parts[2])
      
      # 格子点番号の整合性チェック（1オリジン）
      if grid_index != i
        @warn "格子点番号が期待値と異なります: 期待値=$i, 実際値=$grid_index"
      end
      
      coord[i] = grid_coord
    end
  end
  
  println("格子点データ読み込み完了:")
  @printf("  ファイル: %s\n", filename)
  @printf("  格子点数: %d\n", numNodes)
  @printf("  座標範囲: [%.6f, %.6f]\n", minimum(coord), maximum(coord))
  
  return coord, numNodes
end

#=
@brief 厳密解
@param [out]    e      解ベクトル
@param [in]     SZ     配列長
@param [in]     Δh     セル幅
@param [in]     Z      
=#
function exact_solution!(e::Array{Float64,3}, SZ, ox, Δh, Z::Vector{Float64})
    r2 = sqrt(2.0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        z = Z[k]
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
        p[i,j,2    ] = a
        p[i,j,SZ[3]-1] = a
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_st Zループ開始インデクス
@param [in]     z_ed Zループ終了インデクス
@param [in]     F    ファイルディスクリプタ
=#
function solveSOR!(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, z_st, z_ed, F)

    res0 = resSOR(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, z_st, z_ed)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        res = sor!(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, z_st, z_ed) / res0
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
@brief SOR法の残差
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    加速係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_st Zループ開始インデクス
@param [in]     z_ed Zループ終了インデクス
@ret                 1セルあたりの残差RMS
=#
function resSOR(p::Array{Float64,3}, 
                SZ, 
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64,
                Z::Vector{Float64},
                ΔZ::Vector{Float64},
                z_st::Int, 
                z_ed::Int
                )

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λa(λ0, m[i,j,k+1], m0) / (ΔZ[k]*(Z[k+1]-Z[k]))
        ab = λa(λ0, m[i,j,k-1], m0) / (ΔZ[k]*(Z[k]-Z[k-1]))
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = (dd + ω*(aw+as+ab))*dp / ω
        res += r*r
    end
    
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_st Zループ開始インデクス
@param [in]     z_ed Zループ終了インデクス
@ret                 1セルあたりの残差RMS
=#
function sor!(p::Array{Float64,3}, 
              SZ, 
              λ::Array{Float64,3}, 
              b::Array{Float64,3},
              m::Array{Float64,3}, 
              Δh, 
              ω::Float64,
              Z::Vector{Float64},
             ΔZ::Vector{Float64},
             z_st::Int, 
             z_ed::Int)

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λa(λ0, m[i,j,k+1], m0) / (ΔZ[k]*(Z[k+1]-Z[k]))
        ab = λa(λ0, m[i,j,k-1], m0) / (ΔZ[k]*(Z[k]-Z[k-1]))
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

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_st Zループ開始インデクス
@param [in]     z_ed Zループ終了インデクス
@param [in]     F    ファイルディスクリプタ
=#
function solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, Z, ΔZ, z_st, z_ed, F)

    res0 = resJCB(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, z_st, z_ed)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        res = jacobi!(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, z_st, z_ed, wk) / res0
        #res = rbsor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        
        @printf(F, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_st Zループ開始インデクス
@param [in]     z_ed Zループ終了インデクス
@ret                 1セルあたりの残差RMS
=#
function resJCB(p::Array{Float64,3},
                SZ,
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                m::Array{Float64,3}, 
                Δh, 
                ω::Float64, 
                Z::Vector{Float64}, 
                ΔZ::Vector{Float64}, 
                z_st::Int, 
                z_ed::Int)

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λa(λ0, m[i,j,k+1], m0) / (ΔZ[k]*(Z[k+1]-Z[k]))
        ab = λa(λ0, m[i,j,k-1], m0) / (ΔZ[k]*(Z[k]-Z[k-1]))
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


#=
@brief 緩和Jacobi法
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     λ    熱伝導率
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     Δh   セル幅
@param [in] 　　　　　　　　ω    緩和係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_st Zループ開始インデクス
@param [in]     z_ed Zループ終了インデクス
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
                 Z::Vector{Float64}, 
                ΔZ::Vector{Float64}, 
              z_st::Int, 
              z_ed::Int,
                wk::Array{Float64,3})

    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)

    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λa(λ0, m[i,j,k+1], m0) / (ΔZ[k]*(Z[k+1]-Z[k]))
        ab = λa(λ0, m[i,j,k-1], m0) / (ΔZ[k]*(Z[k]-Z[k-1]))
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




end # end of module