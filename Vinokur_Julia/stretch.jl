using Printf
using Plots

const EPS = 1.0e-5
const ITER_MAX = 20

"""
    stretch(x1, x2, sp1, sp2, numNodes; debug=false)

Vinokurのストレッチング関数を用いて格子点分布を計算する

# Arguments
- `x1::Float64`: 開始点座標
- `x2::Float64`: 終了点座標
- `sp1::Float64`: 開始点側の格子点間隔
- `sp2::Float64`: 終了点側の格子点間隔
- `numNodes::Int`: 分割格子点数
- `debug::Bool`: デバッグ出力の有無（デフォルト: false）

# Returns
- `coord::Vector{Float64}`: 格子点座標配列
- `success::Bool`: 収束成功フラグ
"""
function stretch(x1::Float64, x2::Float64, sp1::Float64, sp2::Float64, numNodes::Int; debug::Bool=false)
  length = abs(x2 - x1)
  n1 = Float64(numNodes - 1)
  
  # 初期化：[0,1]に正規化
  x = zeros(Float64, numNodes)
  wk = zeros(Float64, numNodes)
  
  for j in 1:numNodes
    x[j] = Float64(j-1) / n1
    wk[j] = Float64(j-1) / n1
  end
  
  # [0,1]で分布を計算するためlengthで正規化
  sp1_norm = sp1 / length
  sp2_norm = sp2 / length
  
  if sp1_norm > 1.0 || sp2_norm > 1.0
    @warn "格子間隔が全長より大きいです"
    return zeros(Float64, numNodes), false
  end
  
  s11 = 1.0 / n1 / sp1_norm
  s21 = 1.0 / n1 / sp2_norm
  
  stretching!(wk, x, s11, s21)
  dx1 = wk[2] - wk[1]
  dx2 = wk[numNodes] - wk[numNodes-1]
  
  debug_data = []
  
  itr = 0
  if abs(sp1_norm - dx1) / sp1_norm < EPS && abs(sp2_norm - dx2) / sp2_norm < EPS
    @goto convergence
  end
  
  for iter in 1:ITER_MAX
    itr = iter
    ds1 = s11 * 0.1
    ds2 = s21 * 0.1
    s12 = s11 + ds1
    s22 = s21 + ds2
    
    # 数値微分による勾配計算
    stretching!(wk, x, s12, s21)
    dx1_ds1 = (wk[2] - wk[1] - dx1) / ds1
    dx2_ds1 = (wk[numNodes] - wk[numNodes-1] - dx2) / ds1
    
    stretching!(wk, x, s11, s22)
    dx1_ds2 = (wk[2] - wk[1] - dx1) / ds2
    dx2_ds2 = (wk[numNodes] - wk[numNodes-1] - dx2) / ds2
    
    df1 = sp1_norm - dx1
    df2 = sp2_norm - dx2
    
    # ニュートン法の更新
    d = dx1_ds1 * dx2_ds2 - dx1_ds2 * dx2_ds1
    if abs(d) < 1e-15
      @warn "ヤコビアンが特異です"
      return zeros(Float64, numNodes), false
    end
    
    ds1 = (df1 * dx2_ds2 - df2 * dx1_ds2) / d
    ds2 = (df2 * dx1_ds1 - df1 * dx2_ds1) / d
    
    s12 = s11 + ds1
    s22 = s21 + ds2
    
    stretching!(wk, x, s12, s22)
    dx1 = wk[2] - wk[1]
    dx2 = wk[numNodes] - wk[numNodes-1]
    
    if debug
      push!(debug_data, copy(wk))
      @printf("Iter %2d: dx1=%8.5f dx2=%8.5f error1=%8.5e error2=%8.5e\n", 
              iter, dx1, dx2, abs(sp1_norm-dx1)/sp1_norm, abs(sp2_norm-dx2)/sp2_norm)
    end
    
    if abs(sp1_norm - dx1) / sp1_norm < EPS && abs(sp2_norm - dx2) / sp2_norm < EPS
      break
    end
    
    s11 = s12
    s21 = s22
  end
  
  @label convergence
  # 境界値を正確に設定
  wk[1] = 0.0
  wk[numNodes] = 1.0
  
  # 実座標系にスケールバック
  coord = zeros(Float64, numNodes)
  for j in 1:numNodes
    coord[j] = x1 + wk[j] * length
  end
  
  if debug
    println("\n=== 収束結果 ===")
    @printf("反復回数: %d\n", itr)
    println("格子点座標:")
    for i in 1:numNodes
      if i == 1
        @printf("%3d %8.5f %8.5f\n", i-1, wk[i], coord[i])
      else
        @printf("%3d %8.5f %8.5f %8.5f\n", i-1, wk[i], coord[i], coord[i]-coord[i-1])
      end
    end
  end
  
  success = (itr < ITER_MAX)
  return coord, success
end

"""
両端がs0, s1の間隔の分布を近似する内部関数
"""
function stretching!(wk::Vector{Float64}, x::Vector{Float64}, s0::Float64, s1::Float64)
  epsilon = 0.001
  numNodes = length(x)
  
  b = sqrt(s0 * s1)
  a = b / s1
  
  if b < 1.0 - epsilon
    # asin分枝
    dz = asin(b)
    for j in 1:numNodes
      tanx = tan(dz * x[j])
      wk[j] = tanx / (a * sin(dz) + (1.0 - a * cos(dz)) * tanx)
    end
  elseif b > 1.0 + epsilon
    # asinh分枝
    dz = asinh(b)
    for j in 1:numNodes
      tanhx = tanh(dz * x[j])
      wk[j] = tanhx / (a * sinh(dz) + (1.0 - a * cosh(dz)) * tanhx)
    end
  else
    # 線形近似分枝
    for j in 1:numNodes
      u = x[j] * (1.0 + 2.0 * (b - 1.0) * (x[j] - 0.5) * (1.0 - x[j]))
      wk[j] = u / (a + (1.0 - a) * u)
    end
  end
end

"""
格子分布を可視化する関数
"""
function plot_grid_distribution(coord::Vector{Float64}, title_str::String="Grid Distribution")
  numNodes = length(coord)
  
  # 格子点間隔を計算
  spacings = zeros(Float64, numNodes-1)
  for i in 1:numNodes-1
    spacings[i] = coord[i+1] - coord[i]
  end
  
  # プロット作成
  p = plot(layout=(2,1), size=(800, 600))
  
  # 上段：格子点分布
  scatter!(p[1], coord, ones(numNodes), 
           xlabel="Position", ylabel="", 
           title=title_str, 
           legend=false, 
           markersize=3, 
           ylims=(0.5, 1.5))
  
  # 下段：格子点間隔
  plot!(p[2], 1:numNodes-1, spacings, 
        xlabel="Grid Index", ylabel="Grid Spacing", 
        title="Grid Spacing Distribution", 
        legend=false, 
        linewidth=2, 
        marker=:circle, 
        markersize=3)
  
  return p
end