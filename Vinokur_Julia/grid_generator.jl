include("stretch.jl")
using Printf
using Plots

"""
格子点分布生成・可視化・ファイル入出力機能
"""

"""
    generate_and_save_grid(x1, x2, sp1, sp2, numNodes, filename; debug=false)

Vinokurストレッチング関数を使って格子点分布を生成し、
可視化とファイル出力を行う

# Arguments
- `x1::Float64`: 開始点座標
- `x2::Float64`: 終了点座標
- `sp1::Float64`: 開始点側の格子点間隔
- `sp2::Float64`: 終了点側の格子点間隔
- `numNodes::Int`: 格子点数
- `filename::String`: 出力ファイル名（デフォルト: "grid.dat"）
- `debug::Bool`: デバッグ出力の有無

# Returns
- `coord::Vector{Float64}`: 格子点座標配列
- `success::Bool`: 計算成功フラグ
"""
function generate_and_save_grid(x1::Float64, x2::Float64, sp1::Float64, sp2::Float64, 
                               numNodes::Int, filename::String="grid.txt"; debug::Bool=false)
  
  # 格子点分布計算
  coord, success = stretch(x1, x2, sp1, sp2, numNodes, debug=debug)
  
  if !success
    @warn "格子点分布の計算に失敗しました"
    return coord, success
  end
  
  # ファイル出力
  write_grid_file(coord, filename)
  
  # 可視化
  plot_title = @sprintf("Grid Distribution: %d nodes, sp1=%.4f, sp2=%.4f", 
                       numNodes, sp1, sp2)
  p = plot_grid_distribution(coord, plot_title)
  
  # グラフファイル名を生成
  base_name = splitext(filename)[1]
  plot_filename = "$(base_name)_plot.png"
  savefig(p, plot_filename)
  
  println("格子点分布生成完了:")
  @printf("  格子点数: %d\n", numNodes)
  @printf("  座標範囲: [%.4f, %.4f]\n", x1, x2)
  @printf("  格子間隔: %.4f (開始) → %.4f (終了)\n", sp1, sp2)
  @printf("  データファイル: %s\n", filename)
  @printf("  グラフファイル: %s\n", plot_filename)
  
  return coord, success
end

"""
    write_grid_file(coord, filename)

格子点座標をASCIIファイルに出力する

ファイル形式:
- 1行目: 格子点数
- 2行目以降: 格子点番号（1から開始、Julia式）と座標値のペア

# Arguments
- `coord::Vector{Float64}`: 格子点座標配列
- `filename::String`: 出力ファイル名
"""
function write_grid_file(coord::Vector{Float64}, filename::String)
  open(filename, "w") do f
    # 1行目: 格子点数
    println(f, length(coord))
    
    # 2行目以降: 格子点番号と座標値
    for (i, x) in enumerate(coord)
      @printf(f, "%d %.16e\n", i, x)  # 格子点番号は1から開始（Julia式）
    end
  end
end

"""
    read_grid_file(filename)

ASCIIファイルから格子点座標を読み込む

# Arguments
- `filename::String`: 入力ファイル名

# Returns
- `coord::Vector{Float64}`: 格子点座標配列
- `numNodes::Int`: 格子点数
"""
function read_grid_file(filename::String)
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

"""
    verify_grid_file(original_coord, filename)

元の座標配列とファイルから読み込んだ座標配列を比較検証する

# Arguments
- `original_coord::Vector{Float64}`: 元の格子点座標配列
- `filename::String`: 検証するファイル名

# Returns
- `is_identical::Bool`: 一致判定フラグ
- `max_error::Float64`: 最大誤差
"""
function verify_grid_file(original_coord::Vector{Float64}, filename::String)
  # ファイルから読み込み
  read_coord, numNodes = read_grid_file(filename)
  
  # 配列長チェック
  if length(original_coord) != length(read_coord)
    @warn "配列長が異なります: 元=$(length(original_coord)), 読込=$(length(read_coord))"
    return false, Inf
  end
  
  # 座標値の比較
  errors = abs.(original_coord .- read_coord)
  max_error = maximum(errors)
  is_identical = max_error < 1e-15  # 数値精度を考慮した許容誤差
  
  println("ファイル検証結果:")
  @printf("  最大誤差: %.2e\n", max_error)
  @printf("  一致判定: %s\n", is_identical ? "✓ 一致" : "✗ 不一致")
  
  if !is_identical
    # 誤差の詳細を表示
    worst_index = argmax(errors)
    @printf("  最大誤差位置: インデックス%d\n", worst_index)
    @printf("  元の値: %.16e\n", original_coord[worst_index])
    @printf("  読込値: %.16e\n", read_coord[worst_index])
  end
  
  return is_identical, max_error
end

"""
    compare_grid_spacings(coord1, coord2, label1="Grid 1", label2="Grid 2")

2つの格子点分布の格子間隔を比較する

# Arguments
- `coord1, coord2::Vector{Float64}`: 格子点座標配列
- `label1, label2::String`: ラベル
"""
function compare_grid_spacings(coord1::Vector{Float64}, coord2::Vector{Float64}, 
                              label1::String="Grid 1", label2::String="Grid 2")
  
  if length(coord1) != length(coord2)
    @warn "格子点数が異なります: $(length(coord1)) vs $(length(coord2))"
    return
  end
  
  # 格子間隔を計算
  spacing1 = diff(coord1)
  spacing2 = diff(coord2)
  
  # 比較プロット
  p = plot(layout=(2,1), size=(800, 800))
  
  # 上段: 格子点分布
  plot!(p[1], 1:length(coord1), coord1, label=label1, linewidth=2, marker=:circle, markersize=3)
  plot!(p[1], 1:length(coord2), coord2, label=label2, linewidth=2, marker=:square, markersize=3)
  xlabel!(p[1], "Grid Index")
  ylabel!(p[1], "Coordinate")
  title!(p[1], "Grid Point Distribution Comparison")
  
  # 下段: 格子間隔比較
  plot!(p[2], 1:length(spacing1), spacing1, label="$label1 spacing", linewidth=2, marker=:circle, markersize=3)
  plot!(p[2], 1:length(spacing2), spacing2, label="$label2 spacing", linewidth=2, marker=:square, markersize=3)
  xlabel!(p[2], "Grid Index")
  ylabel!(p[2], "Grid Spacing")
  title!(p[2], "Grid Spacing Comparison")
  
  return p
end

"""
デモンストレーション関数
"""
function demo_grid_generation()
  println("=== Vinokur格子生成デモンストレーション ===\n")
  
  # テストケース1: 左側密集
  println("テストケース1: 左側密集格子")
  coord1, success1 = generate_and_save_grid(0.0, 1.0, 0.01, 0.1, 51, "grid_left_clustering.dat", debug=true)
  
  if success1
    # ファイル検証
    verify_grid_file(coord1, "grid_left_clustering.dat")
  end
  println()
  
  # テストケース2: 右側密集
  println("テストケース2: 右側密集格子")
  coord2, success2 = generate_and_save_grid(0.0, 1.0, 0.1, 0.01, 51, "grid_right_clustering.dat", debug=false)
  println()
  
  # テストケース3: 均等分布
  println("テストケース3: 均等分布格子")
  coord3, success3 = generate_and_save_grid(0.0, 1.0, 0.02, 0.02, 51, "grid_uniform.dat", debug=false)
  println()
  
  # 比較可視化
  if success1 && success2
    println("格子分布比較グラフを作成中...")
    p_compare = compare_grid_spacings(coord1, coord2, "Left Clustering", "Right Clustering")
    savefig(p_compare, "grid_comparison.png")
    println("比較グラフ保存: grid_comparison.png")
  end
  
  # ファイル読み込みテスト
  println("\n=== ファイル読み込みテスト ===")
  if success1
    println("左側密集格子の読み込みテスト:")
    read_coord, numNodes = read_grid_file("grid_left_clustering.dat")
    @printf("読み込み完了: %d点の格子\n", numNodes)
  end
  
  println("\nデモンストレーション完了!")
end

# デモ実行（スクリプトとして実行された場合）
if abspath(PROGRAM_FILE) == @__FILE__
  demo_grid_generation()
end

generate_and_save_grid(0.0, 1.0, 0.02, 0.06, 25)