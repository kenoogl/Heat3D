include("grid_generator.jl")
using Statistics

"""
格子生成・ファイル入出力の詳細テスト
"""

function test_grid_io()
  println("=== 格子生成・ファイル入出力テスト ===\n")
  
  # テストケース1: 軽微な非均等分布（収束しやすいパラメータ）
  println("テストケース1: 軽微な左側密集格子")
  coord1, success1 = generate_and_save_grid(0.0, 1.0, 0.015, 0.025, 21, "test_mild_left.dat")
  
  if success1
    println("✓ 格子生成成功")
    
    # ファイル検証
    is_ok, max_err = verify_grid_file(coord1, "test_mild_left.dat")
    println("✓ ファイル入出力検証:", is_ok ? "成功" : "失敗")
    
    # 格子品質チェック
    spacings = diff(coord1)
    @printf("  最小格子間隔: %.6f\n", minimum(spacings))
    @printf("  最大格子間隔: %.6f\n", maximum(spacings))
    @printf("  格子間隔比: %.2f\n", maximum(spacings)/minimum(spacings))
  else
    println("✗ 格子生成失敗")
  end
  println()
  
  # テストケース2: 均等分布（確実に成功）
  println("テストケース2: 均等分布格子")
  coord2, success2 = generate_and_save_grid(0.0, 2.0, 0.1, 0.1, 21, "test_uniform.dat")
  
  if success2
    println("✓ 格子生成成功")
    
    # ファイル検証
    is_ok, max_err = verify_grid_file(coord2, "test_uniform.dat")
    println("✓ ファイル入出力検証:", is_ok ? "成功" : "失敗")
    
    # 均等性チェック
    spacings = diff(coord2)
    spacing_std = std(spacings)
    @printf("  格子間隔標準偏差: %.2e (均等度指標)\n", spacing_std)
  else
    println("✗ 格子生成失敗")
  end
  println()
  
  # テストケース3: 異なる座標範囲
  println("テストケース3: 異なる座標範囲")
  coord3, success3 = generate_and_save_grid(-1.0, 3.0, 0.2, 0.2, 21, "test_range.dat")
  
  if success3
    println("✓ 格子生成成功")
    @printf("  座標範囲: [%.2f, %.2f]\n", minimum(coord3), maximum(coord3))
    
    # ファイル検証
    is_ok, max_err = verify_grid_file(coord3, "test_range.dat")
    println("✓ ファイル入出力検証:", is_ok ? "成功" : "失敗")
  else
    println("✗ 格子生成失敗")
  end
  println()
  
  # 比較可視化
  if success2 && success3
    println("格子分布比較...")
    p = compare_grid_spacings(coord2, coord3, "Uniform [0,2]", "Uniform [-1,3]")
    savefig(p, "test_comparison.png")
    println("✓ 比較グラフ保存: test_comparison.png")
  end
  
  # ファイル内容の直接確認
  println("\n=== 生成ファイルの内容確認 ===")
  if success2
    println("test_uniform.dat の先頭5行:")
    run(pipeline(`head -6 test_uniform.dat`, stdout))
  end
end

# ファイル読み込み専用テスト
function test_file_reading()
  println("\n=== ファイル読み込み専用テスト ===")
  
  # テスト用ファイルを手動作成（1オリジン）
  test_data = [
    (1, 0.0),
    (2, 0.5),
    (3, 1.0),
    (4, 1.5),
    (5, 2.0)
  ]
  
  open("manual_test.dat", "w") do f
    println(f, length(test_data))
    for (idx, coord) in test_data
      @printf(f, "%d %.16e\n", idx, coord)
    end
  end
  
  println("手動作成ファイル:")
  run(pipeline(`cat manual_test.dat`, stdout))
  
  # 読み込みテスト
  try
    coord, numNodes = read_grid_file("manual_test.dat")
    println("✓ 読み込み成功")
    @printf("  格子点数: %d\n", numNodes)
    @printf("  座標: %s\n", coord)
    
    # 期待値との比較
    expected_coords = [data[2] for data in test_data]
    if coord ≈ expected_coords
      println("✓ データ整合性確認")
    else
      println("✗ データ不整合")
      println("  期待値: $expected_coords")
      println("  実際値: $coord")
    end
  catch e
    println("✗ 読み込み失敗: $e")
  end
end

# 大規模テスト
function test_large_grid()
  println("\n=== 大規模格子テスト ===")
  
  # 201点の大規模格子
  coord, success = generate_and_save_grid(0.0, 10.0, 0.02, 0.02, 201, "large_grid.dat")
  
  if success
    println("✓ 大規模格子生成成功（201点）")
    
    # ファイルサイズ確認
    filesize = stat("large_grid.dat").size
    @printf("  ファイルサイズ: %.2f KB\n", filesize/1024)
    
    # 読み込み性能テスト
    @time read_coord, numNodes = read_grid_file("large_grid.dat")
    println("✓ 大規模格子読み込み完了")
    
    # 検証
    is_ok, max_err = verify_grid_file(coord, "large_grid.dat")
    println("✓ 大規模格子検証:", is_ok ? "成功" : "失敗")
  else
    println("✗ 大規模格子生成失敗")
  end
end

# メイン実行
function main()
  test_grid_io()
  test_file_reading()
  test_large_grid()
  
  println("\n=== 全テスト完了 ===")
  println("生成されたファイル:")
  for file in ["test_mild_left.dat", "test_uniform.dat", "test_range.dat", 
               "manual_test.dat", "large_grid.dat"]
    if isfile(file)
      println("  ✓ $file")
    else
      println("  ✗ $file (未生成)")
    end
  end
end

main()