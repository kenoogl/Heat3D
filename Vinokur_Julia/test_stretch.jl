include("stretch.jl")

"""
Vinokurストレッチング関数のテストスイート
"""

function test_stretch_function()
  println("=== Vinokur Stretching Function Test ===\n")
  
  # テストケース1: 均等分布
  println("テストケース1: 均等分布")
  println("x1=0.0, x2=1.0, sp1=0.1, sp2=0.1, numNodes=11")
  coord1, success1 = stretch(0.0, 1.0, 0.1, 0.1, 11, debug=true)
  
  if success1
    println("✓ 収束成功")
    # 理論値との比較（均等分布の場合）
    expected_spacing = 0.1
    actual_spacing_start = coord1[2] - coord1[1]
    actual_spacing_end = coord1[11] - coord1[10]
    @printf("期待値: %8.5f, 実際値(開始): %8.5f, 実際値(終了): %8.5f\n", 
            expected_spacing, actual_spacing_start, actual_spacing_end)
    
    # 可視化
    p1 = plot_grid_distribution(coord1, "Test Case 1: Uniform Distribution")
    savefig(p1, "test_case1_uniform.png")
    println("プロット保存: test_case1_uniform.png\n")
  else
    println("✗ 収束失敗\n")
  end
  
  # テストケース2: 左側密集
  println("テストケース2: 左側密集")
  println("x1=0.0, x2=1.0, sp1=0.02, sp2=0.2, numNodes=11")
  coord2, success2 = stretch(0.0, 1.0, 0.02, 0.2, 11, debug=true)
  
  if success2
    println("✓ 収束成功")
    actual_spacing_start = coord2[2] - coord2[1]
    actual_spacing_end = coord2[11] - coord2[10]
    @printf("期待値(開始): %8.5f, 実際値(開始): %8.5f\n", 0.02, actual_spacing_start)
    @printf("期待値(終了): %8.5f, 実際値(終了): %8.5f\n", 0.2, actual_spacing_end)
    
    # 可視化
    p2 = plot_grid_distribution(coord2, "Test Case 2: Left Clustering")
    savefig(p2, "test_case2_left_clustering.png")
    println("プロット保存: test_case2_left_clustering.png\n")
  else
    println("✗ 収束失敗\n")
  end
  
  # テストケース3: 右側密集
  println("テストケース3: 右側密集")
  println("x1=0.0, x2=1.0, sp1=0.2, sp2=0.02, numNodes=11")
  coord3, success3 = stretch(0.0, 1.0, 0.2, 0.02, 11, debug=true)
  
  if success3
    println("✓ 収束成功")
    actual_spacing_start = coord3[2] - coord3[1]
    actual_spacing_end = coord3[11] - coord3[10]
    @printf("期待値(開始): %8.5f, 実際値(開始): %8.5f\n", 0.2, actual_spacing_start)
    @printf("期待値(終了): %8.5f, 実際値(終了): %8.5f\n", 0.02, actual_spacing_end)
    
    # 可視化
    p3 = plot_grid_distribution(coord3, "Test Case 3: Right Clustering")
    savefig(p3, "test_case3_right_clustering.png")
    println("プロット保存: test_case3_right_clustering.png\n")
  else
    println("✗ 収束失敗\n")
  end
  
  # テストケース4: 異なる座標範囲
  println("テストケース4: 異なる座標範囲")
  println("x1=2.0, x2=5.0, sp1=0.1, sp2=0.5, numNodes=11")
  coord4, success4 = stretch(2.0, 5.0, 0.1, 0.5, 11, debug=true)
  
  if success4
    println("✓ 収束成功")
    actual_spacing_start = coord4[2] - coord4[1]
    actual_spacing_end = coord4[11] - coord4[10]
    @printf("期待値(開始): %8.5f, 実際値(開始): %8.5f\n", 0.1, actual_spacing_start)
    @printf("期待値(終了): %8.5f, 実際値(終了): %8.5f\n", 0.5, actual_spacing_end)
    
    # 可視化
    p4 = plot_grid_distribution(coord4, "Test Case 4: Different Coordinate Range")
    savefig(p4, "test_case4_different_range.png")
    println("プロット保存: test_case4_different_range.png\n")
  else
    println("✗ 収束失敗\n")
  end
  
  # テストケース5: 極端なケース（収束困難）
  println("テストケース5: 極端なケース")
  println("x1=0.0, x2=1.0, sp1=0.001, sp2=0.5, numNodes=21")
  coord5, success5 = stretch(0.0, 1.0, 0.001, 0.5, 21, debug=true)
  
  if success5
    println("✓ 収束成功")
    actual_spacing_start = coord5[2] - coord5[1]
    actual_spacing_end = coord5[21] - coord5[20]
    @printf("期待値(開始): %8.5f, 実際値(開始): %8.5f\n", 0.001, actual_spacing_start)
    @printf("期待値(終了): %8.5f, 実際値(終了): %8.5f\n", 0.5, actual_spacing_end)
    
    # 可視化
    p5 = plot_grid_distribution(coord5, "Test Case 5: Extreme Case")
    savefig(p5, "test_case5_extreme.png")
    println("プロット保存: test_case5_extreme.png\n")
  else
    println("✗ 収束失敗\n")
  end
  
  # まとめ
  println("=== テスト結果まとめ ===")
  results = [success1, success2, success3, success4, success5]
  test_names = ["均等分布", "左側密集", "右側密集", "異なる座標範囲", "極端なケース"]
  
  for (i, (result, name)) in enumerate(zip(results, test_names))
    status = result ? "✓ 成功" : "✗ 失敗"
    println("テストケース$i ($name): $status")
  end
  
  success_rate = sum(results) / length(results) * 100
  @printf("成功率: %.1f%% (%d/%d)\n", success_rate, sum(results), length(results))
end

# パフォーマンステスト
function performance_test()
  println("\n=== パフォーマンステスト ===")
  
  # 大きなノード数でのテスト
  println("大きなノード数でのテスト (numNodes=101)")
  @time coord, success = stretch(0.0, 1.0, 0.01, 0.1, 101, debug=false)
  
  if success
    println("✓ 収束成功")
    @printf("開始格子間隔: %8.5f, 終了格子間隔: %8.5f\n", 
            coord[2]-coord[1], coord[101]-coord[100])
  else
    println("✗ 収束失敗")
  end
end

# エラーハンドリングテスト
function error_handling_test()
  println("\n=== エラーハンドリングテスト ===")
  
  # 不正な格子間隔（全長より大きい）
  println("不正な格子間隔テスト")
  coord, success = stretch(0.0, 1.0, 2.0, 0.1, 11, debug=false)
  if !success
    println("✓ 適切にエラーハンドリング")
  else
    println("✗ エラーハンドリング失敗")
  end
end

# メイン実行
function main()
  test_stretch_function()
  performance_test()
  error_handling_test()
  
  println("\nすべてのテストが完了しました。")
end

# テスト実行
main()