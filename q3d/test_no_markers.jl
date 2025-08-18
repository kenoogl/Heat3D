include("convergence_history.jl")

# テストデータ作成
conv_data = ConvergenceData("pbicgstab", "gs")
for i in 1:20
    residual = 1.0 * (0.8)^i
    add_residual!(conv_data, i, residual)
end

# マーカー付きとマーカーなしの比較
plot_convergence_curve(conv_data, "test_with_markers.png", show_markers=true)
plot_convergence_curve(conv_data, "test_no_markers.png", show_markers=false)

println("生成されたファイル:")
println("- test_with_markers.png (マーカー付き)")
println("- test_no_markers.png (マーカーなし)")