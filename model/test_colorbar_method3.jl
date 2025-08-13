using Plots
using Printf

println("方法3: カラーバーの指数表示テスト")

# テストデータ生成（指数範囲のデータ）
x = range(0, 1, length=30)
z = range(0, 1, length=30)
# 1e-7から1e-4の範囲のデータを生成
data = [1e-7 * exp(6 * (xi + zi)) for xi in x, zi in z]

println("データ範囲: $(minimum(data)) ~ $(maximum(data))")

# 方法3: 自動計算で指数表示
min_val = minimum(data)
max_val = maximum(data)
n_ticks = 6

# 対数スケールでティック値を計算
log_min = log10(min_val)
log_max = log10(max_val)
log_ticks = range(log_min, log_max, length=n_ticks)
auto_tick_values = [10^x for x in log_ticks]
auto_tick_labels = [@sprintf("%.1E", v) for v in auto_tick_values]

println("自動生成ティック値: ", auto_tick_values)
println("自動生成ラベル: ", auto_tick_labels)

# contourプロット作成（指数表示カラーバー付き）
p = contour(z, x, data', 
    fill=true, 
    c=:thermal,
    colorbar=true,
    colorbar_ticks=(auto_tick_values, auto_tick_labels),
    colorbar_title="Temperature [K]",
    colorbar_titlefontsize=14,
    colorbar_tickfontsize=12,
    xlabel="Z-coordinate [m]",
    ylabel="X-coordinate [m]",
    title="方法3: 自動計算による指数表示カラーバー", 
    size=(800, 500),
    aspect_ratio=:equal)

savefig(p, "colorbar_method3_test.png")
println("テスト完了: colorbar_method3_test.png が生成されました")

# さらに詳細なテスト：温度拡散率データ
println("\n--- 温度拡散率データでのテスト ---")
# 実際の温度拡散率範囲 (m²/s)
alpha_data = [5.52e-7 * (1 + 200*xi*zi) for xi in x, zi in z]  # Resin ~ Silicon程度

min_alpha = minimum(alpha_data)
max_alpha = maximum(alpha_data)
println("温度拡散率範囲: $(min_alpha) ~ $(max_alpha)")

# 温度拡散率用の指数表示
log_min_alpha = log10(min_alpha)
log_max_alpha = log10(max_alpha)
log_ticks_alpha = range(log_min_alpha, log_max_alpha, length=n_ticks)
alpha_tick_values = [10^x for x in log_ticks_alpha]
alpha_tick_labels = [@sprintf("%.1E", v) for v in alpha_tick_values]

p2 = contour(z, x, alpha_data', 
    fill=true, 
    c=:viridis,
    colorbar=true,
    colorbar_ticks=(alpha_tick_values, alpha_tick_labels),
    colorbar_title="α [m²/s]",
    colorbar_titlefontsize=14,
    colorbar_tickfontsize=12,
    xlabel="Z-coordinate [m]",
    ylabel="X-coordinate [m]",
    title="温度拡散率分布（指数表示カラーバー）", 
    size=(800, 500),
    aspect_ratio=:equal)

savefig(p2, "alpha_distribution_exp.png")
println("温度拡散率テスト完了: alpha_distribution_exp.png が生成されました")

# 追加テスト: heatmapでも指数表示
println("\n--- heatmapでの指数表示テスト ---")
p3 = heatmap(z, x, alpha_data', 
    c=:plasma,
    colorbar=true,
    colorbar_ticks=(alpha_tick_values, alpha_tick_labels),
    colorbar_title="α [m²/s]",
    colorbar_titlefontsize=14,
    colorbar_tickfontsize=12,
    xlabel="Z-coordinate [m]",
    ylabel="X-coordinate [m]",
    title="heatmap: 指数表示カラーバー", 
    size=(800, 500),
    aspect_ratio=:equal)

savefig(p3, "heatmap_exp_colorbar.png")
println("heatmapテスト完了: heatmap_exp_colorbar.png が生成されました")