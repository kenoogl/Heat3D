using Printf
include("convergence_history.jl")

# 疑似的な収束履歴データで動作テスト
function test_convergence_history()
    println("=== 収束履歴機能のテスト ===")
    
    # テストデータ1: PBiCGSTAB
    conv_pbicg = ConvergenceData("pbicgstab", "gs")
    
    # 実際のPBiCGSTABに近い残差データを生成
    initial_res = 7.64
    residuals = [initial_res]
    
    for i in 2:50
        # 典型的な収束パターンを模擬
        if i < 10
            res = residuals[end] * (0.8 + 0.2 * rand())
        elseif i < 30
            res = residuals[end] * (0.7 + 0.3 * rand())
        else
            res = residuals[end] * (0.5 + 0.1 * rand())
        end
        push!(residuals, res)
        add_residual!(conv_pbicg, i, res)
    end
    
    # テストデータ2: SOR
    conv_sor = ConvergenceData("sor", "")
    
    # SORの典型的な線形収束
    sor_res = 5.0
    for i in 1:100
        add_residual!(conv_sor, i, sor_res)
        sor_res *= 0.95  # 一定の収束率
    end
    
    # テストデータ3: Jacobi
    conv_jacobi = ConvergenceData("jacobi", "")
    
    # Jacobiの遅い収束
    jac_res = 8.0
    for i in 1:150
        add_residual!(conv_jacobi, i, jac_res)
        jac_res *= 0.98  # 遅い収束
    end
    
    # 個別プロット
    println("個別プロット生成中...")
    plot_convergence_curve(conv_pbicg, "test_pbicgstab.png", target_tol=1.0e-6)
    plot_convergence_curve(conv_sor, "test_sor.png", target_tol=1.0e-6)
    plot_convergence_curve(conv_jacobi, "test_jacobi.png", target_tol=1.0e-6)
    
    # 比較プロット
    println("比較プロット生成中...")
    compare_convergence([conv_pbicg, conv_sor, conv_jacobi], "test_comparison.png", target_tol=1.0e-6)
    
    # CSV出力
    println("CSV出力中...")
    export_convergence_csv(conv_pbicg, "test_pbicgstab.csv")
    export_convergence_csv(conv_sor, "test_sor.csv")
    export_convergence_csv(conv_jacobi, "test_jacobi.csv")
    
    # 収束情報表示
    println("\n=== 収束情報 ===")
    for (name, data) in [("PBiCGSTAB", conv_pbicg), ("SOR", conv_sor), ("Jacobi", conv_jacobi)]
        info = get_convergence_info(data)
        println("\n$name:")
        println("  反復回数: $(info["iterations"])")
        println("  初期残差: $(@sprintf("%.4E", info["initial_residual"]))")
        println("  最終残差: $(@sprintf("%.4E", info["final_residual"]))")
        println("  収束率: $(@sprintf("%.4E", info["convergence_rate"]))")
        println("  桁数低下: $(@sprintf("%.2f", info["reduction_factor"]))")
    end
    
    println("\n=== テスト完了 ===")
    println("生成ファイル:")
    println("  - test_pbicgstab.png")
    println("  - test_sor.png") 
    println("  - test_jacobi.png")
    println("  - test_comparison.png")
    println("  - test_pbicgstab.csv")
    println("  - test_sor.csv")
    println("  - test_jacobi.csv")
end

# テスト実行
test_convergence_history()