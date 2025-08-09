using Printf
using Plots

# 簡易版のテスト用モジュール
module TestNonUniform
    using Printf
    
    function read_grid_file()
        # テスト用のNonUniform格子データを生成
        # 左側密集の格子を模擬
        numNodes = 25
        coord = zeros(Float64, numNodes)
        
        # 指数関数的な格子分布を生成（左側密集）
        for i in 1:numNodes
            t = (i - 1) / (numNodes - 1)
            coord[i] = (exp(2*t) - 1) / (exp(2) - 1)  # [0,1]区間での指数分布
        end
        
        println("テスト用NonUniform格子生成:")
        @printf("  格子点数: %d\n", numNodes)
        @printf("  座標範囲: [%.6f, %.6f]\n", minimum(coord), maximum(coord))
        @printf("  最小間隔: %.6f\n", minimum(diff(coord)))
        @printf("  最大間隔: %.6f\n", maximum(diff(coord)))
        
        return coord, numNodes
    end
    
    function exact_solution!(exact::Array{Float64,3}, SZ, ox, Δh, Z::Vector{Float64})
        # テスト用厳密解（sin-sin分布）
        for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
            x = ox[1] + Δh[1] * (i - 1.5)
            y = ox[2] + Δh[2] * (j - 1.5)
            z = 0.5 * (Z[k] + Z[k+1])  # NonUniform格子でのZ座標
            
            exact[i,j,k] = sin(π*x) * sin(π*y) * sin(π*z)
        end
    end
end

# テスト用プロッタ関数のみを使用
include("plotter_nu.jl")

function demo_nonuniform_visualization()
    println("=== NonUniform格子可視化デモ ===\n")
    
    # パラメータ設定
    NXY = NZ = 25
    MX = MY = MZ = NXY + 2
    dh = 1.0 / NXY
    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0)
    
    println("計算条件:")
    println("  格子数: ", SZ)
    println("  格子間隔: ", Δh)
    
    # NonUniform格子生成
    read_coord, numNodes = TestNonUniform.read_grid_file()
    
    # Z座標配列の構築
    Z = zeros(Float64, MZ)
    for k in 1:NZ
        Z[k+1] = read_coord[k]
    end
    Z[1] = 2*Z[2] - Z[3]      # 境界条件
    Z[MZ] = 2*Z[MZ-1] - Z[MZ-2]  # 境界条件
    
    println("\nZ座標配列構築完了:")
    @printf("  Z範囲: [%.6f, %.6f]\n", Z[1], Z[MZ])
    
    # 格子分布可視化
    println("\n1. 格子分布可視化...")
    plot_grid_distribution_nu(SZ, ox, Δh, Z, "demo_grid_dist.png")
    
    # 厳密解の計算
    println("\n2. 厳密解計算...")
    exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    TestNonUniform.exact_solution!(exact, SZ, ox, Δh, Z)
    
    # 従来法でのプロット（問題のあるバージョン）
    println("\n3. 従来法プロット（格子座標系）...")
    j = div(SZ[2], 2)
    s = exact[2:SZ[1]-1, j, 2:SZ[3]-1]
    p_old = contour(s, fill=true, c=:thermal, 
                   xlabel="Z-index", ylabel="X-index", 
                   title="従来法: 格子インデックス軸", 
                   size=(600, 600))
    savefig(p_old, "demo_exact_old.png")
    println("  保存: demo_exact_old.png")
    
    # 修正法でのプロット（物理座標系）
    println("\n4. 修正法プロット（物理座標系）...")
    plot_slice_nu(exact, SZ, ox, Δh, Z, "demo_exact_new.png")
    println("  保存: demo_exact_new.png")
    
    # Cartesian格子との比較
    println("\n5. Cartesian格子参考プロット...")
    exact_cart = zeros(Float64, SZ[1], SZ[2], SZ[3])
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1] * (i - 1.5)
        y = ox[2] + Δh[2] * (j - 1.5)
        z = ox[3] + Δh[3] * (k - 1.5)  # Cartesian格子でのZ座標
        
        exact_cart[i,j,k] = sin(π*x) * sin(π*y) * sin(π*z)
    end
    
    s_cart = exact_cart[2:SZ[1]-1, j, 2:SZ[3]-1]
    x_coords = [ox[1] + Δh[1] * (i - 1.5) for i in 2:SZ[1]-1]
    z_coords_cart = [ox[3] + Δh[3] * (k - 1.5) for k in 2:SZ[3]-1]
    
    p_cart = contour(z_coords_cart, x_coords, s_cart, 
                    fill=true, c=:thermal, 
                    xlabel="Z-coordinate", ylabel="X-coordinate", 
                    title="Cartesian参考", 
                    size=(600, 600),
                    aspect_ratio=:equal)
    savefig(p_cart, "demo_exact_cartesian.png")
    println("  保存: demo_exact_cartesian.png")
    
    # 詳細比較プロット
    println("\n6. 詳細比較プロット...")
    plot_comparison_nu_cart(exact, exact_cart, SZ, ox, Δh, Z, "demo_comparison.png")
    println("  保存: demo_comparison.png")
    
    println("\n=== デモ完了 ===")
    println("生成ファイル:")
    for fname in ["demo_grid_dist.png", "demo_exact_old.png", "demo_exact_new.png", 
                  "demo_exact_cartesian.png", "demo_comparison.png"]
        if isfile(fname)
            println("  ✓ ", fname)
        else
            println("  ✗ ", fname)
        end
    end
    
    println("\n問題点と修正点:")
    println("  ❌ 従来法: 格子インデックスを軸として使用 → 物理的に不正確")
    println("  ✅ 修正法: 実際のZ座標値を軸として使用 → 物理的に正確")
    println("  ✅ アスペクト比を1:1に設定 → 幾何学的に正確")
end

# デモ実行
demo_nonuniform_visualization()