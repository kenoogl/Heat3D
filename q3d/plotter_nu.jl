using Plots

"""
NonUniform格子対応のプロッタ関数
"""

#=
@brief XZ断面（内部セル）- NonUniform格子対応
@param [in] d      解ベクトル
@param [in] SZ     配列長
@param [in] ox     原点座標
@param [in] Δh     X,Y方向格子間隔（Uniform）
@param [in] Z      Z方向格子点座標（NonUniform）
@param [in] fname  ファイル名
=#
function plot_slice_nu(d::Array{Float64,3}, SZ, ox, Δh, Z::Vector{Float64}, fname)
    j = div(SZ[2], 2)
    s = d[2:SZ[1]-1, j, 2:SZ[3]-1]
    
    # X方向座標軸（内部セル中心）
    x_coords = [ox[1] + Δh[1] * (i - 1.5) for i in 2:SZ[1]-1]
    
    # Z方向座標軸（内部セル中心、NonUniform）
    z_coords = Z
    
    # 物理座標軸でプロット
    p = contour(z_coords, x_coords, s, 
                fill=true, 
                c=:thermal, 
                xlabel="Z-coordinate [physical]", 
                ylabel="X-coordinate [physical]", 
                title="XZ-section (Y=$(ox[2] + Δh[2]*(j-1.5)), NonUniform)", 
                size=(600, 600),
                aspect_ratio=:equal)
    
    savefig(p, fname)
    
    # 格子情報を出力
    println("NonUniform格子プロット情報:")
    @printf("  Z座標範囲: [%.6f, %.6f]\n", minimum(z_coords), maximum(z_coords))
    @printf("  X座標範囲: [%.6f, %.6f]\n", minimum(x_coords), maximum(x_coords))
    @printf("  Z格子点数: %d (NonUniform)\n", length(z_coords))
    @printf("  X格子点数: %d (Uniform)\n", length(x_coords))
    
    return p
end

#=
@brief XZ断面（全セル）- NonUniform格子対応
@param [in] d      解ベクトル
@param [in] SZ     配列長
@param [in] ox     原点座標
@param [in] Δh     X,Y方向格子間隔（Uniform）
@param [in] Z      Z方向格子点座標（NonUniform）
@param [in] fname  ファイル名
=#
function plot_slice2_nu(d::Array{Float64,3}, SZ, ox, Δh, Z::Vector{Float64}, fname)
    j = div(SZ[2], 2)
    s = d[1:SZ[1], j, 1:SZ[3]]
    
    # X方向座標軸（境界含む全セル）
    x_coords = [ox[1] + Δh[1] * (i - 1.5) for i in 1:SZ[1]]
    
    # Z方向座標軸（境界含む全セル、NonUniform）
    z_coords = Z
    
    # 物理座標軸でプロット
    p = contour(z_coords, x_coords, s, 
                fill=true, 
                c=:thermal, 
                xlabel="Z-coordinate [physical]", 
                ylabel="X-coordinate [physical]", 
                title="XZ-section full (Y=$(ox[2] + Δh[2]*(j-1.5)), NonUniform)", 
                size=(600, 600),
                aspect_ratio=:equal)
    
    savefig(p, fname)
    
    return p
end

#=
@brief 格子分布の可視化
@param [in] SZ     配列長
@param [in] ox     原点座標
@param [in] Δh     X,Y方向格子間隔（Uniform）
@param [in] Z      Z方向格子点座標（NonUniform）
@param [in] fname  ファイル名
=#
function plot_grid_distribution_nu(SZ, ox, Δh, Z::Vector{Float64}, fname)
    
    # X方向格子点（Uniform）
    x_coords = [ox[1] + Δh[1] * (i - 1.5) for i in 1:SZ[1]]
    
    # Z方向格子点（NonUniform）
    z_coords = Z
    
    # プロット作成
    p = plot(layout=(2,1), size=(800, 800))
    
    # 上段：格子点分布
    scatter!(p[1], z_coords, ones(length(z_coords)), 
             label="Z-direction (NonUniform)", 
             markersize=4, 
             markershape=:circle,
             xlabel="Z-coordinate",
             ylabel="",
             title="Grid Distribution: NonUniform in Z")
    
    scatter!(p[1], ones(length(x_coords)), x_coords, 
             label="X-direction (Uniform)", 
             markersize=4, 
             markershape=:square,
             xlabel="",
             ylabel="X-coordinate")
    
    # 下段：格子間隔分布
    z_spacings = diff(z_coords)
    x_spacings = diff(x_coords)
    
    plot!(p[2], 2:length(z_coords), z_spacings, 
          label="Z spacing (NonUniform)", 
          linewidth=2, 
          marker=:circle,
          xlabel="Grid Index",
          ylabel="Grid Spacing",
          title="Grid Spacing Distribution")
    
    plot!(p[2], 2:length(x_coords), x_spacings, 
          label="X spacing (Uniform)", 
          linewidth=2, 
          marker=:square)
    
    savefig(p, fname)
    
    # 統計情報出力
    println("格子分布統計:")
    @printf("  Z方向 - 最小間隔: %.6f, 最大間隔: %.6f, 比率: %.2f\n", 
            minimum(z_spacings), maximum(z_spacings), maximum(z_spacings)/minimum(z_spacings))
    @printf("  X方向 - 間隔: %.6f (均等)\n", x_spacings[1])
    
    return p
end

#=
@brief Cartesian格子との比較プロット
@param [in] d_nu     NonUniform格子での解
@param [in] d_cart   Cartesian格子での解
@param [in] SZ       配列長
@param [in] ox       原点座標
@param [in] Δh       格子間隔
@param [in] Z        Z方向格子点座標（NonUniform）
@param [in] fname    ファイル名
=#
function plot_comparison_nu_cart(d_nu::Array{Float64,3}, d_cart::Array{Float64,3}, 
                                 SZ, ox, Δh, Z::Vector{Float64}, fname)
    j = div(SZ[2], 2)
    
    # NonUniform格子データ
    s_nu = d_nu[2:SZ[1]-1, j, 2:SZ[3]-1]
    x_coords = [ox[1] + Δh[1] * (i - 1.5) for i in 2:SZ[1]-1]
    z_coords_nu = [0.5 * (Z[k] + Z[k+1]) for k in 2:SZ[3]-1]
    
    # Cartesian格子データ
    s_cart = d_cart[2:SZ[1]-1, j, 2:SZ[3]-1]
    z_coords_cart = [ox[3] + Δh[3] * (k - 1.5) for k in 2:SZ[3]-1]
    
    # 比較プロット
    p = plot(layout=(1,3), size=(1500, 500))
    
    # NonUniform格子
    contour!(p[1], z_coords_nu, x_coords, s_nu, 
             fill=true, c=:thermal, 
             title="NonUniform Grid",
             xlabel="Z-coordinate", ylabel="X-coordinate")
    
    # Cartesian格子
    contour!(p[2], z_coords_cart, x_coords, s_cart, 
             fill=true, c=:thermal, 
             title="Cartesian Grid",
             xlabel="Z-coordinate", ylabel="X-coordinate")
    
    # 差分
    # 補間が必要な場合は省略し、同一格子点数の場合のみ差分計算
    if size(s_nu) == size(s_cart)
        diff_data = s_nu - s_cart
        contour!(p[3], z_coords_nu, x_coords, diff_data, 
                 fill=true, c=:RdBu, 
                 title="Difference (NU - Cart)",
                 xlabel="Z-coordinate", ylabel="X-coordinate")
    else
        plot!(p[3], title="Size mismatch\nCannot compute difference")
    end
    
    savefig(p, fname)
    
    return p
end