using Plots

function chk_i(i, nx)
    if i<1
        i = 1
    end
    if i>nx
        i = nx
    end
    return i
end

function chk_j(j, ny)
    if j<1
        j = 1
    end
    if j>ny
        j = ny
    end
    return j
end

function find_i(x::Float64, x0, dx, nx)
    i = floor( Int32, (x-x0)/dx+1.5 )
    chk_i(i, nx)
    return i
end

function find_j(y::Float64, y0, dy, ny)
    j = floor( Int32, (y-y0)/dy+1.5 )
    chk_j(j, ny)
    return j
end

#=
@brief XZ断面（内部セル）
@param [in]     ｄ      解ベクトル
@param [in]     SZ     配列長
@param [in]     ox     原点座標
@param [in]     Δh     X,Y方向格子間隔（Uniform）
@param [in]     fname  ファイル名
=#
function plot_slice(d::Array{Float64,3}, SZ, ox, Δh, fname)
    j = div(SZ[2],2)
    s = d[2:SZ[1]-1,j,2:SZ[3]-1]
    y = ox[2] + Δh[2]*(j-1.5)

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$y (j=$j)", size=(600, 600))
    savefig(p, fname)
end


#=
@brief XZ断面(全セル)
@param [in]     ｄ      解ベクトル
@param [in]     SZ     配列長
@param [in]     ox     原点座標
@param [in]     Δh     X,Y方向格子間隔（Uniform）
@param [in]     fname  ファイル名
=#
function plot_slice2(d::Array{Float64,3}, SZ, ox, Δh, fname)
    j = div(SZ[2],2)
    s = d[1:SZ[1],j,1:SZ[3]]
    y = ox[2] + Δh[2]*(j-1.5)

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$y (j=$j)", size=(600, 600))
    savefig(p, fname)
end

#=
@brief XZ断面（内部セル）- NonUniform格子対応
@param [in] d      解ベクトル
@param [in] y      Y座標
@param [in] SZ     配列長
@param [in] ox     原点座標
@param [in] Δh     X,Y方向格子間隔（Uniform）
@param [in] Z      Z方向格子点座標（NonUniform）
@param [in] fname  ファイル名
=#
function plot_slice_nu(d::Array{Float64,3}, y, SZ, ox, Δh, Z::Vector{Float64}, fname, label::String="")
    j = find_j(y, ox[2], Δh[2], SZ[2])
    #j = div(SZ[2], 2)
    s = d[2:SZ[1]-1, j, 2:SZ[3]-1]
    #y = ox[2] + Δh[2]*(j-1.5)
    
    # X方向座標軸（内部セル中心）
    x_coords = [ox[1] + Δh[1] * (i - 1.5) for i in 2:SZ[1]-1]
    
    # Z方向座標軸（内部セル中心、NonUniform）
    z_coords = [Z[k] for k in 2:SZ[3]-1]
    
    # 物理座標軸でプロット
    p = contour(z_coords, x_coords, s, 
                fill=true, 
                c=:thermal, 
                xlabel="Z-coordinate [physical]", 
                ylabel="X-coordinate [physical]", 
                title="XZ-section (Y=$y, j=$j, NonUniform) $label", 
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
@param [in] y      Y座標
@param [in] SZ     配列長
@param [in] ox     原点座標
@param [in] Δh     X,Y方向格子間隔（Uniform）
@param [in] Z      Z方向格子点座標（NonUniform）
@param [in] fname  ファイル名
=#
function plot_slice2_nu(d::Array{Float64,3}, y, SZ, ox, Δh, Z::Vector{Float64}, fname, label::String="")
    j = find_j(y, ox[2], Δh[2], SZ[2])
    #j = div(SZ[2], 2)
    s = d[1:SZ[1], j, 1:SZ[3]]
    #y = ox[2] + Δh[2]*(j-1.5)
    
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
                title="XZ-section full (Y=$y, j=$j, NonUniform) $label", 
                size=(600, 600),
                aspect_ratio=:equal)
    
    savefig(p, fname)
    
    return p
end