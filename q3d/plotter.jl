using Plots


#=
@brief XZ断面（内部セル）
@param [in]     ｄ      解ベクトル
@param [in]     SZ     配列長
@param [in]     fname  ファイル名
=#
function plot_slice(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[2],2)
    s = d[2:SZ[1]-1,j,2:SZ[3]-1]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end


#=
@brief XZ断面(全セル)
@param [in]     ｄ      解ベクトル
@param [in]     SZ     配列長
@param [in]     fname  ファイル名
=#
function plot_slice2(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[2],2)
    s = d[1:SZ[1],j,1:SZ[3]]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end

