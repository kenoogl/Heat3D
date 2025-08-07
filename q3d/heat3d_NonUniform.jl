module NonUniform
export exact_solution!

# Harmonic mean
# @param a left value
# @param b right value
# @param ma mask for left
# @param mb mask for right
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))

#=
@brief 厳密解
@param [in]     e      解ベクトル
@param [in]     SZ     配列長
@param [in]     Δh     セル幅
@param [in,out] Z      
=#
function exact_solution!(e::Array{Float64,3}, SZ, Δh, Z::Float64)
    r2 = sqrt(2.0)
    ox = (0.0, 0.0, 0.0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        z = Z[k]
        e[i,j,k] = sin(π*x)*sin(π*y) / sinh(π*r2) * ( sinh(r2*π*z)-sinh(π*r2*(z-1.0)) )
    end
end


end # end of module