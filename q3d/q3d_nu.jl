using Printf
using ProgressMeter
using Plots

# Harmonic mean
# @param a left value
# @param b right value
# @param ma mask for left
# @param mb mask for right
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))


function residual(p::Array{Float64,3}, SZ, λ::Array{Float64,3}, b::Array{Float64,3},
             m::Array{Float64,3}, Δh, ω::Float64, k::Int64)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    #dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    #dz1 = 1.0 / dz0

    for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        #at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz1
        #ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz1
        #dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        dd = (1.0-m0) + (ae + aw + an + as)*m0
        ss =   ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
            # + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1]
        dp = (((ss-bb)/dd - pp)) * m0
        #r = (dd + ω*(aw+as+ab))*dp / ω
        r = (dd + ω*(aw+as))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2))
end

function sor2d!(p::Array{Float64,3}, SZ, λ::Array{Float64,3}, b::Array{Float64,3},
             m::Array{Float64,3}, Δh, ω::Float64, k::Int64)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz1 = 1.0 / dz0
    for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz1
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz1
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        pn = pp + ω * dp
        p[i,j,k] = pn
        r = (dd + ω*(aw+as+ab))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2))
end


function boundary_condition!(p::Array{Float64,3}, SZ, O, Δh)

    for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = O[1] + Δh[1]*(i-1.5)
        y = O[2] + Δh[2]*(j-1.5)
        a = sin(π*x)*sin(π*y)
        p[i,j,1    ] = a
        p[i,j,SZ[3]] = a
    end

    for k in 2:SZ[3]-1, i in 2:SZ[1]-1
        p[i,1    ,k] = 0.0
        p[i,SZ[2],k] = 0.0
    end

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        p[1    ,j,k] = 0.0
        p[SZ[1],j,k] = 0.0
    end
end


function setMask!(m::Array{Float64,3}, SZ)

    for j in 1:SZ[2], i in 1:SZ[1]
        m[i,j,1    ] = 0.0
        m[i,j,SZ[3]] = 0.0
    end

    for k in 1:SZ[3], i in 1:SZ[1]
        m[i,1    ,k] = 0.0
        m[i,SZ[2],k] = 0.0
    end

    for k in 1:SZ[3], j in 1:SZ[2]
        m[1    ,j,k] = 0.0
        m[SZ[1],j,k] = 0.0
    end
end



function setMaterial!(λ::Array{Float64,3}, SZ)
    λ .= 1.0
#=
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ[i,j,k] = 1.0
    end
=#
end


function exact_solution!(e::Array{Float64,3}, SZ, O, Δh)
    r2 = sqrt(2.0)
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = O[1] + Δh[1]*(i-1.5)
        y = O[2] + Δh[2]*(j-1.5)
        z = O[3] + Δh[3]*(k-1.5)
        e[i,j,k] = sin(π*x)*sin(π*y) / sinh(π*r2) * ( sinh(r2*π*z)-sinh(π*r2*(z-1.0)) )
    end
end


function plot_slice(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[3],2)
    s = d[2:SZ[1]-1,j,2:SZ[3]-1]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end

function plot_slice2(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[3],2)
    s = d[1:SZ[1],j,1:SZ[3]]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end

#=
function calRHS!(b::Array{Float64,3}, SZ, θ::Array{Float64,3}, 
                 λ::Array{Float64,3}, m::Array{Float64,3}, Δh, k::Int64)
    dz0 = Δh[3]
    dz1 = 1.0 / dz0
    sb = 0.0

    for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        pp = θ[i,j,k]
        qr = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz1 * (θ[i  ,j  ,k+1] - pp)
        ql = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz1 * (pp - θ[i  ,j  ,k-1])
        Δq = qr - ql
        Δb = Δq - b[i,j,k]
        b[i,j,k] = Δq
        sb += Δb*Δb
    end

    return sqrt(sb)/((SZ[1]-2)*(SZ[2]-2))
end

function calRHS2!(b::Array{Float64,3}, SZ, θ::Array{Float64,3}, Δh, k::Int64)
    sb = 0.0
    dz0 = Δh[3]
    dz1 = 1.0 / dz0

    for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        Δθ = θ[i,j,k] - 0.5*(θ[i,j,k+1]+θ[i,j,k-1])
        b[i,j,k] = -Δθ * dz1
        sb += Δθ*Δθ
    end

    return sqrt(sb)/((SZ[1]-2)*(SZ[2]-2))
end
=#

function plot_temp(θ::Array{Float64,3}, e::Array{Float64,3}, SZ, Δh)
    
    x = [Δh[3]*(m-1.5) for m in 1:SZ[3]]
    x[1] = 0.0
    x[SZ[3]] = 1.0
    y1= θ[div(SZ[1],2), div(SZ[2],2), :]
    y2= e[div(SZ[1],2), div(SZ[2],2), :]


    p=plot(x, y1, linewidth=2, label="Quasi3D", 
           xlabel="Z", ylabel="θ", title="Temerature distribution",
           marker=:circle, markersize=5,
           legend=:top )
    plot!(x, y2, linewidth=2, label="Exact", 
           xlabel="Z", ylabel="θ", title="Temerature distribution",
           marker=:triangle, markersize=5 )
    display(p)
end


function solve_slice!(θ::Array{Float64,3}, SZ, λ::Array{Float64,3}, b::Array{Float64,3},
             mask::Array{Float64,3}, Δh, ω, k::Int64, F)

    res0 = residual(θ, SZ, λ, b, mask, Δh, ω, k)
    #if res0==0.0
    #    res0 = 1.0
    #end

    @printf(F, "k=%3d res0=%24.14E\n", k, res0)

    for n in 1:Itr_inner
        res = sor2d!(θ, SZ, λ, b, mask, Δh, ω, k) / res0
        @printf(F, "%10d %24.14E\n", n, res)
        if res < tol_inner
            return n, res
        end
    end

    return n, res
end


function solve!(θ, SZ, λ, b, mask, Δh, ω, F)
    for loop in 1:Itr_outer

        @printf(stdout, "Loop=%4d\n", loop)

        sum_r = 0.0
        for k in 2:SZ[3]-1
            #sb = calRHS!(b, SZ, θ, λ, mask, Δh, k)
            #sb = calRHS2!(b, SZ, θ, Δh, k)
            #sum_b += sb
            n, r = solve_slice!(θ, SZ, λ, b, mask, Δh, ω, k, F)
            sum_r += r
            @printf(stdout, "     k=%3d Itr=%5d r(k)= %18.8e sumR= %18.8e\n", k, n, r, sum_r)
        end

        if sum_r < tol_outer
            return
        end 
    end
end

function main(SZ, O, Δh, ω, mode, exact)

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    setMaterial!(λ, SZ)
    #plot_slice2(λ, SZ, "lambda.png")

    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    boundary_condition!(θ, SZ, O, Δh)
    #plot_slice2(θ, SZ, "p0.png")

    b = zeros(Float64, SZ[1], SZ[2], SZ[3])

    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    setMask!(mask, SZ)
    #plot_slice2(mask, SZ, "mask.png")

    F = open("qres.txt", "w")
    solve!(θ, SZ, λ, b, mask, Δh, ω, F)
    close(F)
    
    plot_slice(θ, SZ, "qp.png")
    plot_slice(θ-exact, SZ, "qdiff.png")
    plot_temp(θ, exact, SZ, Δh)
end

#=
function writeSPH(size, org, pch, step, time, var)
    ccall((:write_sph_d, "./iosph.so"), Nothing,
    (Ref{Int64},    # size
     Ref{Float64},  # org
     Ref{Float64},  # pch
     Ref{Int64},    # step
     Ref{Float64},  # time
     Ref{Float64}), # var
     size, org, pch, step, time, var)
end
=#

# Z軸座標の生成
# mode==1のとき等間隔、Δz=(top-bottom)/(SZ[3]-2)
function genZ!(Z::Vector{Float64}, SZ, mode::Int64, ox, dz::Float64)
    if mode==1 
        for k in 1:SZ[3]+1
            Z[k] = ox[3] + (k-2)*dz
        end
    else
        for k in 1:SZ[3]+1
            Z[k] = ox[3] + (k-2)*dz
        end
    end
    
    #=
    for k in 1:SZ[3]+1
        @printf(stdout, "%3d : %6.3f\n", k, Z[k])
    end
    =#
end

#=
@param NN_inner::Int64  Number of 
@param ω::Float64  coef. of Accel.
@param mode  1-uniform, 2-non-uniform(uniform dist), 3-non-uniform
1では、直交等間隔を計算、NN_inner=NN_Zかつmode=1
2ではZ方向の不等間隔コードを等間隔データでテスト
3では不等間隔データでテスト
=#

function diffq3d(NN_inner::Int64, NN_Z::Int64, mode::Int64=1, ω::Float64=1.0)
    NN = NN_inner + 2  # Number of CVs including boundaries
    NZ = NN_Z + 2

    global Itr_outer = 1000
    global Itr_inner = 2000
    global tol_inner = 1.0e-8
    global tol_outer = 1.0e-8

    if mode==1
        dh::Float64 = 1.0 / NN_inner
        SZ = (NN, NN, NN)
    else if mode==2
        SZ = (NN, NN, NZ)
    else if mode==3
        SZ = (NN, NN, NZ)
    else
        println("mode error")
        return
    end
    
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0)
    Z  = zeros(Float64, SZ[3]+1)

    if mode>=2
        genZ!(Z, SZ, mode, ox, Δh[3])
    end

    exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    if mode==1
        exact_solution!(exact, SZ, ox, Δh)
    else
        exact_solution!(exact, SZ, ox, Δh, Z)
    end
    boundary_condition!(exact, SZ, ox, Δh)
    #plot_slice(exact, SZ, "exact.png")
    #writeSPH(SZ, ox, Δh, 0, 0.0, exact)

    @time main(SZ, ox, Δh, ω, mode, exact)
    println(" ")
end

diffq3d(1, 1, 1) # just compile
diffq3d(25, 25, 1)