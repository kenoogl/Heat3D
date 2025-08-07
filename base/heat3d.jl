using Printf
using Plots
using LinearAlgebra

# Harmonic mean
# @param a left value
# @param b right value
# @param ma mask for left
# @param mb mask for right
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))


function residual(p::Array{Float64,3}, SZ, 
                  λ::Array{Float64,3}, 
                  b::Array{Float64,3},
                  m::Array{Float64,3}, Δh, ω::Float64)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
        dd = (1.0-m0) + (ae + aw + an + as + at + ab)*m0
        ss = ( ae * p[i+1,j  ,k  ] + aw * p[i-1,j  ,k  ]
             + an * p[i  ,j+1,k  ] + as * p[i  ,j-1,k  ]
             + at * p[i  ,j  ,k+1] + ab * p[i  ,j  ,k-1] )
        dp = (((ss-bb)/dd - pp)) * m0
        r = (dd + ω*(aw+as+ab))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end

function sor!(p::Array{Float64,3}, SZ, 
              λ::Array{Float64,3}, 
              b::Array{Float64,3},
              m::Array{Float64,3}, Δh, ω::Float64)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = p[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
        aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
        an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
        as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
        at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
        ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
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

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
end


function rbsor_core!(p::Array{Float64,3}, SZ, 
                     λ::Array{Float64,3}, 
                     b::Array{Float64,3},
                     m::Array{Float64,3}, Δh, ω::Float64, color::Int64)
    res::Float64 = 0.0
    dx0 = Δh[1]
    dy0 = Δh[2]
    dz0 = Δh[3]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dz2 = 1.0 / (dz0*dz0)

    for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            pp = p[i,j,k]
            bb = b[i,j,k]
            λ0 = λ[i,j,k]
            m0 = m[i,j,k]
            ae = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2
            aw = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2
            an = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2
            as = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2
            at = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dz2
            ab = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dz2
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
    end

    return res
end


function rbsor!(θ::Array{Float64,3}, SZ, 
                λ::Array{Float64,3}, 
                b::Array{Float64,3},
                mask::Array{Float64,3}, Δh, ω::Float64)
    res::Float64 = 0.0
    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        r = rbsor_core!(θ, SZ, λ, b, mask, Δh, ω, c)
        res += r
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
    
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
    p = contour(s, fill=true, c=:thermal, xlabel="X-axis", ylabel="Z-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end

function plot_slice2(d::Array{Float64,3}, SZ, fname)
    j = div(SZ[3],2)
    s = d[1:SZ[1],j,1:SZ[3]]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="X-axis", ylabel="Z-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end


function calRHS!(b::Array{Float64,3}, SZ)
    b .= 0.0
end


function solve!(θ, SZ, λ, b, mask, Δh, ω, tol, F)

    res0 = residual(θ, SZ, λ, b, mask, Δh, ω)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:Itr_max
        res = sor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        #res = rbsor!(θ, SZ, λ, b, mask, Δh, ω) / res0
        #println(n, " ", res)
        @printf(F, "%10d %24.14E\n", n, res)
        if res < tol
            println("Converge at ", n)
            return
        end
    end
end

function main(SZ, O, Δh, ω, tol, exact, θ)

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    setMaterial!(λ, SZ)
    #plot_slice2(λ, SZ, "lambda.png")

    boundary_condition!(θ, SZ, O, Δh)
    #plot_slice2(θ, SZ, "p0.png")

    b = zeros(Float64, SZ[1], SZ[2], SZ[3])
    #calRHS!(b, SZ)

    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    setMask!(mask, SZ)
    #plot_slice2(mask, SZ, "mask.png")

    F = open("res.txt", "w")
    solve!(θ, SZ, λ, b, mask, Δh, ω, tol, F)
    nrm = sqrt(norm(θ-exact,2))
    @printf(F, "Sum of norm = %24.14E\n", nrm)
    close(F)
end

#=
# Visio用のsphフォーマット出力
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

# @param NN_inner::Int64  Number of 
# @param ω::Float64  coef. of Accel.
function diff3d(NN_inner::Int64, ω::Float64=1.0)
    NN = NN_inner + 2  # Number of CVs including boundaries

    global Itr_max = 5000

    dh::Float64 = 1.0 / NN_inner
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0)
    SZ = (NN, NN, NN)
    tol = 1.0e-8

    exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    exact_solution!(exact, SZ, ox, Δh)
    plot_slice(exact, SZ, "exact.png")
    #writeSPH(SZ, ox, Δh, 0, 0.0, exact)

    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])

    @time main(SZ, ox, Δh, ω, tol, exact, θ)

    plot_slice2(θ, SZ, "p.png")
    plot_slice2(θ-exact, SZ, "diff.png")
end

diff3d(1) # just compile　JITコンパイルを行うためにパラメータは1
diff3d(25) # ここで本実行し、計測