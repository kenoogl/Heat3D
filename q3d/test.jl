using Printf

const zm0 = 0.0
const zm1 = 0.1
const zm2 = 0.25
const zm3 = 0.4
const zm4 = 0.55
const pg_dpth = 0.005
const s_dpth = 0.1

silicon_1 = Dict(
    "x0" => 0.1, "y0" => 0.1, "z0" => zm1,
    "Lx" => 1.0, "Ly" => 1.0, "Lz" => 0.1, "depth" => 0.1
)
silicon_2 = Dict(
    "x0" => 0.1, "y0" => 0.1, "z0" => zm2,
    "Lx" => 1.0, "Ly" => 1.0, "Lz" => 0.1, "depth" => 0.1
)
silicon_3 = Dict(
    "x0" => 0.1, "y0" => 0.1, "z0" => zm3,
    "Lx" => 1.0, "Ly" => 1.0, "Lz" => 0.1, "depth" => 0.1
)

geom = Dict{String, Float64}[]
push!(geom, silicon_1)
push!(geom, silicon_2)
push!(geom, silicon_3)

function is_in_rect(x::Vector{Float64}, b::Vector{Float64}, L::Vector{Float64})
    ((b[1] ≤ x[1] ≤ b[1] + L[1]) &&
     (b[2] ≤ x[2] ≤ b[2] + L[2]) &&
     (b[3] ≤ x[3] ≤ b[3] + L[3]) )
end

function is_in_silicon2(x::Vector{Float64})
    b = zeros(Float64, 3)
    L = zeros(Float64, 3)

    for i in 1:3
        b[1] = geom[i]["x0"]
        b[2] = geom[i]["y0"]
        b[3] = geom[i]["z0"]
        L[1] = geom[i]["Lx"]
        L[2] = geom[i]["Ly"]
        L[3] = geom[i]["Lz"]
        
        if is_in_rect(x, b, L) 
            #@printf(stdout, "%4.2f %4.2f %4.2f : %4.2f %4.2f %4.2f\n", b[1], b[2], b[3], L[1], L[2], L[3])
            return true
        end
    end
    return false
end

# mode==1 : uniform , ==2 : non-uniform fixed
function coordZ!(Z::Vector{Float64}, SZ, mode::Int64, b, dz)
    if mode==1 
        for k in 1:SZ[3]+1
            Z[k] = b[3] + (k-2)*dz
        end
    else
        if SZ[3] != 14
            return
        end

        # mode != 1, NZ[3]=14
        Z[1] = -0.05
        Z[2] = 0.0 #zm0
        Z[3] = 0.05
        Z[4] = 0.1 #zm1
        Z[5] = 0.2-pg_dpth
        Z[6] = 0.2
        Z[7] = 0.25 #zm2
        Z[8] = 0.35-pg_dpth
        Z[9] = 0.35
        Z[10]= 0.4 #zm3
        Z[11]= 0.5-pg_dpth
        Z[12]= 0.5
        Z[13]= 0.55 #zm4
        Z[14]= 0.6
        Z[15]= 0.65
    end

    for k in 1:SZ[3]+1
        @printf(stdout, "%3d %f\n", k, Z[k])
    end
end

function modelA(SZ, Z::Vector{Float64}, b, dx)
    c = zeros(Float64, 3)
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        c[1] = b[1] + (i-1.5)*dx[1]
        c[2] = b[2] + (j-1.5)*dx[2]
        c[3] = 0.5*(Z[k]+Z[k+1])

        in_silicon = is_in_silicon2(c)

        if in_silicon
            @printf(stdout, "%d %d %d\n", i,j,k)
        end
    end
end

function main(NN_inner::Int64)
    NN = NN_inner + 2
    dh::Float64 = 0.6 / NN_inner
    dx = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0)
    SZ = (NN, NN, NN)

    Z = zeros(Float64, SZ[3]+1)

    coordZ!(Z, SZ, 1, ox, dx[3])

    modelA(SZ, Z, ox, dx)
end

main(20)