# Build functions for geometric model
# Ver. 1.0  2025-07-23

using Printf
using LinearAlgebra
using Plots

const zm0 = 0.0
const zm1 = 0.1
const zm2 = 0.25
const zm3 = 0.4
const zm4 = 0.55
const zm5 = 0.6
const pg_dpth = 0.005
const s_dpth = 0.1
const d_ufill = 0.05
const r_bump = 0.03
const r_tsv = 0.02
const h_tsv = 0.1

#  k=SZ[3]                           // 15
#  +----- top            0.6  = zm5  // 14
#  Heat sink (depth=0.05)
#  +-----                0.55 = zm4  // 13
#  underfill (depth=0.05)
#  +-----                0.5         // 12
#  Silcon 3 (depth=0.1)              // 11
#  +-----                0.4 = zm3   // 10
#  underfill (depth=0.05)
#  +-----                0.35        // 9
#  Silcon 2 (depth=0.1)              // 8
#  +-----                0.25 = zm2  // 7
#  underfill (depth=0.05)
#  +-----                0.2         // 6
#  Silcon 1 (depth=0.1)              // 5
#  +-----                0.1 = zm1   // 4
#  underfill (depth=0.05)
#  +------               0.05        // 3
#  Subtrate (depth=0.05)
#  +------ bottom        0.0 = zm0   // 2
#  k=1                               // 1

# boundary box of each element [mm]
substrate = Dict(
    "x0" => 0.0, "y0" => 0.0, "z0" => zm0,
    "Lx" => 1.2, "Ly" => 1.2, "Lz" => 0.05, "mat_id" => 4
)
silicon_1 = Dict(
    "x0" => 0.1, "y0" => 0.1, "z0" => zm1,
    "Lx" => 1.0, "Ly" => 1.0, "Lz" => 0.1, "mat_id" => 2
)
silicon_2 = Dict(
    "x0" => 0.1, "y0" => 0.1, "z0" => zm2,
    "Lx" => 1.0, "Ly" => 1.0, "Lz" => 0.1, "mat_id" => 2
)
silicon_3 = Dict(
    "x0" => 0.1, "y0" => 0.1, "z0" => zm3,
    "Lx" => 1.0, "Ly" => 1.0, "Lz" => 0.1, "mat_id" => 2
)
heatsink = Dict(
    "x0" => 0.0, "y0" => 0.0, "z0" => zm4,
    "Lx" => 1.2, "Ly" => 1.2, "Lz" => 0.05, "mat_id" => 5
)

# プライオリティ順に登録
geom = Dict{String, Any}[]
push!(geom, silicon_1)
push!(geom, silicon_2)
push!(geom, silicon_3)
push!(geom, substrate)
push!(geom, heatsink)

# 物性値
# λ 熱伝導率 [W/mK]
# ρ 密度 [kg/m^3]
# C 比熱 [J/KgK]
# α 温度拡散率 [m^2/s]

# TSV : Yellow
cupper = Dict(
    "id" => 1, "λ" => 386.0, "ρ" => 8960.0, "C" => 383.0, "α" => 1.12e-4
)
# Silicon : Green
silicon = Dict(
    "id" => 2, "λ" => 149.0, "ρ" => 2330.0, "C" => 720.0, "α" => 8.88e-5
)
# bump ハンダ : Purple
solder = Dict(
    "id" => 3, "λ" => 50.0, "ρ" => 8500.0, "C" => 197.0, "α" => 2.99e-5
)
# PCB subtrate : Orange
FR4 = Dict(
    "id" => 4, "λ" => 0.4, "ρ" => 1850.0, "C" => 1000.0, "α" => 2.16e-7
)
# Heatsink : Blue
A1060 = Dict(
    "id" => 5, "λ" => 222.0, "ρ" => 2700.0, "C" => 921.0, "α" => 8.93e-5
)
# Underfill : Grey
Resin = Dict(
    "id" => 6, "λ" => 1.5, "ρ" => 2590.0, "C" => 1050.0, "α" => 5.52e-7
)
# Power grid, Silicon : Red
pwrsrc = Dict(
    "id" => 7, "λ" => 149.0, "ρ" => 2330.0, "C" => 720.0, "α" => 8.88e-5
)

mat = Dict{String, Any}[]
push!(mat, cupper)
push!(mat, silicon)
push!(mat, solder)
push!(mat, FR4)
push!(mat, A1060)
push!(mat, Resin)
push!(mat, pwrsrc)

# =======================================

function searchMat(m::Int64)
    for i in 1:length(mat)
        if mat[i]["id"] == m
            return i
        end
    end
    # if exit for-loop
    println("search material error")
    exit(0)
end


function FillTSV!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    r = r_tsv
    h = h_tsv

    for z in [zm1, zm2, zm3], y in [0.3, 0.5, 0.7, 0.9], x in [0.3, 0.5, 0.7, 0.9]
        is, ie = find_Ri(x, r, ox[1], Δh[1], SZ[1])
        js, je = find_Rj(y, r, ox[2], Δh[2], SZ[2])
        ks = find_k(Z, z, SZ[3])
        ke = find_k(Z, z+h, SZ[3])
        #@printf(stdout, "TSV : [%d - %d]\n",ks, ke)

        for k in ks:ke, j in js:je, i in is:ie
            xc = ox[1] + (i-1.5)*Δh[1]
            yc = ox[2] + (j-1.5)*Δh[2]
            zc = 0.5*(Z[k]+Z[k+1])
            rx = xc - x
            ry = yc - y
            d = rx * rx + ry * ry
            if d ≤ r*r
                if 0 == ID[i,j,k]
                    ID[i,j,k] = cupper["id"]
                end
            end
        end
    end
end


function FillSolder!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    r = r_bump # ball radius

    for z in [zm1-d_ufill, zm2-d_ufill, zm3-d_ufill, zm4-d_ufill], y in [0.3, 0.5, 0.7, 0.9], x in [0.3, 0.5, 0.7, 0.9]
        is, ie = find_Ri(x, r, ox[1], Δh[1], SZ[1])
        js, je = find_Rj(y, r, ox[2], Δh[2], SZ[2])
        ks = find_k(Z, z, SZ[3])
        ke = find_k(Z, z+d_ufill, SZ[3])
        #@printf(stdout, "BUMP: [%d %d %d - %d %d %d]\n", is, js, ks, ie, je, ke)

        for k in ks:ke, j in js:je, i in is:ie
            xc = ox[1] + (i-1.5)*Δh[1]
            yc = ox[2] + (j-1.5)*Δh[2]
            #zc = 0.5*(Z[k]+Z[k+1])
            rx = xc - x
            ry = yc - y
            d = rx * rx + ry * ry
            if d ≤ r*r
                if 0 == ID[i,j,k]
                    ID[i,j,k] = solder["id"]
                end
            end
        end
    end
end


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

function find_Ri(x::Float64, r, x0, dx, nx)
    is = floor( Int64, (x-r-x0)/dx+1.5 )
    ie = floor( Int64, (x+r-x0)/dx+1.5 )
    chk_i(is, nx)
    chk_i(ie, nx)
    return is, ie
end

function find_Rj(y::Float64, r, y0, dy, ny)
    js = floor( Int64, (y-r-y0)/dy+1.5 )
    je = floor( Int64, (y+r-y0)/dy+1.5 )
    chk_i(js, ny)
    chk_i(je, ny)
    return js, je
end

function find_i(x::Float64, x0, dx, nx)
    i = floor( Int64, (x-x0)/dx+1.5 )
    chk_i(i, nx)
    return i
end

function find_j(y::Float64, y0, dy, ny)
    j = floor( Int64, (y-y0)/dy+1.5 )
    chk_j(j, ny)
    return j
end

function find_k(Z::Vector{Float64}, zc, nz)
    if zc<Z[1] || zc>Z[nz+1]
        println("out of scope in Z : find_z()")
        exit()
    end

    for k in 1:nz
        if Z[k] < zc ≤ Z[k+1]
            return k
        end
    end
end


# ジオメトリのbboxを計算
function find_index(b, L, ox, Δh, SZ, Z::Vector{Float64})
    st = zeros(Int64, 3)
    ed = zeros(Int64, 3)

    st[1] = find_i(b[1], ox[1], Δh[1], SZ[1])
    st[2] = find_i(b[2], ox[2], Δh[2], SZ[2])
    st[3] = find_k(Z, b[3], SZ[3])
    ed[1] = find_i(b[1]+L[1], ox[1], Δh[1], SZ[1])
    ed[2] = find_i(b[2]+L[2], ox[2], Δh[2], SZ[2])
    ed[3] = find_k(Z, b[3]+L[3], SZ[3])
    return st, ed
end


# ジオメトリのbboxをフィル
function FillPlate!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    b = zeros(Float64, 3)
    L = zeros(Float64, 3)

    for m in 1:5
        b[1] = geom[m]["x0"]
        b[2] = geom[m]["y0"]
        b[3] = geom[m]["z0"]
        L[1] = geom[m]["Lx"]
        L[2] = geom[m]["Ly"]
        L[3] = geom[m]["Lz"]
        st, ed = find_index(b, L, ox, Δh, SZ, Z)
        #@printf(stdout, "SILICON : [%d %d %d] - [%d %d %d]\n",st[1],st[2],st[3], ed[1],ed[2], ed[3])
        for k in st[3]:ed[3], j in st[2]:ed[2], i in st[1]:ed[1]
            l = searchMat( geom[m]["mat_id"] )
            if 0 == ID[i,j,k]
                ID[i,j,k] = mat[l]["id"]
            end
        end
    end
end


# 厚さ5µの領域
function FillPowerGrid!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64}, lx=0.2, ly=0.2)
    s = s_dpth-pg_dpth
    for z in [zm1+s, zm2+s, zm3+s], y in [0.3, 0.7], x in [0.3, 0.7]
        b = [x, y, z]
        L = [lx, ly, pg_dpth]
        st, ed = find_index(b, L, ox, Δh, SZ, Z)
        #@printf(stdout, "PG : [%d %d %d] - [%d %d %d]\n",st[1],st[2],st[3], ed[1],ed[2], ed[3])

        for k in st[3]:ed[3], j in st[2]:ed[2], i in st[1]:ed[1]
            if 0 == ID[i,j,k]
                ID[i,j,k] = pwrsrc["id"]
            end
        end
    end
end 


function FillResin!(ID::Array{Int64,3}, SZ)
    for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
        if 0 == ID[i,j,k]
            ID[i,j,k] = Resin["id"]
        end
    end
end


# Z軸座標の生成
# mode==1のとき等間隔、Δz=(top-bottom)/(SZ[3]-2)
function genZ!(Z::Vector{Float64}, SZ, mode::Int64, b, dz::Float64)
    if mode==1 
        for k in 1:SZ[3]+1
            Z[k] = b[3] + (k-2)*dz
        end
    else
        # mode != 1, NZ[3]=15
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
        Z[14]= 0.6  #zm5
        Z[15]= 0.65
    end
    
    #=
    for k in 1:SZ[3]+1
        @printf(stdout, "%3d : %6.3f\n", k, Z[k])
    end
    =#
end


# ここでλは温度拡散率
function setLambda!(λ::Array{Float64,3}, SZ, ID::Array{Int64,3})
    for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
        t = ID[i,j,k]
        if t == pwrsrc["id"]
            λ[i,j,k] = pwrsrc["α"]
        elseif t == cupper["id"]
            λ[i,j,k] = cupper["α"]
        elseif t == silicon["id"]
            λ[i,j,k] = silicon["α"]
        elseif t == FR4["id"]
            λ[i,j,k] = FR4["α"]
        elseif t == A1060["id"]
            λ[i,j,k] = A1060["α"]
        elseif t == solder["id"]
            λ[i,j,k] = solder["α"]
        else
            λ[i,j,k] = Resin["α"]
        end
    end
end


function fillID!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    FillPowerGrid!(ID, ox, Δh, SZ, Z)
    FillTSV!(ID, ox, Δh, SZ, Z)
    FillPlate!(ID, ox, Δh, SZ, Z)
    FillSolder!(ID, ox, Δh, SZ, Z)
    FillResin!(ID, SZ)
end

# @param mode (1--Uniform in Z-dir, 2--Non-uniform in Z), while Uniform for X&Y
# @param NXY  Number of inner cells for X&Y dir.
# @param NZ   Number of inner cells for Z dir.
function model_test(mode::Int64, NXY::Int64, NZ::Int64=14)
    MX = MY = NXY + 2  # Number of CVs including boundaries
    if mode==2
        NZ = 14
    end
    MZ = NZ + 2

    dx::Float64 = 1.2 / NXY
    dy::Float64 = 1.2 / NXY
    dz::Float64 = zm5 / NZ
    Δh = (dx, dy, dz) 
    ox = (0.0, 0.0, 0.0)
    SZ = (MX, MY, MZ)
    println(SZ)
    println(Δh)
    
    ID = zeros(Int64, SZ[1], SZ[2], SZ[3])
    Z  = zeros(Float64, SZ[3]+1)
    genZ!(Z, SZ, mode, ox, Δh[3])

    @time fillID!(ID, ox, Δh, SZ, Z)


    #plot_slice(ID, SZ, "id.png")
    #=
    id_xy(ID, 37, "id_k=37.png")
    id_xy(ID, 39, "id_k=39.png")
    id_xy(ID, 36, "id_k=36.png")
    id_xy(ID, 38, "id_k=38.png")
    id_xy(ID, 41, "id_k=41.png")
    id_xy(ID, 42, "id_k=42.png")
    id_xy(ID, 50, "id_k=50.png")
    id_xy(ID, 51, "id_k=51.png")
    id_xy(ID, 52, "id_k=52.png")

    id_yz(ID, 59, "id_i=59.png")
    id_yz(ID, 60, "id_i=60.png")
    id_yz(ID, 61, "id_i=61.png")
    id_yz(ID, 63, "id_i=63.png")
    =#
    id_xy(ID, 0.2, SZ, Z, "id_z=0.2.png")
    id_yz(ID, 0.3, ox, Δh, SZ, "id_x=0.3.png")
 
end

fixed_colors = [:yellow, :green, :purple, :orange, :blue, :gray, :red]
Fcolor = palette(fixed_colors)

function id_xy(d::Array{Int64,3}, z, SZ, Z, fname)
    k = find_k(Z, z, SZ[3])
    p = heatmap( d[:, :, k], 
        clims=(1,length(Fcolor)), 
        title="ID z_index=$k",
        c = Fcolor,
        colorbar=false,
        size = (600,600) )
    savefig(p, fname)
end

function id_yz(d::Array{Int64,3}, x, ox, Δh, SZ, fname)
    i = find_i(x, ox[1], Δh[1], SZ[1])
    p = heatmap( d[i, :, :], 
        clims=(1,length(Fcolor)), 
        title="ID x_index=$i",
        c = Fcolor,
        colorbar=false,
        size = (300,600) )
    savefig(p, fname)
end

function plot_slice(d::Array{Int64,3}, SZ, fname)
    j = 16 #div(SZ[3],2)
    s = d[1:SZ[1],j,1:SZ[3]]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
end

model_test(1,240,120)
#model_test(1,120,60)
#model_test(1,480,240)