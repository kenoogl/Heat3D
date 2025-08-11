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

# geomに登録
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
    if mode==1
        if zc<Z[1] || zc>Z[nz+1]
            println("out of scope in Z : find_z()")
            println(zc, Z)
            exit()
        end

        for k in 1:nz
            if Z[k] ≤ zc < Z[k+1]
                return k
            end
        end
    else
        if zc<Z[1] || zc>Z[nz]
            println("out of scope in Z : find_z()")
            println(zc, Z)
            exit()
        end

        for k in 1:nz-1
            if Z[k] ≤ zc < Z[k+1]
                return k
            end
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


"""
    is_included_rect(a1, a2, b1, b2) -> Bool

セルA(a1,a2) が直方体のジオメトリ領域B(b1,b2) に 50%以上含まれている場合 true を返す。
a1,a2,b1,b2 は対角をなす2点の座標 (x,y,z) タプル。
"""
function is_included_rect(a1, a2, b1, b2)
    # Aの体積
    volA = abs((a2[1] - a1[1]) * (a2[2] - a1[2]) * (a2[3] - a1[3]))

    # 重なり体積
    axlo, aylo, azlo = min.(a1, a2)
    axhi, ayhi, azhi = max.(a1, a2)
    bxlo, bylo, bzlo = min.(b1, b2)
    bxhi, byhi, bzhi = max.(b1, b2)

    ox = max(0.0, min(axhi, bxhi) - max(axlo, bxlo))
    oy = max(0.0, min(ayhi, byhi) - max(aylo, bylo))
    oz = max(0.0, min(azhi, bzhi) - max(azlo, bzlo))

    overlap_vol = ox * oy * oz

    return overlap_vol >= 0.5 * volA
end

"""
    is_half_or_more_included_cuboid_in_cylinder(a1, a2, cyl_center, cyl_radius, cyl_zmin, cyl_zmax; samples=50)

直方体(a1,a2)の体積のうち、円柱（Z軸方向）に含まれる割合が50%以上ならtrue。
円柱は中心xy座標 `cyl_center`、半径 `cyl_radius`、
高さ範囲 [cyl_zmin, cyl_zmax] で指定。
"""
function is_included_cyl(a1, a2, cyl_ctr, cyl_r, cyl_zmin, cyl_zmax; samples=50)
    # 直方体の境界
    xlo, ylo, zlo = min.(a1, a2)
    xhi, yhi, zhi = max.(a1, a2)

    volA = (xhi - xlo) * (yhi - ylo) * (zhi - zlo)

    inside_count = 0
    total_count = 0

    for i in 1:samples, j in 1:samples, k in 1:samples
        # 小セルの中心座標
        x = xlo + (i - 0.5) * (xhi - xlo) / samples
        y = ylo + (j - 0.5) * (yhi - ylo) / samples
        z = zlo + (k - 0.5) * (zhi - zlo) / samples

        # 円柱内判定
        dx = x - cyl_ctr[1]
        dy = y - cyl_ctr[2]
        r2 = dx^2 + dy^2
        if r2 <= cyl_r^2 && cyl_zmin <= z <= cyl_zmax
            inside_count += 1
        end
        total_count += 1
    end

    overlap_vol = volA * inside_count / total_count
    return overlap_vol >= 0.5 * volA
end

"""
    is_half_or_more_included_cuboid_in_sphere(a1, a2, center, radius; samples=50) -> Bool

直方体 A (対角点 a1, a2) の体積のうち、球 (center, radius) に含まれる割合が50%以上なら true を返す。
評価は一様グリッドサンプリングで行い、samples^3 点で近似する（samplesを増やすと精度↑, 計算量↑）。

引数:
- a1, a2 :: (x,y,z) 直方体の対角点（順不同）
- center :: (cx,cy,cz) 球の中心
- radius :: 半径
- samples :: 各軸の分割数（既定 50）

注意:
- 直方体は軸平行（AABB）を仮定
- 早期判定: 直方体の8頂点がすべて球内なら即 true（完全包含）
"""
function is_included_sph(a1, a2, center, radius; samples::Int=50)
    # 直方体の境界（昇順にそろえる）
    xlo, ylo, zlo = min.(a1, a2)
    xhi, yhi, zhi = max.(a1, a2)

    # 体積
    volA = (xhi - xlo) * (yhi - ylo) * (zhi - zlo)
    volA ≤ 0 && return false  # 退避：変な入力

    cx, cy, cz = center
    r2 = radius^2

    # 早期 true 判定：8頂点がすべて球内なら完全包含
    corners = ((xlo,ylo,zlo),(xlo,ylo,zhi),(xlo,yhi,zlo),(xlo,yhi,zhi),
               (xhi,ylo,zlo),(xhi,ylo,zhi),(xhi,yhi,zlo),(xhi,yhi,zhi))
    all_inside = all(((x - cx)^2 + (y - cy)^2 + (z - cz)^2 ≤ r2) for (x,y,z) in corners)
    if all_inside
        return true
    end

    # サンプリング（セル中心）
    inside_count = 0
    total_count  = samples^3

    dx = (xhi - xlo) / samples
    dy = (yhi - ylo) / samples
    dz = (zhi - zlo) / samples

    x = xlo + dx/2
    for i in 1:samples
        y = ylo + dy/2
        for j in 1:samples
            z = zlo + dz/2
            for k in 1:samples
                # 点が球内か
                if (x - cx)^2 + (y - cy)^2 + (z - cz)^2 ≤ r2
                    inside_count += 1
                end
                z += dz
            end
            y += dy
        end
        x += dx
    end

    overlap_vol_est = volA * (inside_count / total_count)
    return overlap_vol_est ≥ 0.5 * volA
end

# ジオメトリのbboxをフィル
function FillPlate!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    b = zeros(Float64, 3)
    L = zeros(Float64, 3)
    c1= zeros(Float64, 3)
    c2= zeros(Float64, 3)
    d1= zeros(Float64, 3)
    d2= zeros(Float64, 3)

    for m in 1:5
        b[1] = geom[m]["x0"]
        b[2] = geom[m]["y0"]
        b[3] = geom[m]["z0"]
        L[1] = geom[m]["Lx"]
        L[2] = geom[m]["Ly"]
        L[3] = geom[m]["Lz"]
        st, ed = find_index(b, L, ox, Δh, SZ, Z)
        #@printf(stdout, "SILICON : [%d %d %d] - [%d %d %d]\n",st[1],st[2],st[3], ed[1],ed[2], ed[3])
        d1 = b
        d2 .= b + L
        for k in st[3]-1:ed[3], j in st[2]-1:ed[2], i in st[1]-1:ed[1]
            c1[1] = ox[1] + Δh[1]*(i-1)
            c1[2] = ox[2] + Δh[2]*(j-1)
            c1[3] = Z[k]
            c2[1] = c1[1] + Δh[1]
            c2[2] = c1[2] + Δh[2]
            c2[3] = Z[k+1]
            #@printf(stdout, "cell : [%f %f %f] - [%f %f %f]\n",c1[1],c1[2],c1[3], c2[1],c2[2], c2[3])
            if is_included_rect(c1, c2, d1, d2) && (0 == ID[i,j,k])
                l = searchMat( geom[m]["mat_id"] )
                ID[i,j,k] = mat[l]["id"]
                #println("$i $j $k")
            end
        end
    end
end


# 厚さ5µの領域
function FillPowerGrid!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64}, lx=0.2, ly=0.2)
    c1= zeros(Float64, 3)
    c2= zeros(Float64, 3)
    d1= zeros(Float64, 3)
    d2= zeros(Float64, 3)
    s = s_dpth-pg_dpth
    pg_count = 0
    for z in [zm1+s, zm2+s, zm3+s], y in [0.3, 0.7], x in [0.3, 0.7]
        b = [x, y, z]
        L = [lx, ly, pg_dpth]
        st, ed = find_index(b, L, ox, Δh, SZ, Z)
        #@printf(stdout, "PG mode=%d: z=%.3f, range=[%.3f,%.3f] -> k=[%d,%d]\n", mode, z, z, z+pg_dpth, st[3], ed[3])
        d1 = b
        d2 .= b + L
        #@printf(stdout, "target : [%f %f %f] - [%f %f %f]\n",d1[1],d1[2],d1[3], d2[1],d2[2],d2[3])
        for k in st[3]-1:ed[3]+1, j in st[2]-1:ed[2]+1, i in st[1]-1:ed[1]+1
            c1[1] = ox[1] + Δh[1]*(i-1)
            c1[2] = ox[2] + Δh[2]*(j-1)
            c1[3] = Z[k]
            c2[1] = c1[1] + Δh[1]
            c2[2] = c1[2] + Δh[2]
            c2[3] = Z[k+1]
            if is_included_rect(c1, c2, d1, d2) && (0 == ID[i,j,k])
                ID[i,j,k] = pwrsrc["id"]
                pg_count += 1
                #println("$i $j $k")
            end
        end
    end
    #@printf(stdout, "Total power grid cells filled (mode=%d): %d\n", mode, pg_count)
end 


function FillResin!(ID::Array{Int64,3}, SZ)
    for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
        if 0 == ID[i,j,k]
            ID[i,j,k] = Resin["id"]
        end
    end
end


function FillTSV!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    c1= zeros(Float64, 3)
    c2= zeros(Float64, 3)
    cyl_ctr= zeros(Float64, 2)
    cyl_r = r_tsv

    for z in [zm1, zm2, zm3], y in [0.3, 0.5, 0.7, 0.9], x in [0.3, 0.5, 0.7, 0.9]
        cyl_ctr[1] = x
        cyl_ctr[2] = y
        cyl_zmin = z
        cyl_zmax = z + h_tsv
        is, ie = find_Ri(x, cyl_r, ox[1], Δh[1], SZ[1])
        js, je = find_Rj(y, cyl_r, ox[2], Δh[2], SZ[2])
        ks = find_k(Z, cyl_zmin, SZ[3])
        ke = find_k(Z, cyl_zmax, SZ[3])
        #@printf(stdout, "TSV : [%d - %d]\n",ks, ke)

        for k in ks-1:ke+1, j in js-1:je+1, i in is-1:ie+1
            c1[1] = ox[1] + Δh[1]*(i-1)
            c1[2] = ox[2] + Δh[2]*(j-1)
            c1[3] = Z[k]
            c2[1] = c1[1] + Δh[1]
            c2[2] = c1[2] + Δh[2]
            c2[3] = Z[k+1]
            if (0 == ID[i,j,k]) && is_included_cyl(c1, c2, cyl_ctr, cyl_r, cyl_zmin, cyl_zmax)
                ID[i,j,k] = cupper["id"]
            end
        end
    end
end


function FillSolder!(ID::Array{Int64,3}, ox, Δh, SZ, Z::Vector{Float64})
    r = r_bump # ball radius
    c1= zeros(Float64, 3)
    c2= zeros(Float64, 3)
    ctr= zeros(Float64, 3)
    dp = d_ufill*0.5

    for z in [zm1-dp, zm2-dp, zm3-dp, zm4-dp], y in [0.3, 0.5, 0.7, 0.9], x in [0.3, 0.5, 0.7, 0.9]
        ctr[1] = x
        ctr[2] = y
        ctr[3] = z
        is, ie = find_Ri(x, r, ox[1], Δh[1], SZ[1])
        js, je = find_Rj(y, r, ox[2], Δh[2], SZ[2])
        ks = find_k(Z, z-dp, SZ[3])
        ke = find_k(Z, z+dp, SZ[3])
        #@printf(stdout, "BUMP: [%d %d %d - %d %d %d]\n", is, js, ks, ie, je, ke)

        for k in ks-1:ke+1, j in js-1:je+1, i in is-1:ie+1
            c1[1] = ox[1] + Δh[1]*(i-1)
            c1[2] = ox[2] + Δh[2]*(j-1)
            c1[3] = Z[k]
            c2[1] = c1[1] + Δh[1]
            c2[2] = c1[2] + Δh[2]
            c2[3] = Z[k+1]
            if (0 == ID[i,j,k]) && is_included_sph(c1, c2, ctr, r)
                ID[i,j,k] = solder["id"]
            end
        end
    end
end


# Z軸座標の生成
# mode==1のとき等間隔、Δz=(top-bottom)/(SZ[3]-2)
function genZ!(Z::Vector{Float64}, SZ, b, dz::Float64)
    if mode==1 
        for k in 1:SZ[3]+1
            Z[k] = b[3] + (k-2)*dz
        end
    else
        # mode = 2, NZ[3]=15
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
    println("FillPowerGrid")
    FillPowerGrid!(ID, ox, Δh, SZ, Z)
    println("FillTSV")
    FillTSV!(ID, ox, Δh, SZ, Z)
    println("FillPlate")
    FillPlate!(ID, ox, Δh, SZ, Z)
    println("FillSolder")
    FillSolder!(ID, ox, Δh, SZ, Z)
    println("FillResin")
    FillResin!(ID, SZ)
end

# @param m_mode (1--Uniform in Z-dir, 2--Non-uniform in Z), while Uniform for X&Y
# @param NXY  Number of inner cells for X&Y dir.
# @param NZ   Number of inner cells for Z dir.
function model_test(m_mode::Int64, NXY::Int64, NZ::Int64=13)
    global mode = m_mode
    MX = MY = NXY + 2  # Number of CVs including boundary cells
    if mode==2
        NZ = 13
    end
    MZ = NZ + 2 # mode=1のときNZはセル数、mode=2のときNZはセル界面数

    dx::Float64 = 1.2 / NXY
    dy::Float64 = 1.2 / NXY
    dz::Float64 = zm5 / NZ
    Δh = (dx, dy, dz) 
    ox = (0.0, 0.0, 0.0)
    SZ = (MX, MY, MZ)
    println(SZ)
    println(Δh)
    
    ID = zeros(Int64, SZ[1], SZ[2], SZ[3])
    if mode==1
        Z = zeros(Float64, SZ[3]+1)
    else
        Z = zeros(Float64, SZ[3])
    end
    genZ!(Z, SZ, ox, Δh[3])

    @time fillID!(ID, ox, Δh, SZ, Z)


    #plot_slice(ID, SZ, "id.png")

    id_xy(ID, 0.2, SZ, Z, "id_z=0.2.png")
    
    # mode別の可視化関数切り替え
    if mode == 1
        # Cartesian格子：従来の関数
        id_yz(ID, 0.3, ox, Δh, SZ, Z, "id_x=0.3_$(NXY)x$(NZ).png")
    else
        # NonUniform格子：物理座標系対応関数
        id_yz_nu(ID, 0.3, ox, Δh, SZ, Z, "id_x=0.3_$(NXY)x$(NZ).png", NXY)
        plot_slice3(ID, 0.5, SZ, ox, Δh, Z, "slice.png", "NonUniform")
    end
 
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

function id_yz(d::Array{Int64,3}, x, ox, Δh, SZ, Z, fname)
    i = find_i(x, ox[1], Δh[1], SZ[1])
    s = d[i, 1:SZ[2], 1:SZ[3]]
    y_coords = [ox[2] + Δh[2] * (j - 1.5) for j in 1:SZ[2]]
    z_coords = [Z[k] for k in 1:SZ[3]]
    p = heatmap( z_coords, y_coords, s, 
        clims=(1,length(Fcolor)), 
        title="ID x=$x (x_index=$i)",
        c = Fcolor,
        colorbar=false,
        size = (300,600),
        aspect_ratio=:equal )
    savefig(p, fname)
end

function id_yz2(d::Array{Int64,3}, x, ox, Δh, SZ, Z, fname)
    i = find_i(x, ox[1], Δh[1], SZ[1])
    s = vec(d[i, 1:SZ[2], 1:SZ[3]])
    y_coords = repeat([ox[2] + Δh[2] * (j - 1.5) for j in 1:SZ[2]], SZ[3])
    z_coords = repeat([Z[k] for k in 1:SZ[3]], inner=SZ[2])
    p = scatter(z_coords, y_coords, marker = (:square, 2), 
         color = s, cgrad = Fcolor, colorbar_title = "Z value",
         size = (300,600) )
    savefig(p, fname)
end

# NonUniform格子対応のid_yz関数（物理座標系を正しく反映?）
function id_yz_nu(d::Array{Int64,3}, x, ox, Δh, SZ, Z, fname, NXY)
    i = find_i(x, ox[1], Δh[1], SZ[1])
    s = d[i, 1:SZ[2], 1:SZ[3]]
    
    # Y方向座標軸（物理座標）
    y_coords = [ox[2] + Δh[2] * (j - 1.5) for j in 1:SZ[2]]
    
    # Z方向座標軸（NonUniform物理座標）
    z_coords = [Z[k] for k in 1:SZ[3]]
    
    p = heatmap( z_coords, y_coords, s, 
        clims=(1,length(Fcolor)), 
        title="ID x=$x (NonUniform, $(NXY)x$(length(Z)-2))",
        xlabel="Z-coordinate [physical]",
        ylabel="Y-coordinate [physical]",
        c = Fcolor,
        colorbar=false,
        size = (300,600),
        aspect_ratio=:equal )
    savefig(p, fname)
end



function plot_slice(d::Array{Int64,3}, SZ, fname)
    j = 16 #div(SZ[3],2)
    s = d[1:SZ[1],j,1:SZ[3]]

    #heatmap(s, xlabel="X-axis", ylabel="Z-axis", title="Y=$j")
    p = contour(s, fill=true, c=:thermal, xlabel="Z-axis", ylabel="X-axis", title="Y=$j", size=(600, 600))
    savefig(p, fname)
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
function plot_slice3(d::Array{Int64,3}, x, SZ, ox, Δh, Z::Vector{Float64}, fname, label::String="")
    i = find_i(x, ox[1], Δh[1], SZ[1])
    s = convert(Array{Float64}, d[i, 1:SZ[2], 1:SZ[3]])
    
    y_coords = [ox[2] + Δh[2] * (j - 1.5) for j in 1:SZ[2]]
    z_coords = [Z[k] for k in 1:SZ[3]]
    
    # 物理座標軸でプロット
    p = contour(z_coords, y_coords, s, 
                fill=true, 
                c=:thermal, 
                xlabel="Z", 
                ylabel="Y", 
                title="YZ-sect. (X=$x, i=$i, $label) $(SZ[2])x$(length(Z)-2)", 
                size=(300, 600),
                aspect_ratio=:equal)
    savefig(p, fname)
    
    return p
end


model_test(1,240,120)
model_test(2,240)
#model_test(1,480,240)
