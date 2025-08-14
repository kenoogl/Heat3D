using Printf
using Random
using LinearAlgebra

include("heat3d_Cartesian.jl")
include("heat3d_CartesianII.jl")
include("heat3d_NonUniform.jl")
include("../model/modelA.jl")
include("plotter.jl")
include("const.jl")


#=
@brief 境界条件 熱源項の設定
@param [in,out] b    右辺項
@param [in]     SZ   配列長
@param [in]     ID   識別子配列
=#
function calRHS!(b::Array{Float64,3}, ID::Array{UInt8,3}, SZ)
    for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
        if ID[i,j,k] == modelA.pwrsrc["id"]
            b[i,j,k] = -Constant.Q_src
        end
    end
end

#=
@brief マスク指定
@param [in]     m      マスク
@param [in]     SZ     配列長
@note  非計算セルは0.0、計算セルには1.0
        ディリクレ型　温度指定 : mask_BC=0.0, θ=θ_BC, λ_BC=λ_inner
        ノイマン型　断熱境界 : λ_BC=0.0
=#
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

function setMask3!(m::Array{Float64,3}, SZ)

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

function conditions(F, SZ, Δh, solver, smoother)
    if mode==1
        @printf(F, "Problem : CUBE on Cartesian grid\n")
    elseif mode==2
        @printf(F, "Problem : CUBE on NonUniform grid from file\n")
    elseif mode==3
        @printf(F, "Problem : IC on NonUniform grid (Opt. 13 layers)\n")
    elseif mode==4
        @printf(F, "Problem : IC on Cartesian grid\n")
    end

    @printf(F, "Grid  : %d %d %d\n", SZ[1], SZ[2], SZ[3])
    @printf(F, "Pitch : %6.4e %6.4e %6.4e\n", Δh[1], Δh[2], Δh[3])
    if solver=="pbicgstab"
        @printf(F, "Solver: %s with smoother %s\n", solver, smoother)
    else
        @printf(F, "Solver: %s\n", solver)
    end
    @printf(F, "ItrMax : %e\n", Constant.ItrMax)
    @printf(F, "ε      : %e\n", Constant.tol)
    @printf(F, "θ_amb  : %e\n", Constant.θ_amb)
    @printf(F, "θ_pcb  : %e\n", Constant.θ_pcb)
    @printf(F, "HT_top : %e\n", Constant.HT_top)
    @printf(F, "HT_side: %e\n", Constant.HT_side)
    @printf(F, "Q_src  : %e\n", Constant.Q_src)
end

#=
@brief 熱伝導率の設定
@param [in]     λ      熱伝導率
@param [in]     SZ     配列長
=#
function setMatOuter!(λ::Array{Float64,3}, SZ)
    for j in 1:SZ[2], i in 1:SZ[1]
        λ[i,j,1    ] = λ[i,j,2    ]
        λ[i,j,SZ[3]] = 0.0
    end

    for k in 1:SZ[3], i in 1:SZ[1]
        λ[i,1    ,k] = 0.0 #λ[i,2    ,k]
        λ[i,SZ[2],k] = 0.0 #λ[i,SZ[2]-1,k]
    end

    for k in 1:SZ[3], j in 1:SZ[2]
        λ[1    ,j,k] = 0.0 #λ[2    ,j,k]
        λ[SZ[1],j,k] = 0.0 #λ[SZ[1]-1,j,k]
    end
end



#= 
@param [in] SZ       内部セル数
@param [in] Δh       セル幅
@param [in] θ        解ベクトル
@param [in] solver   ["jacobi", "sor", "pbicgstab"]
@param [in] smoother ["jacobi", "gs", ""]
=#
function main(SZ, ox, Δh, θ, Z, ΔZ, solver, smoother)

    z_st::Int64=0
    z_ed::Int64=0

    if mode==1 || mode==4
        z_st = 2
        z_ed = SZ[3]-1
    elseif mode==2
        z_st = 3
        z_ed = SZ[3]-2
    # mode=3の設定では上面は熱伝達、下面は定温
    elseif mode==3 
        z_st = 3
        z_ed = SZ[3]-1
    end

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0 # default
    #=
    if example==1
        setMat_1!(λ, SZ)
    elseif example==2
        setMat_2!(λ, SZ)
    elseif example==3
        setMat_3!(λ, SZ)
    else
        println("example error")
        return
    end
    plot_slice2(λ, SZ, "lambda.png")
    =#

    ID = zeros(UInt8, SZ[1], SZ[2], SZ[3]) # mode=3のときのみ有効

    if mode==3
        modelA.fillID!(mode, ID, ox, Δh, SZ, Z)
        modelA.setLambda!(λ, SZ, ID)
        setMatOuter!(λ, SZ)
        plot_slice_xz_nu(1, λ, 0.3e-3, SZ, ox, Δh, Z, "alpha3.png", "α")
    elseif mode==4
        modelA.fillID!(mode, ID, ox, Δh, SZ, Z)
        modelA.setLambda!(λ, SZ, ID)
        setMatOuter!(λ, SZ)
        plot_slice_xz(1, λ, Z, 0.3e-3, SZ, ox, Δh, "alpha4.png", "α")
    end


    if mode==1
        Cartesian.boundary_condition!(θ, SZ, ox, Δh)
    elseif mode==2
        NonUniform.boundary_condition!(θ, SZ, ox, Δh)
    elseif mode==3
        NonUniform.boundary_condition3!(θ, SZ)
    else
        CartesianII.boundary_condition4!(θ, SZ)
    end

    b = zeros(Float64, SZ[1], SZ[2], SZ[3])
    if mode==3 || mode==4
        calRHS!(b, ID, SZ)
    end

    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    if mode==1 || mode==2 || mode==4
        setMask!(mask, SZ)
    else
        setMask3!(mask, SZ)
    end
    #plot_slice2(mask, SZ, "mask.png")


    wk = zeros(Float64, SZ[1], SZ[2], SZ[3])

    if solver=="pbicgstab"
        pcg_p  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_p_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r0 = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_q  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_t_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end

    F = open("log.txt", "w")
    conditions(F, SZ, Δh, solver, smoother)

    if solver=="sor"
        if mode==1
            Cartesian.solveSOR!(θ, SZ, λ, b, mask, Δh, ω, F)
        elseif mode==4
            CartesianII.solveSOR!(θ, SZ, λ, b, mask, Δh, ω, F)
        else
            NonUniform.solveSOR!(θ, SZ, λ, b, mask, Δh, ω, Z, ΔZ, z_st, z_ed, F)
        end
    elseif solver=="jacobi"
        if mode==1
            Cartesian.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, F)
        elseif mode==4
            CartesianII.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, F)
        else
            NonUniform.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, ω, Z, ΔZ, z_st, z_ed, F)
        end
    elseif solver=="pbicgstab"
        if mode==1
            Cartesian.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, ox, Δh, SZ, smoother, F, mode)
        elseif mode==4
            CartesianII.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, ox, Δh, SZ, smoother, F, mode)
        else
            NonUniform.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, 
                ox, Δh, SZ, Z, ΔZ, z_st, z_ed, smoother, F, mode)
        end
    else
        println("solver error")
    end

    min_val = minimum(θ)
    max_val = maximum(θ)
    println("θmin=", min_val, " θmax=", max_val, " L2 norm of θ=",sqrt(norm(θ,2)))

    close(F)
end

# Z軸座標の生成
# mode==1のとき等間隔なのでZ座標不要、Δz=(top-bottom)/(SZ[3]-2)
function genZ!(Z::Vector{Float64}, ΔZ::Vector{Float64}, SZ, ox, dz)
    mz=SZ[3]
    if mode==1 || mode==4
        for k in 1:mz+1
            Z[k] = ox[3] + (k-2)*dz
        end
    elseif mode==2
        read_coord, numNodes = NonUniform.read_grid_file()
        if numNodes != mz-2
            println("Number of genereted grid is not match with parameter NZ")
            exit(0)
        end
        for k in 1:mz-2
            Z[k+1] = read_coord[k]
        end
        Z[1] = 2*Z[2] - Z[3]
        Z[mz] = 2*Z[mz-1] - Z[mz-2]
    else
        # mode = 2, NZ[3]=15
        Z[1] = -0.05e-3
        Z[2] = 0.0 #zm0
        Z[3] = 0.05e-3
        Z[4] = 0.1e-3 #zm1
        Z[5] = 0.2e-3-modelA.pg_dpth
        Z[6] = 0.2e-3
        Z[7] = 0.25e-3 #zm2
        Z[8] = 0.35e-3-modelA.pg_dpth
        Z[9] = 0.35e-3
        Z[10]= 0.4e-3 #zm3
        Z[11]= 0.5e-3-modelA.pg_dpth
        Z[12]= 0.5e-3
        Z[13]= 0.55e-3 #zm4
        Z[14]= 0.6e-3  #zm5
        Z[15]= 0.65e-3
    end

    for k in 2:mz-1
        ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])
    end
    ΔZ[2] = 0.5*ΔZ[2]
    ΔZ[mz-1] = 0.5*ΔZ[mz-1]
    #println(Z)
    #println(ΔZ)
end

#=
# Visio用のsphフォーマット出力
function writeSPH(size, org, pch, step, time, var)
    ccall((:write_sph_d, "./iosph.so"), Nothing,
    (Ref{Int},    # size
     Ref{Float64},  # org
     Ref{Float64},  # pch
     Ref{Int},    # step
     Ref{Float64},  # time
     Ref{Float64}), # var
     size, org, pch, step, time, var)
end
=#

#=
@param [in] mode (
                    1--Uniform in Z-dir, (Cartesian) : CUBE
                    2--Uniform for X&Y and Non-uniform in Z, but data is uniform : CUBE
                    3--Uniform for X&Y and Non-uniform in Z : modelA
                    4--Uniform : modelA
@param NXY  Number of inner cells for X&Y dir.
@param NZ   Number of inner cells for Z dir.
@param [in] solver    ["jacobi", "sor", "pbicgstab"]
@param [in] smoother  ["jacobi", "gs", ""]
=#
function q3d(m_mode::Int, NXY::Int, NZ::Int, solver::String="sor", smoother::String="")
    global mode = m_mode

    MX = MY = NXY + 2  # Number of CVs including boundaries
    if mode==3
        NZ = 13
    end
    MZ = NZ + 2

    if !(mode==1 ||mode==2 || mode==3 || mode==4)
        println("mode error")
        return
    end

    if mode==1
        if NXY != NZ
            println("NXY must be equal to NZ")
            return
        end
    end

    dh::Float64 = 1.0

    if mode==1 || mode==2
        dh = 1.0 / NXY
    else
        dh = 5e-6 # 5μ [m]
    end
    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0) #原点を仮定

    println(SZ)
    println(Δh)

    if mode==1 || mode==4
        Z = zeros(Float64, SZ[3]+1)
    else
        Z = zeros(Float64, SZ[3])
    end

    ΔZ= zeros(Float64, SZ[3]-1)
    #@show typeof(Z)

    genZ!(Z, ΔZ, SZ, ox, Δh[3])

    if mode<=2
        exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end 

    if mode==1
        Cartesian.exact_solution!(exact, SZ, ox, Δh)
        plot_slice_xz(1, exact, Z, 0.5e-3, SZ, ox, Δh, "exact.png", "Exact")
    elseif mode==2
        NonUniform.exact_solution!(exact, SZ, ox, Δh, Z)
        plot_slice_xz_nu(1, exact, 0.5e-3, SZ, ox, Δh, Z, "exact_nu.png", "Exact")
    end
    
    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    θ .= Constant.θ_amb # 初期温度設定

    @time main(SZ, ox, Δh, θ, Z, ΔZ, solver, smoother)

    
    if mode==1 || mode==2
        nrm = sqrt(norm(θ-exact,2))
        F = open("log.txt", "a")
        @printf(F, "L2 norm of error = %24.14E\n", nrm)
        close(F)
    end
    

    if mode==1
        plot_slice_xz(2, θ, Z, 0.5e-3, SZ, ox, Δh, "p.png", "solution")
        plot_slice_xz(2, θ-exact, Z, 0.5e-3, SZ, ox, Δh, "diff.png", "diff")
    elseif mode==2
        plot_slice_xz_nu(2, θ, 0.5e-3, SZ, ox, Δh, Z, "p_nu.png", "solution")
        plot_slice_xz_nu(2, θ-exact, 0.5e-3, SZ, ox, Δh, Z, "diff_nu.png", "diff")
    elseif mode==3
        plot_slice_xz_nu(1, θ, 0.5e-3, SZ, ox, Δh, Z, "temp3_xz_nu.png", "temperature[K]")
        plot_slice_xy_nu(1, θ, 0.195e-3, SZ, ox, Δh, Z, "temp3_xy_nu.png", "temperature[K]")
    else
        plot_slice_xz(2, θ, Z, 0.3e-3, SZ, ox, Δh, "temp4_xz.png")
        plot_slice_xy(2, θ, 0.48e-3, SZ, ox, Δh, Z, "temp4_xy.png")
        plot_line_z(θ, SZ, ox, Δh, "tempZ.png")
    end

end

#q3d(3, 240, 13, "pbicgstab") # ここで本実行し、計測
#q3d(1, 25, 25, "pbicgstab")
q3d(4, 240, 120, "pbicgstab", "gs") 