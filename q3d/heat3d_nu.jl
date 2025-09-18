using Printf
using LinearAlgebra
using FLoops
using ThreadsX

include("common.jl")
include("heat3d_Cartesian.jl")
include("heat3d_NonUniform.jl")
include("../model/modelA.jl")
include("plotter.jl")
include("convergence_history.jl")
include("parse_log_residuals.jl")
include("boundary_conditions.jl")
include("Zcoord.jl")


"""
Mode1用の境界条件（Cartesian格子の立方体問題）
Z面に三角関数分布の等温条件、X,Y面は等温0K
"""
function set_mode1_bc_parameters()
    # 各面の境界条件を定義
    x_minus_bc = BoundaryConditions.isothermal_bc(0.0)           # X軸負方向面: 等温0K
    x_plus_bc  = BoundaryConditions.isothermal_bc(0.0)           # X軸正方向面: 等温0K
    y_minus_bc = BoundaryConditions.isothermal_bc(0.0)           # Y軸負方向面: 等温0K
    y_plus_bc  = BoundaryConditions.isothermal_bc(0.0)           # Y軸正方向面: 等温0K
    z_minus_bc = BoundaryConditions.isothermal_bc(0.0)           # Z軸負方向面: 等温0K (後で三角関数分布で上書き)
    z_plus_bc  = BoundaryConditions.isothermal_bc(0.0)           # Z軸正方向面: 等温0K (後で三角関数分布で上書き)
    
    # 境界条件セットを作成
    return BoundaryConditions.create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
Mode2用の境界条件（Non-uniform格子の立方体問題）
Z面に三角関数分布の等温条件、X,Y面は等温0K
"""
function set_mode2_bc_parameters()
    # 各面の境界条件を定義
    x_minus_bc = BoundaryConditions.isothermal_bc(0.0)           # X軸負方向面: 等温0K
    x_plus_bc  = BoundaryConditions.isothermal_bc(0.0)           # X軸正方向面: 等温0K
    y_minus_bc = BoundaryConditions.isothermal_bc(0.0)           # Y軸負方向面: 等温0K
    y_plus_bc  = BoundaryConditions.isothermal_bc(0.0)           # Y軸正方向面: 等温0K
    z_minus_bc = BoundaryConditions.isothermal_bc(0.0)           # Z軸負方向面: 等温0K (後で三角関数分布で上書き)
    z_plus_bc  = BoundaryConditions.isothermal_bc(0.0)           # Z軸正方向面: 等温0K (後で三角関数分布で上書き)
    
    # 境界条件セットを作成
    return BoundaryConditions.create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
Mode3用の境界条件（NonUniform格子のIC問題）  
Z下面: PCB温度、Z上面: 熱伝達、側面: 断熱
"""
function set_mode3_bc_parameters()
    θ_amb = 300.0 # [K]
    θ_pcb = 300.0 # [K]
    HT_top = 2.98e-4 # 5 [W/(m^2 K)] / (\rho C)_silicon > [m/s]
    HT_side = 2.98e-6 # 5 [W/(m^2 K)] / (\rho C)_silicon > [m/s]
    
    # 各面の境界条件を定義
    x_minus_bc = BoundaryConditions.convection_bc(HT_side, θ_amb)
    x_plus_bc  = BoundaryConditions.convection_bc(HT_side, θ_amb)
    y_minus_bc = BoundaryConditions.convection_bc(HT_side, θ_amb)  
    y_plus_bc  = BoundaryConditions.convection_bc(HT_side, θ_amb)
    z_minus_bc = BoundaryConditions.isothermal_bc(θ_pcb)                  # Z軸負方向面: PCB温度
    z_plus_bc  = BoundaryConditions.convection_bc(HT_top, θ_amb)          # Z軸正方向面: 熱伝達
    
    # 境界条件セットを作成
    return BoundaryConditions.create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
Mode4用の境界条件（IC package問題）
Z下面: PCB温度、Z上面: 周囲温度、側面: 周囲温度
"""
function set_mode4_bc_parameters()
    θ_amb = 300.0 # [K]
    θ_pcb = 300.0 # [K]
    HT_top = 2.98e-4 # 5 [W/(m^2 K)] / (\rho C)_silicon > [m/s]
    HT_side = 2.98e-6 # 5 [W/(m^2 K)] / (\rho C)_silicon > [m/s]
    # 各面の境界条件を定義
    x_minus_bc = BoundaryConditions.convection_bc(HT_side, θ_amb)
    x_plus_bc  = BoundaryConditions.convection_bc(HT_side, θ_amb)
    y_minus_bc = BoundaryConditions.convection_bc(HT_side, θ_amb)  
    y_plus_bc  = BoundaryConditions.convection_bc(HT_side, θ_amb)
    z_minus_bc = BoundaryConditions.isothermal_bc(θ_pcb)         # Z軸負方向面: PCB温度
    z_plus_bc  = BoundaryConditions.convection_bc(HT_top, θ_amb)
    
    # 境界条件セットを作成
    return BoundaryConditions.create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end


"""
@brief 計算領域内部の熱源項の設定、ガイドセル部分は境界条件定数で利用
@param [in,out] b    右辺項
@param [in]     ID   識別子配列
"""
function HeatSrc!(b::Array{Float64,3}, ID::Array{UInt8,3}, par)
    backend = get_backend(par)
    SZ = size(b)
    
    @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        if ID[i,j,k] == modelA.pwrsrc["id"]
            b[i,j,k] = -Constant.Q_src
        end
    end
end

"""
@brief CUBE問題の境界条件 Cartesian
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     Δh   セル幅
"""
function bc_cube!(p::Array{Float64,3}, ox, Δh, par)
    backend = get_backend(par)
    SZ = size(p)

    @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        a = sin(π*x)*sin(π*y)
        p[i,j,1    ] = a
        p[i,j,SZ[3]] = a
    end
end

"""
@brief CUBE問題の境界条件 NonUniform
@param [in,out] p    解ベクトル
@param [in]     SZ   配列長
@param [in]     Δh   セル幅
"""
function bc_cube_nu!(p::Array{Float64,3}, ox, Δh, par)
    backend = get_backend(par)
    SZ = size(p)
    
    @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        x = ox[1] + Δh[1]*(i-1.5)
        y = ox[2] + Δh[2]*(j-1.5)
        a = sin(π*x)*sin(π*y)
        p[i,j,2    ] = a
        p[i,j,SZ[3]-1] = a
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
    @printf(F, "ε      : %e\n", itr_tol)
end

"""
@param z_range Z方向の計算内部領域インデクス
"""
function preprocess!(λ, Z, ΔZ, ox, Δh, ID)
    SZ = size(λ)
    Zcoordinate.genZ!(Z, ΔZ, SZ, ox, Δh[3], mode)

    if mode==3 || mode==4
        modelA.fillID!(mode, ID, ox, Δh, Z)
        modelA.setLambda!(λ, ID)
    end
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

function dif(θ, exact)
    SZ = size(θ)
    d=0.0
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
       s = (θ[i,j,k]-exact[i,j,k])
       d = d + s*s
    end
    return sqrt(d)
end

"""
@param [in] Δh       セル幅
@param [in] θ        解ベクトル
@param [in] solver   ["sor", "pbicgstab", "cg"]
@param [in] smoother ["gs", ""]
"""
function main(ox, Δh, θ, b, mask, Z, ΔZ, ID, λ, solver, smoother, z_range, bc_set, par)
    # 収束履歴の初期化
    conv_data = ConvergenceData(solver, smoother)

    SZ = size(θ)

    HF = zeros(Float64, 6)
    HT = zeros(Float64, 6)

    HF[1] = bc_set.x_minus.heat_flux
    HF[2] = bc_set.x_plus.heat_flux
    HF[3] = bc_set.y_minus.heat_flux
    HF[4] = bc_set.y_plus.heat_flux
    HF[5] = bc_set.z_minus.heat_flux
    HF[6] = bc_set.z_plus.heat_flux
    #println(HF)

    HT[1] = bc_set.x_minus.heat_transfer_coefficient
    HT[2] = bc_set.x_plus.heat_transfer_coefficient
    HT[3] = bc_set.y_minus.heat_transfer_coefficient
    HT[4] = bc_set.y_plus.heat_transfer_coefficient
    HT[5] = bc_set.z_minus.heat_transfer_coefficient
    HT[6] = bc_set.z_plus.heat_transfer_coefficient
    #println(HT)

    if mode==3 || mode==4
        HeatSrc!(b, ID, par)
    end

    if solver=="pbicgstab"
        pcg_p  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_p_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_r0 = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_q  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_s_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
        pcg_t_ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    elseif solver=="cg"
        cg_p  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        cg_r  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        cg_ax = zeros(Float64, SZ[1], SZ[2], SZ[3])
        cg_ap = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end

    F = open("log.txt", "w")
    conditions(F, SZ, Δh, solver, smoother)

    if solver=="sor"
        if mode==1 || mode==4
            Cartesian.solveSOR!(θ, λ, b, mask, Δh, Constant.ω, F, itr_tol, HF, HT, par)
        elseif mode==2 || mode==3
            NonUniform.solveSOR!(θ, λ, b, mask, Δh, Constant.ω, Z, ΔZ, z_range, F, itr_tol, HF, HT, par)
        end
    elseif solver=="pbicgstab"
        if mode==1 || mode==4
            Cartesian.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, ox, Δh, smoother, F, mode, itr_tol, HF, HT, par)
        elseif mode==2 || mode==3
            NonUniform.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, 
                ox, Δh, Z, ΔZ, z_range, smoother, F, mode, itr_tol, HF, HT, par)
        end
    elseif solver=="cg"
        if mode==1 || mode==4
            Cartesian.CG!(θ, b, cg_p, cg_r, cg_ax, cg_ap, λ, mask, Δh, F, itr_tol, HF, HT, par)
        end
    else
        println("solver error")
    end

    s = θ[2:SZ[1]-1, 2:SZ[2]-1, z_range[1]:z_range[2]]
    min_val = minimum(s)
    max_val = maximum(s)
    @printf(F, "θmin=%e  θmax=%e  L2 norm of θ=%e\n", min_val, max_val, norm(s,2))

    close(F)
    
    # ログファイルから残差データを解析してconv_dataに追加
    if solver == "pbicgstab" || solver == "sor" || solver == "cg" # || solver == "jacobi"
        parse_residuals_from_log!(conv_data, "log.txt")
    end
    
    # 収束履歴データを返す
    return conv_data
end

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
function q3d(m_mode::Int, NXY::Int, NZ::Int, 
         solver::String="sor", smoother::String=""; 
         epsilon::Float64=1.0e-6, par::String="thread")
    global mode = m_mode
    global itr_tol = epsilon

    println("Julia version: $(VERSION)")
    
    if par=="sequential"
        println("Sequential execution")
    elseif par=="thread"
        println("Available num. of threads: ", Threads.nthreads())
    else
        println("Invalid paralle mode")
        exit()
    end

    println("="^60)

    MX = MY = NXY + 2  # Number of CVs including boundaries
    MZ = NZ + 2

    if !(mode==1 ||mode==2 || mode==3 || mode==4)
        println("mode error")
        return
    end

    if mode==1 || mode==2
        if NXY != NZ
            println("NXY must be equal to NZ")
            return
        end
    end

    dh::Float64 = 1.0

    if mode==1 || mode==2
        dh = 1.0 / NXY
    elseif mode==3 || mode==4
        dh = 1.2e-3 / NXY
    end
    dh = round(dh,digits=8) #4.9999.... >> 5.0にしたい

    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0) #原点を仮定
    z_range = zeros(Int64, 2)

    println(SZ, "  Itr.ε= ", itr_tol)
    #println(Δh)

    if mode==1 || mode==4
        Z = zeros(Float64, SZ[3]+1)
    else
        Z = zeros(Float64, SZ[3])
    end

    ΔZ= zeros(Float64, SZ[3]-1)
    #@show typeof(Z)

    λ = Array{Float64}(undef, SZ[1], SZ[2], SZ[3])
    λ .= 1.0 # default

    ID   = zeros(UInt8, SZ[1], SZ[2], SZ[3]) # mode=3のときのみ有効
    θ    = zeros(Float64, SZ[1], SZ[2], SZ[3])
    b    = zeros(Float64, SZ[1], SZ[2], SZ[3])
    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    
    @time preprocess!(λ, Z, ΔZ, ox, Δh, ID)
    #plot_slice2(mask, SZ, "mask.png")

    if mode==3
        plot_slice_xz_nu(1, mode, λ, 0.3e-3, SZ, ox, Δh, Z, "alpha3.png", "α")
    elseif mode==4
        plot_slice_xz(1, mode, λ, Z, 0.3e-3, SZ, ox, Δh, "alpha4.png", "α")
    end

    if mode<=2
        exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end 

    if mode==1
        Cartesian.exact_solution!(exact, ox, Δh, par)
        plot_slice_xz(1, mode, exact, Z, 0.5, SZ, ox, Δh, "exact.png", "Exact")
    elseif mode==2
        NonUniform.exact_solution!(exact, ox, Δh, Z, par)
        plot_slice_xz_nu(1, mode, exact, 0.5, SZ, ox, Δh, Z, "exact_nu.png", "Exact")
    end
    
    # Boundary condition
    if mode==1 # Cartesian, CUBE
        bc_set = set_mode1_bc_parameters()
        z_range[1] = 2
        z_range[2] = SZ[3]-1
        θ_init = 0.0
    elseif mode==2 # Non-uniform, CUBE
        bc_set = set_mode2_bc_parameters()
        z_range[1] = 3
        z_range[2] = SZ[3]-2
        θ_init = 0.0
    elseif mode==3 # Non-uniform, 3D-IC
        bc_set = set_mode3_bc_parameters()
        z_range[1] = 3
        z_range[2] = SZ[3]-1
        θ_init = 300.0
    elseif mode==4 # Cartesian, 3D-IC
        bc_set = set_mode4_bc_parameters()
        z_range[1] = 2
        z_range[2] = SZ[3]-1
        θ_init = 300.0
    end

    θ .= θ_init # 初期温度設定

    BoundaryConditions.print_boundary_conditions(bc_set)
    BoundaryConditions.apply_boundary_conditions!(θ, λ, mask, bc_set, mode)
    if mode==1
        bc_cube!(θ, ox, Δh, par) # Z方向の上下面の分布を上書き
    elseif mode==2
        bc_cube_nu!(θ, ox, Δh, par)
    end

    tm = @elapsed conv_data = main(ox, Δh, θ, b, mask, Z, ΔZ, ID, λ, solver, smoother, z_range, bc_set, par)

    
    if mode==1 || mode==2
        nrm = dif(θ, exact)
        F = open("log.txt", "a")
        @printf(F, "L2 norm of serror = %24.14E\n", nrm)
        close(F)
    end
    

    if mode==1
        plot_slice_xz(2, mode, θ, Z, 0.5, SZ, ox, Δh, "p.png", "solution")
        plot_slice_xz(2, mode, abs.(θ-exact), Z, 0.5, SZ, ox, Δh, "diff.png", "diff")
        println(" L2 norm of θ-exact=",dif(θ, exact))
    elseif mode==2
        plot_slice_xz_nu(2, mode, θ, 0.5, SZ, ox, Δh, Z, "p_nu.png", "solution")
        plot_slice_xz_nu(2, mode, abs.(θ-exact), 0.5, SZ, ox, Δh, Z, "diff_nu.png", "diff")
        println(" L2 norm of θ-exact=",dif(θ, exact))
    elseif mode==3
        plot_slice_xz_nu(2, mode, θ, 0.3e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.3.png")
        plot_slice_xz_nu(2, mode, θ, 0.4e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.4.png")
        plot_slice_xz_nu(2, mode, θ, 0.5e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.5.png")
        plot_slice_xy_nu(2, mode, θ, 0.18e-3, SZ, ox, Δh, Z, "temp3_xy_nu_z=0.18.png")
        plot_slice_xy_nu(2, mode,θ, 0.33e-3, SZ, ox, Δh, Z, "temp3_xy_nu_z=0.33.png")
        plot_slice_xy_nu(2, mode, θ, 0.48e-3, SZ, ox, Δh, Z, "temp3_xy_nu_z=0.48.png")
        plot_line_z_nu(θ, SZ, ox, Δh, Z, 0.6e-3, 0.6e-3,"temp3Z_ctr", "Center")
        plot_line_z_nu(θ, SZ, ox, Δh, Z, 0.4e-3, 0.4e-3,"temp3Z_tsv", "TSV")
    else
        plot_slice_xz(2, mode, θ, Z, 0.3e-3, SZ, ox, Δh, "temp4_xz_y=0.3.png")
        plot_slice_xz(2, mode, θ, Z, 0.4e-3, SZ, ox, Δh, "temp4_xz_y=0.4.png")
        plot_slice_xz(2, mode, θ, Z, 0.5e-3, SZ, ox, Δh, "temp4_xz_y=0.5.png")
        plot_slice_xy(2, mode, θ, 0.18e-3, SZ, ox, Δh, Z, "temp4_xy_z=0.18.png")
        plot_slice_xy(2, mode, θ, 0.33e-3, SZ, ox, Δh, Z, "temp4_xy_z=0.33.png")
        plot_slice_xy(2, mode, θ, 0.48e-3, SZ, ox, Δh, Z, "temp4_xy_z=0.48.png")
        plot_line_z(θ, SZ, ox, Δh, 0.6e-3, 0.6e-3, "temp4Z_ctr", "Center")
        plot_line_z(θ, SZ, ox, Δh, 0.4e-3, 0.4e-3, "temp4Z_tsv", "TSV")
    end
    
    # 収束履歴の出力（反復解法の場合のみ）
    if solver == "pbicgstab" || solver == "sor" || solver == "cg" # || solver == "jacobi"
        # 収束グラフとCSV出力
        conv_filename = "convergence_$(solver)_mode$(mode)_$(NXY)x$(NZ)"
        if !isempty(smoother)
            conv_filename *= "_$(smoother)"
        end
        
        # プロットとCSV出力
        try
            plot_convergence_curve(conv_data, "$(conv_filename).png", target_tol=itr_tol, show_markers=false)
            export_convergence_csv(conv_data, "$(conv_filename).csv")
        catch e
            println("Error in convergence history output: $e")
        end
        
        # 収束情報の表示
        info = get_convergence_info(conv_data)
        if !isempty(info)
            println("\n=== Convergence Information ===")
            println("Mode: $mode, Grid: $(NXY)x$(NZ)")
            println("Solver: $(info["solver"]), Smoother: $(info["smoother"])")
            println("Iterations: $(info["iterations"])")
            initial_res_str = @sprintf("%.6E", info["initial_residual"])
            final_res_str = @sprintf("%.6E", info["final_residual"])
            conv_rate_str = @sprintf("%.6E", info["convergence_rate"])
            reduction_str = @sprintf("%.2f", info["reduction_factor"])
            println("Initial residual: $initial_res_str")
            println("Final residual: $final_res_str")
            println("Residual reduction factor: $conv_rate_str")
            println("Order reduction: $reduction_str")
            println("===============================")
        end
    end

    println(tm, "[sec]")
    println(" ")
end

if abspath(PROGRAM_FILE) == @__FILE__
  #q3d(1, 25, 25, "sor", epsilon=1.0e-4, par="thread")
  #q3d(1, 25, 25, "pbicgstab", epsilon=1.0e-4, par="sequential")
  #q3d(1, 25, 25, "pbicgstab", "gs", epsilon=1.0e-4, par="thread")
  #q3d(1, 25, 25, "cg", epsilon=1.0e-4, par="thread")

  #q3d(2, 25, 25, "sor", epsilon=1.0e-4, par="sequential")
  #q3d(2, 25, 25, "pbicgstab", "gs", epsilon=1.0e-4, par="sequential")

  q3d(3, 240, 31, "pbicgstab", "gs", epsilon=1.0e-4, par="sequential")
  #q3d(3, 120, 31, "pbicgstab", "gs", epsilon=1.0e-4, par="sequential")
  
  #q3d(4, 240, 120, "sor", epsilon=1.0e-4, par="sequential")
  #q3d(4, 240, 120, "pbicgstab", "gs", epsilon=1.0e-4, par="thread") 
end
#q3d(1, 25, 25, "pbicgstab")
