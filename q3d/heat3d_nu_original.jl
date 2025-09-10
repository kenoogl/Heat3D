using Printf
using Random
using LinearAlgebra

include("heat3d_Cartesian.jl")
include("heat3d_CartesianII.jl")
include("heat3d_NonUniform.jl")
include("heat3d_NonUniformII.jl")
include("../model/modelA.jl")
include("plotter.jl")
include("const.jl")
include("convergence_history.jl")
include("parse_log_residuals.jl")


#=
@brief 境界条件 熱源項の設定
@param [in,out] b    右辺項
@param [in]     ID   識別子配列
=#
function calRHS!(b::Array{Float64,3}, ID::Array{UInt8,3})
    SZ = size(b)
    for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
        if ID[i,j,k] == modelA.pwrsrc["id"]
            b[i,j,k] = -Constant.Q_src
        end
    end
end

#=
@brief マスク指定
@param [in]     m      マスク
@note  非計算セルは0.0、計算セルには1.0
        ディリクレ型　温度指定 : mask_BC=0.0, θ=θ_BC, λ_BC=λ_inner
        ノイマン型　断熱境界 : λ_BC=0.0
=#
function setMask!(m::Array{Float64,3})
    SZ = size(m)

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
    @printf(F, "ε      : %e\n", itr_tol)
    @printf(F, "θ_amb  : %e\n", Constant.θ_amb)
    @printf(F, "θ_pcb  : %e\n", Constant.θ_pcb)
    @printf(F, "HT_top : %e\n", Constant.HT_top)
    @printf(F, "HT_side: %e\n", Constant.HT_side)
    @printf(F, "Q_src  : %e\n", Constant.Q_src)
end

#=
@brief 熱伝導率の設定
@param [in]     λ      熱伝導率
=#
function setMatOuter!(λ::Array{Float64,3})
    SZ = size(λ)
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


function preprocess!(λ, Z, ΔZ, ox, Δh, ID)
    SZ = size(λ)
    genZ!(Z, ΔZ, SZ, ox, Δh[3])

    if mode==3 || mode==4
        modelA.fillID!(mode, ID, ox, Δh, SZ, Z)
        modelA.setLambda!(λ, SZ, ID)
        setMatOuter!(λ)
    end
end

#= 
@param [in] Δh       セル幅
@param [in] θ        解ベクトル
@param [in] solver   ["jacobi", "sor", "pbicgstab"]
@param [in] smoother ["jacobi", "gs", ""]
=#
function main(ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother)
    # 収束履歴の初期化
    conv_data = ConvergenceData(solver, smoother)

    z_st::Int64=0
    z_ed::Int64=0
    SZ = size(θ)

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

    if mode==1
        Cartesian.boundary_condition!(θ, ox, Δh)
    elseif mode==2
        NonUniform.boundary_condition!(θ, SZ, ox, Δh)
    elseif mode==3
        NonUniformII.boundary_condition3!(θ, SZ)
    elseif mode==4
        CartesianII.boundary_condition4!(θ)
    end

    b = zeros(Float64, SZ[1], SZ[2], SZ[3])
    if mode==3 || mode==4
        calRHS!(b, ID)
    end

    mask = ones(Float64, SZ[1], SZ[2], SZ[3])
    setMask!(mask)
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
    elseif solver=="cg"
        cg_p  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        cg_r  = zeros(Float64, SZ[1], SZ[2], SZ[3])
        cg_ax = zeros(Float64, SZ[1], SZ[2], SZ[3])
        cg_ap = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end

    F = open("log.txt", "w")
    conditions(F, SZ, Δh, solver, smoother)

    if solver=="sor"
        if mode==1
            Cartesian.solveSOR!(θ, λ, b, mask, Δh, Constant.ω, F, itr_tol)
        elseif mode==2
            NonUniform.solveSOR!(θ, SZ, λ, b, mask, Δh, Constant.ω, Z, ΔZ, z_st, z_ed, F, itr_tol)
        elseif mode==3
            NonUniformII.solveSOR!(θ, SZ, λ, b, mask, Δh, Constant.ω, Z, ΔZ, z_st, z_ed, F, itr_tol)
        elseif mode==4
            CartesianII.solveSOR!(θ, λ, b, mask, Δh, Constant.ω, F, itr_tol)
        end
    elseif solver=="jacobi"
        if mode==1
            Cartesian.solveJACOBI!(θ, λ, b, mask, wk, Δh, Constant.ω, F, itr_tol)
        elseif mode==2
            NonUniform.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, Constant.ω, Z, ΔZ, z_st, z_ed, F, itr_tol)
        elseif mode==3
            NonUniformII.solveJACOBI!(θ, SZ, λ, b, mask, wk, Δh, Constant.ω, Z, ΔZ, z_st, z_ed, F, itr_tol)
        elseif mode==4
            CartesianII.solveJACOBI!(θ, λ, b, mask, wk, Δh, Constant.ω, F, itr_tol)
        end
    elseif solver=="pbicgstab"
        if mode==1
            Cartesian.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, ox, Δh, smoother, F, mode, itr_tol)
        elseif mode==2
            NonUniform.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, 
                ox, Δh, SZ, Z, ΔZ, z_st, z_ed, smoother, F, mode, itr_tol)
        elseif mode==3
            NonUniformII.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, 
                ox, Δh, SZ, Z, ΔZ, z_st, z_ed, smoother, F, mode, itr_tol)
        elseif mode==4
            CartesianII.PBiCGSTAB!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, ox, Δh, smoother, F, mode, itr_tol)
        end
    elseif solver=="cg"
        if mode==1
            Cartesian.CG!(θ, b, cg_p, cg_r, cg_ax, cg_ap, 
                λ, mask, Δh, F, mode, itr_tol)
        elseif mode==2
            NonUniform.CG!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, 
                ox, Δh, SZ, Z, ΔZ, z_st, z_ed, smoother, F, mode, itr_tol)
        elseif mode==3
            NonUniformII.CG!(θ, b, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, 
                pcg_s_, pcg_t_, λ, mask, wk, 
                ox, Δh, SZ, Z, ΔZ, z_st, z_ed, smoother, F, mode, itr_tol)
        elseif mode==4
            CartesianII.CG!(θ, b, cg_p, cg_r, cg_ax, cg_ap, 
                λ, mask, Δh, F, mode, itr_tol)
        end
    else
        println("solver error")
    end

    s = θ[2:SZ[1]-1, 2:SZ[2]-1, z_st:z_ed]
    min_val = minimum(s)
    max_val = maximum(s)
    @printf(F, "θmin=%e  θmax=%e  L2 norm of θ=%e\n", min_val, max_val, norm(s,2))

    close(F)
    
    # ログファイルから残差データを解析してconv_dataに追加
    if solver == "pbicgstab" || solver == "sor" || solver == "jacobi" || solver == "cg"
        parse_residuals_from_log!(conv_data, "log.txt")
    end
    
    # 収束履歴データを返す
    return conv_data
end

function Zcase1!(Z::Vector{Float64}, SZ)
    if SZ[3]!=15
        println("MZ must be 15")
        exit(0)
    end
    p = 0.005e-3
    Z[1] = 2.0*modelA.zm0-modelA.zm1
    Z[2] = modelA.zm0
    Z[3] = modelA.zm1
    Z[4] = modelA.zm2
    Z[5] = modelA.zm3
    Z[6] = modelA.zm4
    Z[7] = modelA.zm5
    Z[8] = modelA.zm6
    Z[9] = modelA.zm7
    Z[10]= modelA.zm8
    Z[11]= modelA.zm9
    Z[12]= modelA.zm10
    Z[13]= modelA.zm11
    Z[14]= modelA.zm12
    Z[15]= 2.0*modelA.zm12-modelA.zm11
end

function Zcase2!(Z::Vector{Float64}, SZ)
    if SZ[3]!=33
        println("MZ must be 33")
        exit(0)
    end
    p = 0.005e-3
    Z[1] = modelA.zm0 - p
    Z[2] = modelA.zm0
    Z[3] = modelA.zm0 + p
    Z[4] = modelA.zm1 - p
    Z[5] = modelA.zm1
    Z[6] = modelA.zm1 + p
    Z[7] = modelA.zm2 - p
    Z[8] = modelA.zm2
    Z[9] = modelA.zm2 + p
    Z[10]= modelA.zm3 - p
    Z[11]= modelA.zm3
    Z[12]= modelA.zm4
    Z[13]= modelA.zm4 + p
    Z[14]= modelA.zm5 - p
    Z[15]= modelA.zm5
    Z[16]= modelA.zm5 + p
    Z[17]= modelA.zm6 - p
    Z[18]= modelA.zm6
    Z[19]= modelA.zm7
    Z[20]= modelA.zm7 + p
    Z[21]= modelA.zm8 - p
    Z[22]= modelA.zm8
    Z[23]= modelA.zm8 + p
    Z[24]= modelA.zm9 - p
    Z[25]= modelA.zm9 
    Z[26]= modelA.zm10
    Z[27]= modelA.zm10 + p
    Z[28]= modelA.zm11 - p
    Z[29]= modelA.zm11
    Z[30]= modelA.zm11 + p
    Z[31]= modelA.zm12 - p
    Z[32]= modelA.zm12
    Z[33]= modelA.zm12 + p
end

function Zcase3!(Z::Vector{Float64}, SZ, ox, dz)
    for k in 1:SZ[3]
        Z[k] = ox[3] + (k-2)*dz
    end
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
        #Zcase1!(Z, SZ)
        Zcase2!(Z, SZ)
        #Zcase3!(Z, SZ, ox, dz)
    end

    for k in 2:mz-1
        ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])
    end

    if mode==3 || mode==2
        ΔZ[2] = 0.5*ΔZ[2]
        ΔZ[mz-1] = 0.5*ΔZ[mz-1]
    end
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
function dif(θ, exact)
    SZ = size(θ)
    d=0.0
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
       s = (θ[i,j,k]-exact[i,j,k])
       d = d + s*s
    end
    return sqrt(d)
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
function q3d(m_mode::Int, NXY::Int, NZ::Int, solver::String="sor", smoother::String=""; epsilon::Float64=1.0e-6)
    global mode = m_mode
    global itr_tol = epsilon

    MX = MY = NXY + 2  # Number of CVs including boundaries
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
    elseif mode==3 || mode==4
        dh = 1.2e-3 / NXY
    end
    dh = round(dh,digits=8) #4.9999.... >> 5.0にしたい

    SZ = (MX, MY, MZ)
    Δh = (dh, dh, dh)
    ox = (0.0, 0.0, 0.0) #原点を仮定

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

    ID = zeros(UInt8, SZ[1], SZ[2], SZ[3]) # mode=3のときのみ有効

    @time preprocess!(λ, Z, ΔZ, ox, Δh, ID)

    if mode==3
        plot_slice_xz_nu(1, λ, 0.3e-3, SZ, ox, Δh, Z, "alpha3.png", "α")
    elseif mode==4
        plot_slice_xz(1, mode, λ, Z, 0.3e-3, SZ, ox, Δh, "alpha4.png", "α")
    end

    if mode<=2
        exact = zeros(Float64, SZ[1], SZ[2], SZ[3])
    end 

    if mode==1
        Cartesian.exact_solution!(exact, ox, Δh)
        plot_slice_xz(1, mode, exact, Z, 0.5, SZ, ox, Δh, "exact.png", "Exact")
    elseif mode==2
        NonUniform.exact_solution!(exact, SZ, ox, Δh, Z)
        plot_slice_xz_nu(1, exact, 0.5, SZ, ox, Δh, Z, "exact_nu.png", "Exact")
    end
    
    θ = zeros(Float64, SZ[1], SZ[2], SZ[3])
    θ .= Constant.θ_amb # 初期温度設定

    @time conv_data = main(ox, Δh, θ, Z, ΔZ, ID, λ, solver, smoother)

    
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
        plot_slice_xz_nu(2, θ, 0.5, SZ, ox, Δh, Z, "p_nu.png", "solution")
        plot_slice_xz_nu(2, abs.(θ-exact), 0.5, SZ, ox, Δh, Z, "diff_nu.png", "diff")
        println(" L2 norm of θ-exact=",dif(θ, exact))
    elseif mode==3
        plot_slice_xz_nu(2, θ, 0.3e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.3.png")
        plot_slice_xz_nu(2, θ, 0.4e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.4.png")
        plot_slice_xz_nu(2, θ, 0.5e-3, SZ, ox, Δh, Z, "temp3_xz_nu_y=0.5.png")
        plot_slice_xy_nu(2, θ, 0.18e-3, SZ, ox, Δh, Z, "temp3_xy_nu_z=0.18.png")
        plot_slice_xy_nu(2, θ, 0.33e-3, SZ, ox, Δh, Z, "temp3_xy_nu_z=0.33.png")
        plot_slice_xy_nu(2, θ, 0.48e-3, SZ, ox, Δh, Z, "temp3_xy_nu_z=0.48.png")
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
    if solver == "pbicgstab" || solver == "sor" || solver == "jacobi" || solver == "cg"
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

end

if abspath(PROGRAM_FILE) == @__FILE__
  q3d(3, 240, 31, "pbicgstab", "gs", epsilon=1.0e-4)
  #q3d(3, 120, 31, "pbicgstab", "gs", epsilon=1.0e-4)
  #q3d(4, 240, 120, "pbicgstab", "gs", epsilon=1.0e-4) 
  #q3d(1, 25, 25, "pbicgstab", "gs", epsilon=1.0e-8)
  #q3d(1, 25, 25, "cg", epsilon=1.0e-8)
  #q3d(2, 25, 25, "pbicgstab", "gs", epsilon=1.0e-8)
end
q3d(1, 25, 25, "pbicgstab")
#q3d(2, 25, 25, "sor")
#q3d(3, 240, 121, "pbicgstab", "gs")
#q3d(4, 240, 120, "pbicgstab", "gs", epsilon=1.0e-4)
#q3d(4, 240, 120, "cg", epsilon=1.0e-4) 
#q3d(3, 240, 31, "pbicgstab", "gs", epsilon=1.0e-4)
#q3d(1, 25, 25, "cg")