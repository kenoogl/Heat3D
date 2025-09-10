# 境界条件設定の使用例

include("boundary_conditions.jl")
using .BoundaryConditions

"""
Mode1用の境界条件（Cartesian格子の立方体問題）
Z面に三角関数分布の等温条件、X,Y面は等温0K
"""
function get_mode1_boundary_conditions()
    # 各面の境界条件を定義
    x_minus_bc = isothermal_bc(0.0)           # X軸負方向面: 等温0K
    x_plus_bc  = isothermal_bc(0.0)           # X軸正方向面: 等温0K
    y_minus_bc = isothermal_bc(0.0)           # Y軸負方向面: 等温0K
    y_plus_bc  = isothermal_bc(0.0)           # Y軸正方向面: 等温0K
    z_minus_bc = isothermal_bc(0.0)           # Z軸負方向面: 等温0K (後で三角関数分布で上書き)
    z_plus_bc  = isothermal_bc(0.0)           # Z軸正方向面: 等温0K (後で三角関数分布で上書き)
    
    # 境界条件セットを作成
    return create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
Mode4用の境界条件（IC package問題）
Z下面: PCB温度、Z上面: 周囲温度、側面: 周囲温度
"""
function get_mode4_boundary_conditions()
    # 定数値を取得（const.jlから）
    θ_pcb = 300.0   # Constant.θ_pcb 
    θ_amb = 300.0   # Constant.θ_amb
    
    # 各面の境界条件を定義
    x_minus_bc = isothermal_bc(θ_amb)         # X軸負方向面: 周囲温度
    x_plus_bc  = isothermal_bc(θ_amb)         # X軸正方向面: 周囲温度
    y_minus_bc = isothermal_bc(θ_amb)         # Y軸負方向面: 周囲温度  
    y_plus_bc  = isothermal_bc(θ_amb)         # Y軸正方向面: 周囲温度
    z_minus_bc = isothermal_bc(θ_pcb)         # Z軸負方向面: PCB温度
    z_plus_bc  = isothermal_bc(θ_amb)         # Z軸正方向面: 周囲温度
    
    # 境界条件セットを作成
    return create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
Mode3用の境界条件（NonUniform格子のIC問題）  
Z下面: PCB温度、Z上面: 熱伝達、側面: 断熱
"""
function get_mode3_boundary_conditions()
    # 定数値を取得
    θ_pcb = 300.0      # Constant.θ_pcb
    θ_amb = 300.0      # Constant.θ_amb
    HT_top = 2.98e-4   # Constant.HT_top
    
    # 各面の境界条件を定義
    x_minus_bc = adiabatic_bc()                        # X軸負方向面: 断熱
    x_plus_bc  = adiabatic_bc()                        # X軸正方向面: 断熱
    y_minus_bc = adiabatic_bc()                        # Y軸負方向面: 断熱
    y_plus_bc  = adiabatic_bc()                        # Y軸正方向面: 断熱
    z_minus_bc = isothermal_bc(θ_pcb)                  # Z軸負方向面: PCB温度
    z_plus_bc  = convection_bc(HT_top, θ_amb)          # Z軸正方向面: 熱伝達
    
    # 境界条件セットを作成
    return create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
カスタム境界条件の例
X面: 熱流束、Y面: 断熱、Z面: 等温と熱伝達の組み合わせ
"""
function get_custom_boundary_conditions()
    # 各面の境界条件を定義
    x_minus_bc = heat_flux_bc(1000.0)                 # X軸負方向面: 1000 W/m²の熱流束
    x_plus_bc  = heat_flux_bc(-500.0)                 # X軸正方向面: -500 W/m²の熱流束（冷却）
    y_minus_bc = adiabatic_bc()                       # Y軸負方向面: 断熱
    y_plus_bc  = adiabatic_bc()                       # Y軸正方向面: 断熱
    z_minus_bc = isothermal_bc(350.0)                 # Z軸負方向面: 350Kの等温
    z_plus_bc  = convection_bc(25.0, 293.15)          # Z軸正方向面: h=25, T∞=293.15Kの熱伝達
    
    # 境界条件セットを作成
    return create_boundary_conditions(x_minus_bc, x_plus_bc,
                                      y_minus_bc, y_plus_bc,
                                      z_minus_bc, z_plus_bc)
end

"""
境界条件を適用する例
"""
function example_apply_boundary_conditions()
    # サンプル配列を作成
    SZ = (10, 10, 10)
    θ = ones(Float64, SZ...) * 300.0    # 初期温度300K
    λ = ones(Float64, SZ...) * 1.0      # 初期熱伝導率1.0
    mask = ones(Float64, SZ...)         # 初期マスク値1.0
    b = zeros(Float64, SZ...)         # 右辺項
    Δh = (0.1, 0.1, 0.1)               # セル幅
    
    # Mode4の境界条件を取得
    bc_set = get_mode4_boundary_conditions()
    
    # 境界条件情報を表示
    print_boundary_conditions(bc_set)
    
    # 境界条件を適用
    apply_boundary_conditions!(θ, λ, mask, b, bc_set, Δh)
    
    println("\n境界条件適用後:")
    println("θ[1,5,5] (X-面) = $(θ[1,5,5])")
    println("θ[10,5,5] (X+面) = $(θ[10,5,5])") 
    println("θ[5,1,5] (Y-面) = $(θ[5,1,5])")
    println("θ[5,10,5] (Y+面) = $(θ[5,10,5])")
    println("θ[5,5,1] (Z-面) = $(θ[5,5,1])")
    println("θ[5,5,10] (Z+面) = $(θ[5,5,10])")
    
    println("\nmask[1,5,5] (X-面) = $(mask[1,5,5])")
    println("mask[5,5,1] (Z-面) = $(mask[5,5,1])")
    println("mask[5,5,10] (Z+面) = $(mask[5,5,10])")
end

# 使用例の実行
if abspath(PROGRAM_FILE) == @__FILE__
    example_apply_boundary_conditions()
end