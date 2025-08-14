# common constant
module Constant

const ItrMax = 500 #8000
const tol    = 1.0e-6
const FloatMin = 1.0e-37
const ω      = 1.0

# 境界条件
const θ_amb = 300.0
const θ_pcb = 300.0
const HT_top = 2.98e-6 # 5 [W/(m^2 K)] / (\rho C)_silicon > [m/s]
const HT_side = 2.98e-6 # 5 [W/(m^2 K)] / (\rho C)_silicon > [m/s]
const Q_src = 9.5374e4# 1.6x10^11 [W/m^3] / (\rho C)_silicon > [K/s]
# 8x10^5 [W/m^2], 厚さ5μm > 1.6 x 10^11 [W/m^3]
# \rho_silicon= 2330.0 [kg/m^3]
# C_silicon= 720.0 [J/(Kg K)]

end # end of module Constant