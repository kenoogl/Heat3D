# common constant
module Constant

const ItrMax = 8000
const tol    = 1.0e-6
const FloatMin = 1.0e-37
const ω      = 1.0

# 境界条件
const Q_src = 9.5374e4 # 1.6x10^11 [W/m^3] / (\rho C)_silicon > [K/s]

end # end of module Constant


"""
Harmonic mean
@param a left value
@param b right value
@param ma mask for left
@param mb mask for right
"""
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))


function get_backend(par::String)
    return (par == "thread") ? ThreadedEx() : SequentialEx()
end
