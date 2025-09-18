module Zcoordinate
export Zcase1!, Zcase2!, Zcase3!, genZ!

include("../model/modelA.jl")
using Printf

"""
    read_grid_file(filename)

ASCIIファイルから格子点座標を読み込む

# Arguments
- `filename::String`: 入力ファイル名

# Returns
- `coord::Vector{Float64}`: 格子点座標配列
- `numNodes::Int`: 格子点数
"""
function read_grid_file(filename::String="grid.txt")
  if !isfile(filename)
    error("ファイルが見つかりません: $filename")
  end
  
  coord = Float64[]
  numNodes = 0
  
  open(filename, "r") do f
    # 1行目: 格子点数を読み込み
    line = readline(f)
    numNodes = parse(Int, strip(line))
    
    # 座標配列を初期化
    coord = zeros(Float64, numNodes)
    
    # 格子点データを読み込み
    for i in 1:numNodes
      line = readline(f)
      parts = split(strip(line))
      
      if length(parts) != 2
        error("ファイル形式エラー（行$(i+1)）: $line")
      end
      
      grid_index = parse(Int, parts[1])
      grid_coord = parse(Float64, parts[2])
      
      # 格子点番号の整合性チェック（1オリジン）
      if grid_index != i
        @warn "格子点番号が期待値と異なります: 期待値=$i, 実際値=$grid_index"
      end
      
      coord[i] = grid_coord
    end
  end
    
  println("格子点データ読み込み完了:")
  @printf("  ファイル: %s\n", filename)
  @printf("  格子点数: %d\n", numNodes)
  @printf("  座標範囲: [%.6f, %.6f]\n", minimum(coord), maximum(coord))
  
  return coord, numNodes
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
function genZ!(Z::Vector{Float64}, ΔZ::Vector{Float64}, SZ, ox, dz, mode)
    mz=SZ[3]
    if mode==1 || mode==4
        for k in 1:mz+1
            Z[k] = ox[3] + (k-2)*dz
        end
    elseif mode==2
        read_coord, numNodes = read_grid_file()
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
    println(ΔZ)
end

end # end of module