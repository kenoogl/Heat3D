# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Heat3Dは3次元熱拡散シミュレーションのプロジェクトです。主にJuliaで実装されており、熱伝導方程式の数値解法を複数のアルゴリズムで提供します。Vinokurストレッチング関数のJulia実装も含まれています。

## Core Architecture

### ディレクトリ構成
- `base/`: 基本的な数値解法アルゴリズム（SOR、BiCGSTAB）
- `q3d/`: 3次元四層構造の熱解析モデル（NonUniform格子対応済み）
- `Vinokur_Julia/`: Vinokurストレッチング関数のJulia実装
- `vinokur_code/`: C++実装の格子生成コード（元実装）
- `vinokur_timing/`: 性能比較結果

### 主要な数値解法
1. **SOR法** (Successive Over-Relaxation): `heat3d.jl`, `heat3d_NonUniform.jl`
2. **Red-Black SOR法** (RB-SOR): 並列化対応のSOR法、NonUniform格子対応済み
3. **BiCGSTAB法** (Bi-Conjugate Gradient Stabilized): `pbicgstab.jl`
4. **CG法** (Conjugate Gradient): 共役勾配法
5. **Vinokurストレッチング関数**: 格子点分布計算（Julia実装）

### データ構造
- `Array{Float64,3}`: 3次元温度場、物性値配列
- `SZ`: 配列サイズ `(nx, ny, nz)`
- `Δh`: 格子間隔 `(dx, dy, dz)`
- `ox`: 原点座標 `(x0, y0, z0)`

## Development Commands

### Julia実行
```bash
julia base/heat3d.jl          # 基本的な熱拡散計算
julia base/heat3d_bicg.jl     # BiCGSTAB法による計算
julia q3d/modelA.jl           # 3次元四層モデルの実行
julia q3d/heat3d_nu.jl        # 統合された3次元熱拡散（境界条件システム対応）
julia q3d/boundary_example.jl # 境界条件設定の使用例
julia Vinokur_Julia/test_stretch.jl    # Vinokurストレッチング関数テスト
julia vinokur_timing/benchmark.jl      # C++とJuliaの性能比較
julia q3d/demo_nonuniform_plot.jl      # NonUniform格子可視化デモ
```

### 新しい実行例（v2.0境界条件システム）
```julia
# Mode別の実行例
q3d(1, 25, 25, "sor", epsilon=1.0e-4)           # Cartesian立方体問題
q3d(2, 25, 25, "pbicgstab", "gs", epsilon=1.0e-8)  # NonUniform立方体問題  
q3d(3, 240, 31, "pbicgstab", "gs", epsilon=1.0e-4)  # NonUniform IC問題
q3d(4, 240, 120, "cg", epsilon=1.0e-4)         # Cartesian IC問題
```

### C++コンパイル（格子生成）
```bash
cd vinokur_code
make CXX=g++                  # g++コンパイラを使用（推奨）
./test 101 0.0 1.0 0.01 0.01  # 実行例（101ノード、均等分布）
make clean                    # クリーンアップ
```

## Code Conventions

### Julia
- 関数名: snake_case（例: `boundary_condition!`）
- 配列操作: in-place更新には`!`を付与
- インデントは2スペース
- 物理パラメータは`const`で定義

### 物性値設定
材料IDによる物性値管理：
- ID=1: 銅（TSV）
- ID=2: シリコン
- ID=3: はんだ（bump）
- ID=4: PCB基板（FR4）
- ID=5: ヒートシンク（A1060）
- ID=6: アンダーフィル樹脂
- ID=7: 電源グリッド

### 境界条件システム（v2.0対応）
**境界条件タイプ:**
- `ISOTHERMAL`: 等温条件（Dirichlet境界）- mask=0.0、温度固定
- `HEAT_FLUX`: 熱流束条件（Neumann境界）- 指定熱流束、断熱条件は流束=0
- `CONVECTION`: 熱伝達条件（Robin境界）- 熱伝達係数と周囲温度指定

**構造体ベース設計:**
```julia
struct BoundaryCondition
    type::BoundaryType
    temperature::Float64
    heat_flux::Float64  
    heat_transfer_coefficient::Float64
    ambient_temperature::Float64
end
```

**Mode別の境界条件パターン:**
- Mode1,2: 立方体問題（6面等温0K + Z面三角関数分布）
- Mode3: IC問題（NonUniform、熱伝達側面 + PCB等温底面 + 熱伝達上面）
- Mode4: IC問題（Cartesian、熱伝達側面 + PCB等温底面 + 熱伝達上面）

## Key Functions

### 基本計算関数（base/）
- `sor!()`: SOR法による反復計算
- `rbsor!()`: Red-Black SOR法（並列化対応）
- `residual()`: 残差計算
- `boundary_condition!()`: 境界条件設定

### NonUniform格子解法関数（q3d/heat3d_NonUniform.jl）
- `solveSOR!()`: NonUniform格子対応SOR法求解
- `rbsor!()`: NonUniform格子対応Red-Black SOR法（新設）
- `rbsor_core!()`: RB-SORカーネル関数（マルチカラー並列対応）
- `resSOR()`: NonUniform格子残差計算
- `PBiCGSTAB!()`: NonUniform格子対応BiCGSTAB法

### 3次元モデル関数（q3d/）
- `fillID!()`: 材料ID配列の生成
- `setLambda!()`: 物性値配列の設定
- `genZ!()`: Z方向非等間隔格子生成
- `model_test()`: 統合テスト実行

### 境界条件管理（q3d/boundary_conditions.jl）
- `isothermal_bc()`, `heat_flux_bc()`, `adiabatic_bc()`, `convection_bc()`: 境界条件作成関数
- `create_boundary_conditions()`: 6面境界条件セット作成
- `apply_boundary_conditions!()`: 境界条件の配列への適用
- `print_boundary_conditions()`: 境界条件情報の表示
- Mode別境界条件設定: `set_mode1_bc_parameters()`, `set_mode2_bc_parameters()`, etc.

### NonUniform格子対応関数（q3d/plotter.jl）
- `plot_slice_nu()`: NonUniform格子断面可視化（物理座標系）
- `plot_slice2_nu()`: NonUniform格子全セル可視化（物理座標系）
- `find_j()`: Y座標からグリッドインデックス変換

### Vinokurストレッチング関数（Vinokur_Julia/）
- `stretch()`: 格子点分布計算（Julia実装）
- `stretching!()`: 内部ストレッチング計算
- `plot_grid_distribution()`: 格子分布可視化

## Testing and Visualization

### プロット出力
- `plot_slice()`: 断面可視化（Cartesian格子用）
- `plot_slice_nu()`: NonUniform格子対応断面可視化（物理座標系）
- `plot_slice2()`: 全セル断面可視化（Cartesian格子用）
- `plot_slice2_nu()`: NonUniform格子対応全セル断面可視化（物理座標系）
- `id_xy()`, `id_yz()`: 材料分布可視化
- 出力形式: PNG画像

**NonUniform格子可視化の重要な修正点:**
- 従来の`plot_slice()`は格子インデックスを座標軸として使用（物理的に不正確）
- 修正版`plot_slice_nu()`は実際の物理座標（Z座標値）を使用
- アスペクト比1:1設定により幾何学的正確性を確保

### 収束判定
- 相対残差: `tol = 1.0e-8`
- 最大反復数: `ItrMax = 5000-8000`
- Vinokurストレッチング: `EPS = 1.0e-5`, `ITER_MAX = 1000`

## Performance Comparison

### C++ vs Julia ベンチマーク結果（Vinokurストレッチング関数）

**実行時間比較（均等分布格子）:**

| ノード数 | C++時間 | Julia時間 | 速度比 |
|---------|---------|-----------|--------|
| 51      | 3.4ms   | 0.001ms   | 3,400x |
| 101     | 1.7ms   | 0.001ms   | 1,700x |
| 201     | 1.2ms   | 0.001ms   | 1,200x |
| 501     | 1.2ms   | 0.005ms   | 240x   |

**詳細分析（101ノード）:**
- Julia平均: 0.0003ms（300ナノ秒）
- C++平均: 1.231ms
- **総合速度比: Juliaが4,364倍高速**

**メモリ効率:**
- Juliaメモリ使用量: 5.34 KB（201ノード）
- アロケーション回数: 4回

## Dependencies

### Julia Packages
- `Printf`: 出力フォーマット
- `Plots`: 可視化
- `BenchmarkTools`: 性能測定（ベンチマーク用）
- `LinearAlgebra`: 線形代数演算（norm計算等）

### C++ Requirements
- g++コンパイラ（pgc++から変更推奨）
- 数学ライブラリ（標準）

## 重要な改善点

### コードベース統合とリファクタリング（v2.0）
**主要な統合作業:**
1. **境界条件システムの統合**: 新しい構造体ベースの境界条件システム
2. **モジュール統合**: CartesianII → Cartesian、NonUniformII → NonUniform
3. **SZ引数除去**: 配列次元情報の直接取得によるコード簡素化

**新しい境界条件アーキテクチャ:**
- `BoundaryConditions`モジュール: 3種類×6面の境界条件管理
- 境界条件タイプ: `ISOTHERMAL`, `HEAT_FLUX`, `CONVECTION`
- Mode別境界条件設定関数: `set_mode1_bc_parameters()`, `set_mode2_bc_parameters()`, etc.

**統合されたファイル:**
- `q3d/boundary_conditions.jl`: 新設された境界条件管理モジュール
- `q3d/boundary_example.jl`: 境界条件使用例
- `q3d/Zcoord.jl`: Z座標管理（genZ!関数統合）
- `q3d/heat3d_nu.jl`: Mode別境界条件設定の統一インターフェース

### Mode2専用境界条件関数の追加
**問題**: Mode2（NonUniform格子）でZ面境界条件の設定位置が不正確
- Cartesian格子: k=1, k=SZ[3]（ゴースト点）に設定
- NonUniform格子: k=2, k=SZ[3]-1（内点）に設定が必要

**解決策**: 
- `bc_cube_nu!()`関数を新設（Mode2専用）
- Mode別の境界条件適用: `apply_isothermal!()`でmode引数追加
- Z面境界でのmode判定による正確な位置設定

**各Mode境界条件詳細:**
- Mode1: Cartesian立方体（6面等温0K → Z面三角関数分布で上書き）
- Mode2: NonUniform立方体（6面等温0K → Z面三角関数分布で上書き）  
- Mode3: NonUniform IC（X,Y面=熱伝達HT_side、Z下面=PCB等温、Z上面=熱伝達HT_top）
- Mode4: Cartesian IC（X,Y面=熱伝達HT_side、Z下面=PCB等温、Z上面=熱伝達HT_top）

### Red-Black SOR法のNonUniform格子対応
**新機能**: NonUniform格子でのRed-Black SOR法実装
- `rbsor!()`: 2色（Red/Black）マルチカラー並列化SOR法
- `rbsor_core!()`: カーネル関数、`@simd`最適化対応
- チェッカーボードパターン: `2+mod(k+j+color,2):2:SZ[1]-1`
- NonUniform格子のZ方向半セル処理に対応

**性能向上**: 従来のSORからRB-SORに切り替え（Line 125）
```julia
# 従来: res = sor!(θ, λ, b, mask, Δh, ω, Z, ΔZ, z_range, HF, HT) / res0
# 新版: res = rbsor!(θ, SZ, λ, b, mask, Δh, ω, z_range, HF, HT) / res0
```

### プロッタ可視化エラー修正
**問題**: 数値解が0の場合にlog10(0)=-Infエラー
**解決策**: `max_val<=0.0`チェックを追加してmin_val同様に処理

**影響するファイル**: 
- `q3d/heat3d_NonUniform.jl`: RB-SOR法追加、デフォルトをRB-SORに変更
- `q3d/plotter.jl`: エラー処理追加、NonUniform対応関数統合
- `q3d/heat3d_nu.jl`: Mode別境界条件とbc_cube_nu!統合
- `q3d/boundary_conditions.jl`: 新設境界条件システム