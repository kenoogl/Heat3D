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
1. **SOR法** (Successive Over-Relaxation): `heat3d.jl`
2. **BiCGSTAB法** (Bi-Conjugate Gradient Stabilized): `pbicgstab.jl`
3. **Red-Black SOR**: 並列化対応のSOR法
4. **Vinokurストレッチング関数**: 格子点分布計算（Julia実装）

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
julia q3d/heat3d_nu.jl        # NonUniform格子3次元熱拡散（修正済み）
julia Vinokur_Julia/test_stretch.jl    # Vinokurストレッチング関数テスト
julia vinokur_timing/benchmark.jl      # C++とJuliaの性能比較
julia q3d/demo_nonuniform_plot.jl      # NonUniform格子可視化デモ
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

### 境界条件
- ディリクレ境界: `mask=0.0`、温度固定
- ノイマン境界: `λ=0.0`、断熱条件

## Key Functions

### 基本計算関数（base/）
- `sor!()`: SOR法による反復計算
- `rbsor!()`: Red-Black SOR法
- `residual()`: 残差計算
- `boundary_condition!()`: 境界条件設定

### 3次元モデル関数（q3d/）
- `fillID!()`: 材料ID配列の生成
- `setLambda!()`: 物性値配列の設定
- `genZ!()`: Z方向非等間隔格子生成
- `model_test()`: 統合テスト実行

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

### NonUniform格子可視化問題の解決（v1.1）
**問題**: `q3d/heat3d_nu.jl`でNonUniform格子の可視化が物理的に不正確
- 従来の`plot_slice()`は格子インデックスを座標軸として使用
- 結果: 物理座標系を反映しない歪んだ可視化

**解決策**: 
1. `plot_slice_nu()`関数を追加（物理座標系対応）
2. 実際のZ座標値を軸として使用
3. アスペクト比1:1設定で幾何学的正確性確保
4. `plotter.jl`に統合済み

**影響するファイル**: 
- `q3d/plotter.jl`: NonUniform対応関数統合
- `q3d/heat3d_nu.jl`: 修正版可視化関数呼び出し