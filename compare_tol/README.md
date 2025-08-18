# 許容誤差比較実験

Heat3Dシミュレーションにおける異なる収束許容誤差値（10^-3 から 10^-8）での計算精度を比較する実験です。

## 実験概要

- **対象許容誤差**: 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8
- **計算パラメータ**: q3d(3, 240, 31, "pbicgstab", "gs", ε)
- **精度評価**: export_zline_csv関数出力のZ方向温度分布を比較
- **評価指標**: 温度差（絶対値・相対値）、統計的指標（L2ノルム等）

## ファイル構成

```
compare_tol/
├── run_tolerance_test.jl      # 個別許容誤差テスト実行スクリプト
├── compare_csv_results.jl     # CSV結果比較・レポート生成スクリプト
├── README.md                  # このファイル
└── [実行結果ファイル]
    ├── temp3Z_ctr_eps*.csv    # 中心線温度分布
    ├── temp3Z_tsv_eps*.csv    # TSV線温度分布
    ├── convergence_*_eps*.csv # 収束履歴
    └── tolerance_comparison_report.md # 比較レポート
```

## 実行方法

### ステップ1: 個別計算の実行

各許容誤差値で計算を実行します：

```bash
cd /Users/Daily/Development/Heat3D/compare_tol

# 許容誤差 1.0e-3 での計算
julia run_tolerance_test.jl 1.0e-3

# 許容誤差 1.0e-4 での計算  
julia run_tolerance_test.jl 1.0e-4

# 許容誤差 1.0e-5 での計算
julia run_tolerance_test.jl 1.0e-5

# 許容誤差 1.0e-6 での計算
julia run_tolerance_test.jl 1.0e-6

# 許容誤差 1.0e-7 での計算
julia run_tolerance_test.jl 1.0e-7

# 許容誤差 1.0e-8 での計算
julia run_tolerance_test.jl 1.0e-8
```

### ステップ2: 結果の比較・レポート生成

全ての計算が完了したら、結果を比較してレポートを生成します：

```bash
julia compare_csv_results.jl
```

## 出力ファイル

### 個別計算の出力

各許容誤差値での計算により以下のファイルが生成されます：

- `temp3Z_ctr_eps<値>.csv`: 中心線（x=0.6mm, y=0.6mm）のZ方向温度分布
- `temp3Z_tsv_eps<値>.csv`: TSV線のZ方向温度分布  
- `convergence_pbicgstab_mode3_240x31_gs_eps<値>.csv`: 収束履歴データ
- `convergence_pbicgstab_mode3_240x31_gs_eps<値>.png`: 収束履歴グラフ

ファイル名例：
- `temp3Z_ctr_eps1_0e-4.csv` (許容誤差 1.0e-4)
- `temp3Z_tsv_eps1_0em6.csv` (許容誤差 1.0e-6)

### 比較レポート

`tolerance_comparison_report.md` に以下の内容が含まれます：

1. **温度統計**: 各許容誤差での最小・最大・平均温度
2. **精度比較**: 参照解（最厳密許容誤差）との温度差
3. **収束解析**: 推奨許容誤差値の提案
4. **結論**: 用途別の許容誤差推奨値

## スクリプトの詳細

### run_tolerance_test.jl

- 指定された許容誤差値でq3d()を実行
- 出力ファイル名に許容誤差値を自動付与
- 実行時間を計測
- エラーハンドリング機能

### compare_csv_results.jl

- 各許容誤差のCSVファイルを自動検出・読み込み
- 統計解析（最小・最大・平均・標準偏差・L2ノルム）
- 参照解との差分計算（絶対差・相対差・RMS差）
- Markdownレポート自動生成

## 注意事項

1. **計算時間**: 厳密な許容誤差ほど計算時間が長くなります
2. **メモリ使用量**: 240x240x31グリッドは大容量メモリを要求します
3. **ファイル管理**: 各実行で複数のファイルが生成されるため注意が必要
4. **実行順序**: 個別計算→比較分析の順序で実行してください

## トラブルシューティング

### エラー: "File not found"
- CSVファイルが生成されていない可能性があります
- run_tolerance_test.jlが正常終了したか確認してください

### エラー: "Memory allocation failed"
- システムメモリが不足している可能性があります
- より小さなグリッドサイズでテストしてください

### 計算が終わらない
- 厳密な許容誤差（1.0e-7, 1.0e-8）は長時間要する場合があります
- Ctrl+Cで中断し、より緩い許容誤差から開始してください

## 期待される結果

- **計算精度**: 厳密な許容誤差ほど高精度な解が得られる
- **計算コスト**: 許容誤差と計算時間は反比例関係
- **推奨値**: 工学解析には1.0e-4〜1.0e-5、研究用途には1.0e-6〜1.0e-7

---

*更新日: 2025-08-18*