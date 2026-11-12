# 追加実験の目次

[English](README.md)

2026年11月の実験をテーマ別にまとめています。各レポートに計算条件、結果、検証方法、保存データの説明があります。4量子ビットの基準実験と8月の追加解析は、それぞれのReleaseから再現用ZIPを取得できます。[7月の基準実験](https://github.com/Riku-Tono/four-qubit-open-vs-ring-dynamics/releases/tag/v0.1.0-n4-fixed) · [8月の追加解析](https://github.com/Riku-Tono/four-qubit-open-vs-ring-dynamics/releases/tag/v0%2C2)

- [端点横磁場への応答](#端点横磁場への応答)
- [継ぎ目近傍の波形と短時間展開](#継ぎ目近傍の波形と短時間展開)
- [連結相関の距離則と相殺](#連結相関の距離則と相殺)
- [一様結合強度と相関差](#一様結合強度と相関差)
- [逆位相ペアと局所回転](#逆位相ペアと局所回転)
- [4サイトの同位相と初期角度誤差](#4サイトの同位相と初期角度誤差)
- [全64グラフの有限ショット識別](#全64グラフの有限ショット識別)

## 実行方法

隣接フォルダのコードや保存済みデータを参照する実験があるため、リポジトリ全体を取得してください。各表の実行順に従い、対応する実験フォルダで各項目の先頭に `python -B` を付けて実行します。スクリプトが結果ファイルを再生成するため、保存済みの結果も参照したい場合は作業用コピーで実行してください。

たとえば、リポジトリのルートから開放鎖の距離則を検証する手順です。

```sh
cd experiments/open-chain-distance-law-2026-11-09
python -B verify.py
```

数値計算には主に NumPy 2 系を使います。`plot.py` と端点磁場の最初の3実験の `run.py` には Matplotlib も必要です。標準ライブラリで実行できる実験は下に記載しています。実験ごとの条件や入力は各レポートを参照してください。

## 端点横磁場への応答

| レポート | 調べたこと | 実行順 |
| --- | --- | --- |
| [4量子ビットの端点磁場と反ユニタリ対称性](endpoint-field-2026-11-09/REPORT.md) | 端点の横磁場の符号と相関差の関係 | `run.py` |
| [逆方向の端点磁場と継ぎ目シェアのスケーリング](endpoint-field-scaling-2026-11-09/REPORT.md) | 磁場強度を変えたときの継ぎ目ペアの寄与 | `run.py` |
| [4量子ビットでの継ぎ目シェアの線形応答](endpoint-field-zero-slope-2026-11-09/REPORT.md) | 弱い磁場に対する一次応答の外挿 | `run.py` |
| [6量子ビットでの継ぎ目シェアの線形応答](endpoint-field-six-qubit-2026-11-09/REPORT.md) | 4量子ビットとの比較と積分刻みの照合 | `run.py` |

## 継ぎ目近傍の波形と短時間展開

| レポート | 調べたこと | 実行順 |
| --- | --- | --- |
| [全ペア積分の4→6量子ビット差と観測時間](seam-pair-time-window-2026-11-09/REPORT.md) | 全体の差に効くペアと観測窓による変化 | `run.py` |
| [継ぎ目をまたぐ2組の相関波形](seam-neighbor-waveforms-2026-11-09/REPORT.md) | 開放鎖とリングの波形および二乗積分の内訳 | `run.py` |
| [継ぎ目をまたぐ2組の短時間展開](seam-neighbor-short-time-2026-11-09/REPORT.md) | 厳密級数の先頭項と経路長の関係 | `run.py`（標準ライブラリ） |

## 連結相関の距離則と相殺

この節のスクリプトは Python の標準ライブラリと保存済みの実験ファイルを使います。

| レポート | 調べたこと | 実行順 |
| --- | --- | --- |
| [8量子ビットでの経路長と連結相関の先頭次数](seam-neighbor-eight-qubit-2026-11-09/REPORT.md) | 8量子ビットの開放鎖・リングでの先頭次数 | `run.py` |
| [開放 Ising 鎖の短時間距離則](open-chain-distance-law-2026-11-09/REPORT.md) | 距離と XX・YY・ZZ の先頭次数の関係 | `verify.py` |
| [弱結合リングの短時間距離則](weak-ring-distance-law-2026-11-09/REPORT.md) | 2経路の係数と弱結合リングでの距離則 | `coefficients.py` → `verify.py` |
| [正反対ペアの先頭項が相殺する弱結合](antipodal-cancellation-2026-11-09/REPORT.md) | 偶数サイズのリングで YY と ZZ それぞれの先頭項が消える条件 | `formulas.py` → `verify.py` |

## 一様結合強度と相関差

| レポート | 調べたこと | 実行順 |
| --- | --- | --- |
| [一様結合強度と直線・輪の連結相関差](uniform-coupling-sweep-2026-11-10/REPORT.md) | 結合強度によるペア別・軸別の相関差の変化 | `run.py` → `plot.py` |
| [相関差が最大のペアが切り替わる結合強度](uniform-coupling-crossover-2026-11-10/REPORT.md) | 外側ペアと中央ペアの相関差の交点 | `run.py` → `plot.py` |
| [外側ペアと中央ペアの首位交代を時間帯で分解](crossover-time-window-2026-11-10/REPORT.md) | 観測終了時刻と相関差の順位の関係 | `run.py` → `plot.py` |

## 逆位相ペアと局所回転

| レポート | 調べたこと | 実行順 |
| --- | --- | --- |
| [初期状態 1000 からの占有率・振幅・位相](phase-alignment-2026-11-10/REPORT.md) | 開放鎖とリングにおける位相差と占有率の補完関係 | `run.py` → `plot.py` |
| [4量子ビットの逆位相制約](antiphase-constraints-2026-11-10/REPORT.md) | ペアの占有率の和が1になる関係と数値残差 | `run.py` |
| [0番サイトの Y 軸回転と逆位相制約](antiphase-y-rotation-2026-11-10/REPORT.md) | 初期状態を回転したときの残差と位相 | `run.py` → `endpoint_audit.py` |
| [正負の Y 軸回転で分ける逆位相ペアの残差](antiphase-y-signed-2026-11-10/REPORT.md) | 回転の符号に対する残差の偶成分・奇成分 | `run.py` → `plot.py` |
| [0番量子ビットの X 軸回転と Y 軸回転](antiphase-xy-rotation-2026-11-10/REPORT.md) | 回転軸による占有率和の残差の違い | `run.py` → `plot.py` |
| [逆位相ペアの残差を抑える回転軸](antiphase-axis-optimization-2026-11-10/REPORT.md) | 有限の観測窓で残差を最小化する軸 | `run.py` → `plot.py` |
| [短窓で決めた回転軸の時間外検証](antiphase-axis-holdout-2026-11-10/REPORT.md) | 短い観測窓で選んだ軸を後続の窓で検証 | `run.py` → `plot.py` |
| [逆位相ペアの残差を抑える軸の長時間極限](antiphase-axis-long-time-2026-11-10/REPORT.md) | 観測窓を延ばしたときの最適軸とその極限 | `run.py` → `plot.py` |

## 4サイトの同位相と初期角度誤差

| レポート | 調べたこと | 実行順 |
| --- | --- | --- |
| [4サイトの局所回転による占有率の同位相探索](four-site-in-phase-2026-11-10/REPORT.md) | 4サイトの占有率が同位相で振動する初期状態の探索 | `baseline.py` → `run.py select` → `run.py evaluate` → `plot.py` |
| [開放4サイトの占有率の同位相はどこまで続くか](open-chain-phase-persistence-2026-11-10/REPORT.md) | 選んだ状態を長時間の別の観測窓で検証 | `run.py evaluate` → `plot.py` |
| [同位相振動を支える周波数と角度誤差への頑健性](open-chain-phase-spectrum-robustness-2026-11-10/REPORT.md) | 周波数成分の分解と初期角度をずらしたときの変化 | `spectrum.py` → `majorana.py` → `robustness.py` → `plot.py` |
| [独立な初期角度誤差に対する同位相振動](open-chain-independent-angle-errors-2026-11-10/REPORT.md) | サイトごとに異なる角度誤差を与えた場合の位相 | `sweep.py` → `diagnostics.py` → `plot.py` |
| [PLV が 0.95 を下回る最初の誤差幅](open-chain-plv-onset-2026-11-10/REPORT.md) | 誤差の格子走査で位相同期指標のしきい値を初めて下回る条件 | `scan.py` → `diagnostics.py` → `plot.py` |

長時間同位相の検証は、公開済みの `selection.json` を使います。初期状態の探索からやり直す場合は、別の作業用コピーで `python -B run.py select` を先に実行してください。

## 全64グラフの有限ショット識別

読み出し誤り5%の条件で、99ショットでの識別を保証する設計と、52ショット以下を除外する下限を収録しています。測定設定と回数は事前に固定します。 [実験レポートと再計算手順](Fable-Codex-2026-11-13/REPORT.md) · [99ショットの監査ZIP](Fable-Codex-2026-11-13/Codex/four_qubit_p5_99_audit.zip)
