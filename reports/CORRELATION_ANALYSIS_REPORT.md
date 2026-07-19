# 二体相関解析：切替なし4bit直列対4bit円環

## 定義

一体量は `<A_i>`、二体raw量は `<A_i A_j>`、connected相関は `C^AA_ij=<A_i A_j>-<A_i><A_j>`（A=X,Y,Z）です。4×4表示では比較対象でない対角を0に固定します。

各チャネルのupper L1/L2は6個の `i<j` 成分に対する距離です。Frobenius距離は対称な4×4ゼロ対角行列に対し `sqrt(2) * upper_L2`、統合距離は `D_corr=sqrt(D_XX^2+D_YY^2+D_ZZ^2)` です。

## 主結果（有限N=4、0<=t<=3.688）

- XX: 最大upper-L2=0.604701641381 at t=3.340
- YY: 最大upper-L2=0.469912653201 at t=3.190
- ZZ: 最大upper-L2=0.665558112764 at t=3.688
- 統合D_corr最大=0.988619981363 at t=3.200
- 状態距離最大=0.546934944036 at t=3.462; この時のD_corr=0.974637690628
- D_corrが1e-08を最初に超える時刻: t=0.012
- 状態距離>0.5区間での平均チャネル距離: XX=0.564273, YY=0.40003, ZZ=0.63102

t=3.462 のconnected相関（open / closed / 差）:

| channel | pair | open | closed | open-closed |
|---|---:|---:|---:|---:|
| XX | 01 | -0.196609683 | 0.00676489627 | -0.203374579 |
| XX | 02 | -0.0235895718 | -0.0297264872 | 0.00613691543 |
| XX | 03 | 0.478921356 | 0.00676489627 | 0.47215646 |
| XX | 12 | -0.231523258 | 0.00676489627 | -0.238288154 |
| XX | 13 | -0.0235895718 | -0.0297264872 | 0.00613691543 |
| XX | 23 | -0.196609683 | 0.00676489627 | -0.203374579 |
| YY | 01 | 0.37779505 | 0.316782279 | 0.0610127711 |
| YY | 02 | 0.213084588 | 0.368699724 | -0.155615136 |
| YY | 03 | 0.0486186388 | 0.316782279 | -0.26816364 |
| YY | 12 | 0.56456502 | 0.316782279 | 0.247782741 |
| YY | 13 | 0.213084588 | 0.368699724 | -0.155615136 |
| YY | 23 | 0.37779505 | 0.316782279 | 0.0610127711 |
| ZZ | 01 | 0.416403264 | 0.464558696 | -0.0481554322 |
| ZZ | 02 | 0.145375256 | 0.456845621 | -0.311470366 |
| ZZ | 03 | 0.0238689717 | 0.464558696 | -0.440689725 |
| ZZ | 12 | 0.391416953 | 0.464558696 | -0.0731417431 |
| ZZ | 13 | 0.145375256 | 0.456845621 | -0.311470366 |
| ZZ | 23 | 0.416403264 | 0.464558696 | -0.0481554322 |

## 距離間の対応

- D_XX vs D_state: Pearson=0.904029, Spearman=0.997396, max=0.604702 at t=3.340, 状態距離最大との時刻差=-0.122
- D_YY vs D_state: Pearson=0.941816, Spearman=0.958606, max=0.469913 at t=3.190, 状態距離最大との時刻差=-0.272
- D_ZZ vs D_state: Pearson=0.952728, Spearman=0.976435, max=0.665558 at t=3.688, 状態距離最大との時刻差=0.226
- D_corr vs D_state: Pearson=0.953794, Spearman=0.990560, max=0.98862 at t=3.200, 状態距離最大との時刻差=-0.262
- D_onebody vs D_state: Pearson=0.969556, Spearman=0.998150, max=0.667066 at t=3.688, 状態距離最大との時刻差=0.226
- D_ZZ_connected vs D_state: Pearson=0.952728, Spearman=0.976435, max=0.665558 at t=3.688, 状態距離最大との時刻差=0.226
- D_all_connected vs D_state: Pearson=0.953794, Spearman=0.990560, max=0.98862 at t=3.200, 状態距離最大との時刻差=-0.262
- D_onebody_all_connected vs D_state: Pearson=0.981143, Spearman=0.998150, max=1.18283 at t=3.688, 状態距離最大との時刻差=0.226
- 局所占有形状差最大はt=2.154。D_corr最大t=3.200との差は1.046。

相関係数は時系列の併行性を記述するだけで、因果関係を意味しません。

## 対称性監査

- closed_nearest_rotation_residual: max 7.438e-15
- XX_closed_nearest_rotation_residual: max 4.996e-15
- closed_opposite_rotation_residual: max 5.773e-15
- XX_closed_opposite_rotation_residual: max 5.176e-15
- open_reflection_01_23_residual: max 7.161e-15
- XX_open_reflection_01_23_residual: max 6.217e-15
- open_reflection_02_13_residual: max 6.078e-15
- XX_open_reflection_02_13_residual: max 3.220e-15
- YY_closed_nearest_rotation_residual: max 6.634e-15
- YY_closed_opposite_rotation_residual: max 5.773e-15
- YY_open_reflection_01_23_residual: max 7.161e-15
- YY_open_reflection_02_13_residual: max 3.761e-15
- ZZ_closed_nearest_rotation_residual: max 7.438e-15
- ZZ_closed_opposite_rotation_residual: max 4.774e-15
- ZZ_open_reflection_01_23_residual: max 6.051e-15
- ZZ_open_reflection_02_13_residual: max 6.078e-15

## 固有空間対再構成

縮退固有空間への射影 `P_E` を使い、一体量とraw二体量を基底不変な空間対寄与へ分解しました。connected量は線形観測量ではないため、raw二体量と一体量を各々再構成した後、raw-minus-productで構成しています。

| branch | channel | total pairs | k90 | k95 | k99 |
|---|---|---:|---:|---:|---:|
| open | XX | 55 | 11 | 14 | 23 |
| open | YY | 55 | 19 | 22 | 37 |
| open | ZZ | 55 | 15 | 25 | 31 |
| closed | XX | 21 | 7 | 7 | 11 |
| closed | YY | 21 | 11 | 12 | 14 |
| closed | ZZ | 21 | 8 | 8 | 14 |

## 問いへの回答

1. 状態距離が0.5を超える区間では平均的に ZZ が最大です。ただしXXとZZは近く、単一チャネルだけの説明ではありません。
2. t=3.462では局所占有の正規化L1差が0.0403536なのに対し、D_corr=0.974638です。最大の単一成分差はXX pair 03 の0.472156です。
3. closedの最近接4対と対向2対は、下記対称性残差の範囲で回転対称性を保ちます。
4. openでは反射対応する対は一致しますが、端点を含む対と内部・離隔対は同一にはならず、その差が相関マップに現れます。
5. 各チャネル最大値はXX=0.604702, YY=0.469913, ZZ=0.665558です。最大値ではZZ、統合距離最大t=3.200ではZZが最大ですがXXも同程度です。
6. D_corr最大はt=3.200で、状態距離最大t=3.462とも局所占有差最大t=2.154とも一致しません。
7. 99%基準ではopenに23--37対、closedに11--14対を要し、特にopenは少数対だけより多モード競合として記述する方が安全です。
8. 一体+全connected特徴距離と状態距離のPearson相関は0.981143です。局所占有形状差では弱く見える時刻にも差を可視化しますが、完全な状態再構成ではありません。

## 数値監査

- exact_onebody: 7.140e-07
- exact_raw_two_body: 6.221e-07
- exact_connected: 5.473e-07
- max_imaginary_exact: 1.405e-16
- max_expectation_bound_excess: 8.882e-16
- max_hermiticity_residual: 0.000e+00
- connected_raw_minus_product_identity: 0.000e+00
- global_phase_expectation_residual: 2.225e-16
- dt_onebody: 5.355e-07
- dt_raw_two_body: 4.666e-07
- dt_connected: 4.104e-07
- mode_raw_max_residual: 6.661e-16
- mode_onebody_max_residual: 9.437e-16
- connected_from_reconstructed_ingredients_max_residual: 1.443e-15
- dt002_max_norm_error: 6.402e-13
- dt001_max_norm_error: 7.814e-13

Rust testsは10/10合格、`cargo fmt --check`合格。global phaseに依存しない期待値・connected相関・距離のみを物理比較に使用しました。

## 結論の範囲

これは有限N=4・有限時間の固定Hamiltonian数値結果です。一般的輸送則、相転移、熱力学極限、N=7 Lindblad逆流との同一機構、または円環が応答移動を因果的に起こすという主張ではありません。
