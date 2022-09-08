# 03_ventilation_layer
外壁通気層（屋根通気層・通気ブロック・瓦を含む）の効果を検討するリポジトリ

## Pythonコード

### global_number.py
- 物性値や固定値を定義しているファイル。

### ventilation_wall_parameters.py
- 総当たりのパラメータと計算結果を取得し、CSVに出力する処理を行うファイル。
- 詳細計算、簡易計算No.1～4、放射熱伝達率、対流熱伝達率の検証に対応。
- 関数dump_csv_all_case_resultを実行すると、全ケースの計算結果をCSVファイルとして出力する。ただし処理に時間がかかるので、不要な処理はコメントアウトする。

### ventilation_wall.py
- 詳細計算（熱収支式を解き、通気層の状態値を取得する）を行う関数を定義しているファイル。
- 戻り値はdataclass（WallStatusValues）で定義。

### ventilation_wall_simplified.py
- 簡易計算No.1～4を行う関数を定義しているファイル。

### envelope_performance_factors.py
- 通気層を有する壁体の熱貫流率や、表面熱流などを計算する関数を定義しているファイル。
- 最終的には使用していない関数があると思われる。
- 当初は通気層を考慮した壁体の相当熱貫流率、日射熱取得率を直接計算する方法を検討していたが、通気層内の放射熱伝達率、対流熱伝達率を固定化して、室内表面熱流を直接計算する方法に変更したため、下記の関数は使用していない。
  - overall_heat_transfer_coefficient
  - get_k_e
  - solar_heat_gain_coefficient

### heat_transfer_coefficient.py
- 放射熱伝達率、対流熱伝達率を計算する関数を定義しているファイル。
- それぞれ、詳細計算と簡易計算（回帰式）の両方を定義している。

### reference_natural_convection.py 
- 対流熱伝達率の計算結果比較のため、既往文献における自然対流熱伝達率を計算する関数を定義しているファイル。

### solar_radiation.py
- 傾斜面日射量を求める関数を定義しているファイル。

### climate_data_editor.py
-  気象データCSVファイルを読み込み、傾斜面日射量を追加して別名のCSファイルで保存する関数を定義しているファイル。

### boundary_condition_creator.py
- 境界条件作成用に、地域区分別の冬期、夏期の平均外気温度、平均傾斜面日射量を計算する関数、通気層内の面1、面2の表面温度を計算する関数を定義しているファイル。

### validation.py
- 通気層を有する壁体の熱貫流率の検証を行うための関数を定義しているファイル。
- 通気層を考慮した壁体の相当熱貫流率、日射熱取得率の計算結果の検証用のため、最終的には使用していない。

### validation0.py
- 通気層を有する壁体の熱貫流率について、通気風量0での検証を行うための関数を定義しているファイル。
- 通気層を考慮した壁体の相当熱貫流率、日射熱取得率の計算結果の検証用のため、最終的には使用していない。

### validation_iterate_error.py
- 熱収支式の計算結果を検証するためのファイル。最終的には使用していない。

### validation_ke_max_case.py
- 熱収支式の計算結果を検証するためのファイル。最終的には使用していない。


## Jpyter Notebook

### valid_radiative_heat_transfer.ipynb
- 放射熱伝達率の計算結果を検証するためのノートブック。
- 通気層が無限の平行面と仮定した場合と、通気胴縁や屋根垂木等の影響を考慮した有限の二次元空間と仮定した場合の計算結果の違いを検証している。
- 報告書の図2.3、図2.4の描画用。
- heat_transfer_coefficient.pyをインポートして使用。

### valid_natural_convective_heat_transfer.ipynb
- 風速0の場合の自然対流熱伝達率の計算結果を検証するためのノートブック。
- 本計算方法による自然対流熱伝達率と、既往文献の計算方法による自然対流熱伝達率の計算結果を比較している。
- 報告書の図2.5の描画用。
- heat_transfer_coefficient.py、reference_natural_convection.py をインポートして使用。

### valid_convective_heat_transfer.ipynb
- 有風状態の場合の対流熱伝達率の計算結果を検証するためのノートブック。
- 本計算方法による計算方法と、既往文献の計算方法による対流熱伝達率の計算結果を比較している。
- 報告書の図2.6の描画用。
- heat_transfer_coefficient.pyをインポートして使用。

### ventilation_wall_analysis_convection_heat_transfer_coefficient .ipynb
- 通気層の対流熱伝達率の推定方法を検証するためのノートブック。
- ventilation_wall_parameters.pyで出力した詳細計算、対流熱伝達率の簡易計算によるパラメトリックスタディ結果のCSVファイルを読み込み、グラフを描画する。
- 報告書の図4.1～図4.8の描画用。
- グラフのサイズは、関数setPltSingle内で設定する。 

### ventilation_wall_analysis_radiative_heat_transfer_coefficient .ipynb
- 通気層の放射熱伝達率の推定方法を検証するためのノートブック。
- ventilation_wall_parameters.pyで出力した詳細計算、放射熱伝達率の簡易計算によるパラメトリックスタディ結果のCSVファイルを読み込み、グラフを描画する。
- 報告書の図4.9～図4.14の描画用。
- グラフのサイズは、関数setPltSingle内で設定する。 

### ventilation_wall_analysis_comparison_simplified_calculation.ipynb
- 詳細計算法と簡易計算No.1～4の計算結果を比較するためのノートブック。
- ventilation_wall_parameters.pyで出力した詳細計算、簡易計算No.1～4の各種パラメトリックスタディ結果のCSVファイルを読み込み、グラフを描画する。
- 報告書の図5.6の描画用。 

### ventilation_wall_casestudy.ipynb
- 簡易計算法No.3による修正U値、修正η値の検証を行うためのノートブック。
- 2021年の空気調和・衛生工学会での発表論文（通気層を有する壁体の熱的性能の評価方法に関する研究（第二報）通気層を有する壁体の簡易計算法の提案）のFig.8～Fig.10を描画。
- 下記のPythonファイルをインポートして使用（envelope_performance_factors.pyはインポートしているが、使用していない）。
   - ventilation_wall.py
   - ventilation_wall_simplified.py
   - global_number.py

### ventilation_wall_boundary_condition.ipynb
- 簡易計算法における境界条件を検討するためのノートブック。
- climate_data_editor.pyで出力した地域区分別の気象データのCSVファイルを読み込み、計算に使用。
- 報告書の表7～表15の算出用。
- 下記のPythonファイルをインポートして使用（ventilation_wall.pyはインポートしているが、使用していない）。
   - boundary_condition_creator.py
   - heat_transfer_coefficient.py 

### ventilation_wall_analysis_tempelature_and_correction_coefficient.ipynb
- 通気層の等価温度の簡易計算法を検証するためのノートブック。最終的には使用していない。
- 当初は通気層を考慮した壁体の相当熱貫流率（U_e）を直接計算する方法を検討していた。相当熱貫流率の計算には、通気層内の等価温度（θ_as_e）が計算できればよいので、通気層内の等価温度簡易に計算する方法を検討している。
- ただし、通気層内の放射熱伝達率、対流熱伝達率を固定化して、室内表面熱流を直接計算する方法に変更したため、この評価は使用していない。

### ventilation_wall_test.ipynb
- 通気層を有する壁体に関する計算結果を確認するためのノートブック。最終的には使用していない。