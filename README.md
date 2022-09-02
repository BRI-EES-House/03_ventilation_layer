# 03_ventilation_layer
外壁通気層（屋根通気層・通気ブロック・瓦を含む）の効果を検討するリポジトリ

## Pythonコード

### global_number.py
	物性値を定義しているファイル。

### ventilation_wall_parameters.py
	総当たりのパラメータと計算結果を取得し、CSVに出力する処理を行うファイル。
	詳細計算、簡易計算No.1～4、放射熱伝達率、対流熱伝達率の検証に対応。

### ventilation_wall.py
	詳細計算（熱収支式を解き、通気層の状態値を取得する）を行う関数を定義しているファイル。
	戻り値はdataclass（WallStatusValues）で定義。

### ventilation_wall_simplified.py
	簡易計算No.1～4を行う関数を定義しているファイル。

### envelope_performance_factors.py
	通気層を有する壁体の熱貫流率や、表面熱流などを計算する関数を定義しているファイル。
	最終的には仕様していない関数があると思われる。。

### heat_transfer_coefficient.py
	放射熱伝達率、対流熱伝達率を計算する関数を定義しているファイル。
	それぞれ、詳細計算と簡易計算（回帰式）の両方を定義している。

### solar_radiation.py
	傾斜面日射量を求める関数を定義しているファイル。

### climate_data_editor.py
	 気象データCSVファイルを読み込み、傾斜面日射量を追加して別名のCSファイルで保存する関数を定義しているファイル。

### boundary_condition_creator.py
	境界条件作成用に、地域区分別の冬期、夏期の平均外気温度、平均傾斜面日射量を計算する関数、通気層内の面1、面2の表面温度を計算する関数を定義しているファイル。

### validation.py
	通気層を有する壁体の熱貫流率の検証を行うための関数を定義しているファイル。

### validation0.py
	通気層を有する壁体の熱貫流率について、通気風量0での検証を行うための関数を定義しているファイル。

### reference_natural_convection.py 
	対流熱伝達率の計算結果比較のため、既往文献における自然対流熱伝達率を計算する関数を定義しているファイル。

### validation_iterate_error.py
	熱収支式の計算結果を検証するためのファイル？

### validation_ke_max_case.py
	熱収支式の計算結果を検証するためのファイル？


## Jpyter Notebook

### valid_radiative_heat_transfer.ipynb
	放射熱伝達率の計算結果を検証するためのノートブック。

### valid_natural_convective_heat_transfer.ipynb
	風速0の場合の対流熱伝達率の計算結果を検証するためのノートブック。

### valid_convective_heat_transfer.ipynb
	風速が0より大きい場合の対流熱伝達率の計算結果を検証するためのノートブック。

### ventilation_wall_analysis_comparison_simplified_calculation.ipynb
	詳細計算法と簡易計算（No.1～4）の計算結果を比較するためのノートブック。

### ventilation_wall_analysis_convection_heat_transfer_coefficient .ipynb
	通気層の対流熱伝達率の推定方法を検証するためのノートブック。

### ventilation_wall_analysis_radiative_heat_transfer_coefficient .ipynb
	通気層の放射熱伝達率の推定方法を検証するためのノートブック。

### ventilation_wall_casestudy.ipynb
	簡易計算法No.3による修正U値、修正η値の検証を行うためのノートブック。

### ventilation_wall_boundary_condition.ipynb
	簡易計算法における境界条件を検討するためのノートブック。

### ventilation_wall_analysis_tempelature_and_correction_coefficient.ipynb
	通気層の等価温度、通気層を考慮した壁体の相当熱貫流率補正係数の検証するためのノートブック。
	最終的には使用していない。