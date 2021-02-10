import pandas as pd
import solar_radiation


def add_inclined_solar_radiation_to_csv(directory_name: str, csv_file_name: str):
    """
    気象データCSVファイルを読み込み、傾斜面日射量を追加して別名のCSファイルで保存する

    :param directory_name:  ディレクトリ名
    :param csv_file_name:   CSVファイル名
    :return:
    """

    # CSVファイルを読み込む
    df = pd.read_csv(directory_name + '/' + csv_file_name, index_col=0, encoding="shift-jis")

    # 不要な列を削除
    df = df.drop("Unnamed: 10", axis=1)

    # 列名を変更（"["や"/"があるとうまくデータを扱えないため）
    df = df.rename(
        columns={'外気温[℃]': '外気温_degree', '外気絶対湿度 [kg/kgDA]': '外気絶対湿度_kg_kgDA',
                 '法線面直達日射量 [W/m2]': '法線面直達日射量_W_m2', '水平面天空日射量 [W/m2]': '水平面天空日射量_W_m2',
                 '水平面夜間放射量 [W/m2]': '水平面夜間放射量_W_m2', '太陽高度角[度]': '太陽高度角_度',
                 '太陽方位角[度]': '太陽方位角_度'})

    # 計算結果格納用配列を用意
    inclined_solar_radiations_0 = []  # 傾斜面日射量（傾斜角0°）, W/m2
    inclined_solar_radiations_30 = []  # 傾斜面日射量（傾斜角30°）, W/m2
    inclined_solar_radiations_90 = []  # 傾斜面日射量（傾斜角90°）, W/m2

    # 気象データの行ループ
    for row in df.itertuples():

        # Note: 傾斜面の方位よらない円柱面の傾斜面日射量とするため、太陽方位角と傾斜面方位角には同じ値を与える
        # 傾斜面日射量（傾斜角0°）を計算
        inclined_solar_radiation_0 = solar_radiation.get_solar_radiation_on_inclined_surfaces(
            normal_surface_direct_radiation=row.法線面直達日射量_W_m2,
            horizontal_surface_sky_radiation=row.水平面天空日射量_W_m2,
            solar_altitude=row.太陽高度角_度, solar_azimuth=row.太陽方位角_度,
            surface_tilt_angle=0.0, surface_azimuth=row.太陽方位角_度
        )

        # 傾斜面日射量（傾斜角30°）を計算
        inclined_solar_radiation_30 = solar_radiation.get_solar_radiation_on_inclined_surfaces(
            normal_surface_direct_radiation=row.法線面直達日射量_W_m2,
            horizontal_surface_sky_radiation=row.水平面天空日射量_W_m2,
            solar_altitude=row.太陽高度角_度, solar_azimuth=row.太陽方位角_度,
            surface_tilt_angle=30.0, surface_azimuth=row.太陽方位角_度
        )

        # 傾斜面日射量（傾斜角90°）を計算
        inclined_solar_radiation_90 = solar_radiation.get_solar_radiation_on_inclined_surfaces(
            normal_surface_direct_radiation=row.法線面直達日射量_W_m2,
            horizontal_surface_sky_radiation=row.水平面天空日射量_W_m2,
            solar_altitude=row.太陽高度角_度, solar_azimuth=row.太陽方位角_度,
            surface_tilt_angle=90.0, surface_azimuth=row.太陽方位角_度
        )

        # 計算結果を配列に格納
        inclined_solar_radiations_0.append(inclined_solar_radiation_0)
        inclined_solar_radiations_30.append(inclined_solar_radiation_30)
        inclined_solar_radiations_90.append(inclined_solar_radiation_90)

    # 計算結果をDataFrameに追加
    df['傾斜面日射量_0度_W_m2'] = inclined_solar_radiations_0
    df['傾斜面日射量_30度_W_m2'] = inclined_solar_radiations_30
    df['傾斜面日射量_90度_W_m2'] = inclined_solar_radiations_90

    # CSVファイル出力
    df.to_csv(directory_name + '/rev_' + csv_file_name, encoding="shift-jis")


def edit_all_climate_data():
    """
    気象データファイルに傾斜面日射量を追加する処理を行う
    :return:
    """

    # 気象データファイル名のリストを設定する
    directory_name = 'climateData'
    csv_file_name_list = ['climateData_1.csv', 'climateData_2.csv', 'climateData_3.csv', 'climateData_4.csv',
                          'climateData_5.csv', 'climateData_6.csv', 'climateData_7.csv', 'climateData_8.csv']

    # 傾斜面日射量を追加し、別名のCSVファイルとして保存する
    for file_name in csv_file_name_list:
        add_inclined_solar_radiation_to_csv(directory_name=directory_name, csv_file_name=file_name)


if __name__ == '__main__':

    edit_all_climate_data()
