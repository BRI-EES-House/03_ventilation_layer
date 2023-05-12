import itertools
import pandas as pd
import numpy as np
import ventilation_layer.global_number as global_number
import ventilation_wall as vw
from ventilation_layer import heat_transfer_coefficient as htc


def get_average_climate(csv_file_path: str, target_month: int, target_solar_radiation_name: str) -> tuple:
    """
    日照時間帯の平均外気温度、平均傾斜面日射量を計算する関数

    :param csv_file_path:       CSVファイルへのパス
    :param target_month:        計算の対象月
    :param target_solar_radiation_name: 対象とする傾斜面日射量の項目名
    :return: 平均外気温度, degree, 平均傾斜面日射量, w/m2
    """
    # CSVファイルを読み込む
    df = pd.read_csv(csv_file_path, index_col=0, encoding="shift-jis")

    # 対象月の日射量が0より大きいデータのみを抽出
    df_target = df.query("月 == " + str(target_month) + " & " + target_solar_radiation_name + " > 0")

    return df_target['外気温_degree'].mean(), df_target[target_solar_radiation_name].mean()


def get_all_region_climate_conditions(angle_list: list) -> pd.DataFrame:
    """
    地域区分別の冬期、夏期の平均外気温度、平均傾斜面日射量を計算する関数

    :param angle_list:  傾斜角リスト, degree
    :return: 地域×季節×傾斜角の総当たりの計算結果（DataFrame）
    """

    # 地域区分別の気象データファイル名のリストを作成
    directory_name = 'climateData'
    csv_file_name = 'rev_climateData_'

    regions = [1, 2, 3, 4, 5, 6, 7, 8]
    seasons = ['winter', 'summer']
    parameter_list = list(itertools.product(regions, seasons, angle_list))
    df = pd.DataFrame(parameter_list, columns=['region', 'season', 'angle'])

    # 計算結果格納用の配列を用意
    theta_e_ave = []
    j_surf_ave = []

    for row in df.itertuples():

        # 気象データファイルへのパスを設定
        csv_file_path = directory_name + '/' + csv_file_name + str(row.region) + '.csv'

        # 対象とする傾斜面日射量の項目名を設定
        target_solar_radiation_name = '傾斜面日射量_' + str(row.angle) + '度_W_m2'

        # 季節により対象となる月を設定
        target_month = 1
        if row.season == 'winter':
            target_month = 1
        elif row.season == 'summer':
            target_month = 8

        # 指定した月の日射のある時間帯の平均外気温、平均傾斜面日射量を求める
        buf_theta_e_ave, buf_j_surf_ave = get_average_climate(csv_file_path=csv_file_path,
                                                              target_month=target_month,
                                                              target_solar_radiation_name=target_solar_radiation_name)

        # 平均外気温、平均傾斜面日射量を配列に格納
        theta_e_ave.append(buf_theta_e_ave)
        j_surf_ave.append(buf_j_surf_ave)

    # 計算結果をDataFrameに追加
    df['theta_e_ave'] = theta_e_ave
    df['j_surf_ave'] = j_surf_ave

    return df


def calc_ventilation_wall_surface_temperatures(angle: float, theta_e: float, j_surf: float, season: str) -> tuple:
    """
    通気層内の面1、面2の表面温度を計算する関数

    :param angle:       通気層の傾斜角, degree
    :param theta_e:     外気温度, degree
    :param j_surf:      日射量, W/m2
    :param season:      季節（'winter' or 'summer'）
    Returns:
        通気層内の面1、面2の表面温度, degrees
    """

    # 固定値の設定
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    # 対流熱伝達、放射熱伝達の計算方法を指定
    calc_mode_h_cv = 'detailed'
    calc_mode_h_rv = 'detailed'

    # the room temperature, degrees
    theta_r = 20.0 if season == 'winter' else 27.0

    # the solar absorption ratio on the exterior surface, -
    a_surf = np.array([0.0, np.median([0.0, 1.0]), 1.0], dtype=float).mean()

    # the thermal conductance of the outside material, W/m2K
    c_1 = np.array([0.5, np.median([0.5, 100.0]), 100.0], dtype=float).mean()

    # the thermal conductance of the inside material, W/m2K
    c_2 = np.array([0.1, np.median([0.1, 5.0]), 5.0], dtype=float).mean()

    # the length of the ventilation layer, m
    l_h = np.array([3.0, np.median([3.0, 12.0]), 12.0], dtype=float).mean()

    # the width of the ventilation layer, m
    l_w = np.array([0.05, np.median([0.05, 10.0]), 10.0], dtype=float).mean()

    # the thickness of the ventilation layer, m
    l_d = np.array([0.05, np.median([0.05, 0.3]), 0.3], dtype=float).mean()

    # the mean air velocity of the ventilation layer, m/s
    v_a = np.array([0.0, np.median([0.0, 1.0]), 1.0], dtype=float).mean()

    # the emissivity of the surface 2 facing the ventilation layer, -
    emissivity_2=np.array([0.1, np.median([0.1, 0.9]), 0.9], dtype=float).mean()

    # 通気層の状態値を取得
    theta_out_surf, theta_1_surf, theta_2_surf, theta_in_surf, theta_as, _ = vw.get_wall_status_values(
        index=0,
        theta_e=theta_e,
        theta_r=theta_r,
        j_surf=j_surf,
        a_surf=a_surf,
        c_1=c_1,
        c_2=c_2,
        l_h=l_h,
        l_w=l_w,
        l_d=l_d,
        angle=angle,
        v_a=v_a,
        eps_1=0.9,
        eps_2=emissivity_2,
        calc_mode_h_cv=calc_mode_h_cv, calc_mode_h_rv=calc_mode_h_rv, h_out=h_out, h_in=h_in
    )

    # the convective thermal transmittance coefficient, W/m2K
    h_cv = htc.get_h_cv(calc_mode=calc_mode_h_cv, v_a=v_a, theta_1=theta_1_surf, theta_2=theta_2_surf, angle=angle, l_h=l_h, l_d=l_d)

    # the heat flow from the outside(include solar irradiance) to the exterior surface, W/m2
    q_outer_flow = vw.get_q_flow_out(theta_e=theta_e, a_surf=a_surf, j_surf=j_surf, theta_0=theta_out_surf, h_out=h_out)

    # the exhausted heat flow from the ventilation layer, W/m2
    q_exhaust_flow = vw.get_q_flow_exhaust(v_a=v_a, l_d=l_d, l_w=l_w, l_h=l_h, theta_1=theta_1_surf, theta_2=theta_2_surf, theta_4=theta_as, theta_as_in=theta_e, h_cv=h_cv)

    # the heat flow to the inside, W/m2
    q_inner_flow = vw.get_q_flow_in(theta_3=theta_in_surf, theta_r=theta_r, h_in=h_in)

    return theta_1_surf, theta_2_surf, theta_as, q_outer_flow, q_exhaust_flow, q_inner_flow


def add_ventilation_wall_temperatures_and_heat_flow(target_df: pd.DataFrame) -> pd.DataFrame:
    """
    通気層内の面1、面2の表面温度を計算し、DataFrameに追加して返す関数

    :param target_df: 平均外気温度、平均日射量のDataFrame
    :return: 通気層内の面1、面2の表面温度を追加したDataFrame
    """

    theta_1_surf = []
    theta_2_surf = []
    theta_as_ave = []
    q_outer_flow = []
    q_exhaust_flow = []
    q_inner_flow = []

    for row in target_df.itertuples():
        _theta_1_surf, _theta_2_surf, _theta_as_ave, buf_q_outer_flow, buf_q_exhaust_flow, buf_q_inner_flow = calc_ventilation_wall_surface_temperatures(
            angle=row.angle, theta_e=row.theta_e_ave, j_surf=row.j_surf_ave, season=row.season
        )
        theta_1_surf.append(_theta_1_surf)
        theta_2_surf.append(_theta_2_surf)
        theta_as_ave.append(_theta_as_ave)
        q_outer_flow.append(buf_q_outer_flow)
        q_exhaust_flow.append(buf_q_exhaust_flow)
        q_inner_flow.append(buf_q_inner_flow)

    target_df['theta_1_surf'] = theta_1_surf
    target_df['theta_2_surf'] = theta_2_surf
    target_df['theta_as_ave'] = theta_as_ave
    target_df['q_outer_flow'] = q_outer_flow
    target_df['q_exhaust_flow'] = q_exhaust_flow
    target_df['q_inner_flow'] = q_inner_flow

    return target_df
