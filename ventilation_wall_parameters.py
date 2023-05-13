import itertools
import pandas as pd
import numpy as np
from typing import List, Tuple

import ventilation_layer.global_number as global_number
import ventilation_layer.ventilation_wall as vw
import ventilation_wall_simplified as vws
import envelope_performance_factors as epf
import ventilation_layer.heat_transfer_coefficient as htc


class Log:
    def write(self, msg):
        """
        np.seterrのエラーログ記録用の処理（エラー発生時はコンソールにメッセージを出力する）
        :param msg:　エラーメッセージ
        :return: なし
        """
        print("LOG: %s" % msg)


def get_parameter_data_frame() -> pd.DataFrame:
    """
    複数のパラメータの総当たりの組み合わせ（直積）のリストを作成する
    :param なし
    :return: 総当たりのパラメータリスト
    """


    # the outdoor temperature, degrees
    # 外気温度は、冬期条件（-10.0～10.0degC）、夏期条件（25.0～35.0degC）を]それぞれ与える
    # [-10.0, 0.0, 10.0, 25.0, 30.0, 35.0]
    theta_e = [-10.0, 0.0, 10.0, 25.0, 30.0, 35.0]

    # 上記以外のパラメータには、一部を除いて想定される上下限値と中央値の3点を与える

    # the solar irradiance, W/m2
    # [0.0, 500.0, 1000.0]
    j_surf = [0.0, 500.0, 1000.0]

    # the solar absorption ratio, -
    # [0.0, 0.5, 1.0]
    a_surf = [0.0, 0.5, 1.0]

    # the thermal conductance of the exterior material, W/m2K
    # [0.5, 50.25, 100.0]
    C_1 = [0.5, 50.25, 100.0]
    
    # the thermal conductance of the interior material, W/m2K
    # [0.1, 2.55, 5.0]
    C_2 = [0.1, 2.55, 5.0]
    
    # the length of the ventilation layer, m
    # [3.0, 7.5, 12.0]
    l_h = [3.0, 7.5, 12.0]
    
    # the width of the ventilation layer, m
    # [0.05, 5.025, 10.0]
    l_w = [0.05, 5.025, 10.0]

    # the thickness of the ventilation layer, m
    # [0.05, 0.175, 0.3]
    l_d = [0.05, 0.175, 0.3]
    
    # the angle of the ventilation layer, degrees
    # [0.0, 45.0, 90.0]
    angle = [0.0, 45.0, 90.0]

    # the air vlocity of the ventilation layer, m/s
    # [0.0, 0.5, 1.0]
    v_a = [0.0, 0.5, 1.0]
    
    # the interval between the furring or the rafter, m
    # [0.45]
    l_s = [0.45]
    
    # the emissivity of the surface 1 facing the ventilation layer, -
    # [0.9]
    emissivity_1 = [0.9]

    # the emissivity of the surface 1 facing the ventilation layer, -
    # [0.1, 0.5, 0.9]
    emissivity_2 = [0.1, 0.5, 0.9]
    
    parameter_list = list(itertools.product(theta_e, j_surf, a_surf, C_1, C_2, l_h, l_w, l_d, angle, v_a, l_s, emissivity_1, emissivity_2))
    
    parameter_name = ['theta_e', 'j_surf', 'a_surf', 'C_1', 'C_2', 'l_h', 'l_w', 'l_d', 'angle', 'v_a', 'l_s', 'emissivity_1', 'emissivity_2']

    df = pd.DataFrame(parameter_list, columns=parameter_name)

    # Give the temperature of 20.0 degrees for winter and 27.0 degrees for summer as the indoor temperature.
    df['theta_r'] = np.where(df.theta_e > 20.0, 27.0, 20.0)

    return df


def get_wall_status_data_by_detailed_calculation(calc_mode_h_cv: str, calc_mode_h_rv: str) -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、各ケースの計算結果を保有するDataFrameを作成する

    :param calc_mode_h_cv: 対流熱伝達率の計算モード
    :param calc_mode_h_rv: 放射熱伝達率の計算モード
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    df = get_parameter_data_frame()

    # 固定値の設定
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
            
        # 通気層の状態値を取得
        results = [
            vw.get_wall_status_values(
                index=row[0],
                theta_e=row.theta_e,
                theta_r=row.theta_r,
                j_surf=row.j_surf,
                a_surf=row.a_surf,
                c_1=row.C_1,
                c_2=row.C_2,
                l_h=row.l_h,
                l_w=row.l_w,
                l_d=row.l_d,
                angle=row.angle,
                v_a=row.v_a,
                eps_1=row.emissivity_1,
                eps_2=row.emissivity_2,
                calc_mode_h_cv=calc_mode_h_cv,
                calc_mode_h_rv=calc_mode_h_rv,
                h_out=h_out,
                h_in=h_in
            ) for row in df.itertuples()
        ]

        # the temperature of the surface of exterior, degrees
        theta_out_surf = np.array([result[0] for result in results])

        # the temperature of the surface of exterior side facing the ventilation layer, degrees
        theta_1 = np.array([result[1] for result in results])

        # the temperature of the surface of interior side facing the ventilation layer, degrees
        theta_2 = np.array([result[2] for result in results])

        # the temperature on the surface of the inside, degrees
        theta_in_surf = np.array([result[3] for result in results])    

        # the air tempearature in the ventilation layer, degrees
        theta_as_ave = np.array([result[4] for result in results])

        # the heat balance of the surface of the exterilr, W/m2
        heat_balance_0 = [result[5][0] for result in results]

        # the heat balance of the surface of the exterior side of the ventilation layer, W/m2
        heat_balance_1 = [result[5][1] for result in results]

        # the heat balance of the surface of the interior side of the ventilation layer, W/m2
        heat_balance_2 = [result[5][2] for result in results]

        # the heat balance of the surface of the inside, W/m2
        heat_balance_3 = [result[5][3] for result in results]

        # the heat balance of the air in the ventilation layer, W/m2
        heat_balance_4 = [result[5][4] for result in results]

        is_optimize_succeed = [result[6].is_optimize_succeed for result in results]
        optimize_message = [result[6].optimize_message for result in results]

        c_2 = df.C_2.to_numpy()
        v_a = df.v_a.to_numpy()
        angle =df.angle.to_numpy()
        l_h = df.l_h.to_numpy()
        l_d = df.l_d.to_numpy()
        theta_r = df.theta_r.to_numpy()
        theta_e = df.theta_e.to_numpy()
        a_surf = df.a_surf.to_numpy()
        j_surf = df.j_surf.to_numpy()
        eps1 = df.emissivity_1.to_numpy()
        eps2 = df.emissivity_2.to_numpy()

        # the effective emissivity, -
        eps_eff = np.vectorize(htc.get_e)(eps1=eps1, eps2=eps2)

        # the convective heat transfer coefficient, W/m2K
        h_cv = np.vectorize(htc.get_h_cv)(calc_mode=calc_mode_h_cv, v_a=v_a, theta_1=theta_1, theta_2=theta_2, angle=angle, l_h=l_h, l_d=l_d)

        # the radiative heat transfer coefficient, W/m2K
        h_rv = np.vectorize(htc.get_h_rv)(eps_eff=eps_eff, calc_mode=calc_mode_h_rv, theta_1=theta_1, theta_2=theta_2)

        # the equivallent temperature of the ventilation layer, degrees
        theta_as_e = np.vectorize(epf.get_theata_as_e)(theta_as_ave=theta_as_ave, theta_1_surf=theta_1, h_cv=h_cv, h_rv=h_rv)

        # SAT temperature, degrees
        theta_sat = epf.get_theta_SAT(theta_e=theta_e, a_surf=a_surf, j_surf=j_surf, h_out=h_out)

        # the thermal resistance of the interior material, m2K/W
        r_i = epf.get_r_i(C_2=c_2)

        # the heat flow of the inside surface, W/m2
        q_room_side = epf.get_heat_flow_room_side_by_vent_layer_heat_resistance(r_i=r_i, theta_2=theta_2, theta_r=theta_r)

        # the correction factor for calculating the equivalent thermal transmission coefficient of the wall with the ventilation layer
        k_e = np.vectorize(epf.get_k_e)(theta_as_e=theta_as_e, theta_r=theta_r, theta_SAT=theta_sat)

        # 計算結果をDataFrameに追加
        df['theta_sat'] = theta_sat
        df['theta_out_surf'] = theta_out_surf
        df['theta_1_surf'] = theta_1
        df['theta_2_surf'] = theta_2
        df['theta_in_surf'] = theta_in_surf
        df['theta_as_ave'] = theta_as_ave
        df['effective_emissivity'] = eps_eff
        df['h_cv'] = h_cv
        df['h_rv'] = h_rv
        df['theta_as_e'] = theta_as_e
        df['k_e'] = k_e
        df['q_room_side'] = q_room_side
        df['heat_balance_0'] = heat_balance_0
        df['heat_balance_1'] = heat_balance_1
        df['heat_balance_2'] = heat_balance_2
        df['heat_balance_3'] = heat_balance_3
        df['heat_balance_4'] = heat_balance_4
        df['is_optimize_succeed'] = is_optimize_succeed
        df['optimize_message'] = optimize_message

    return df


def get_wall_status_data_by_simplified_calculation_no_01() -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、簡易計算法案No.1（簡易版の行列式）による計算結果を保有するDataFrameを作成する

    :param: なし
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    df = get_parameter_data_frame()

    # 固定値の設定
    h_out = global_number.get_h_out()

    # 計算結果格納用配列を用意
    theta_sat = []          # 相当外気温度[℃]
    theta_1_surf = []       # 通気層に面する面1の表面温度[℃]
    theta_2_surf = []       # 通気層に面する面1の表面温度[℃]]
    theta_as_ave = []       # 通気層の平均温度[℃]
    effective_emissivity = []    # 有効放射率[-]
    h_cv = []               # 通気層の対流熱伝達率[W/(m2・K)]
    h_rv = []               # 通気層の放射熱伝達率[W/(m2・K)]
    q_room_side = []        # 室内表面熱流[W/m2]

    # エラーログ出力用の設定
    log = Log()
    saved_handler = np.seterrcall(log)

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
        for row in df.itertuples():
            print(row[0])
            # パラメータを設定
            parms = (vws.Parameters(theta_e=row.theta_e,
                                   theta_r=row.theta_r,
                                   J_surf=row.j_surf,
                                   a_surf=row.a_surf,
                                   C_1=row.C_1,
                                   C_2=row.C_2,
                                   l_h=row.l_h,
                                   l_w=row.l_w,
                                   l_d=row.l_d,
                                   angle=row.angle,
                                   v_a=row.v_a,
                                   l_s=row.l_s,
                                   emissivity_1=row.emissivity_1,
                                   emissivity_2=row.emissivity_2))

            # 通気層の状態値を取得
            temps, h_cv_buf, h_rv_buf, r_i_buf = vws.get_vent_wall_temperature_by_simplified_calculation_no_01(parm=parms, h_out=h_out)
            theta_1_surf.append(temps[0])
            theta_2_surf.append(temps[2])
            theta_as_ave.append(temps[1])
            effective_emissivity.append(htc.get_e(eps1=row.emissivity_1, eps2=row.emissivity_2))
            h_cv.append(h_cv_buf)
            h_rv.append(h_rv_buf)

            # 相当外気温度を計算
            theta_sat_buf = epf.get_theta_SAT(theta_e=row.theta_e, a_surf=row.a_surf, j_surf=row.j_surf, h_out=h_out)
            theta_sat.append(theta_sat_buf)

            # 室内側表面熱流を計算
            q_room_side.append(epf.get_heat_flow_room_side_by_vent_layer_heat_resistance(r_i=r_i_buf, theta_2=temps[2], theta_r=row.theta_r))

    # 計算結果をDataFrameに追加
    df['theta_sat'] = theta_sat
    df['theta_1_surf'] = theta_1_surf
    df['theta_2_surf'] = theta_2_surf
    df['theta_as_ave'] = theta_as_ave
    df['effective_emissivity'] = effective_emissivity
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['q_room_side'] = q_room_side

    return df


def get_wall_status_data_by_simplified_calculation_no_02() -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、簡易計算法案No.2（簡易式）による計算結果を保有するDataFrameを作成する

    :param: なし
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    df = get_parameter_data_frame()
                      
    # 固定値の設定
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    # 計算結果格納用配列を用意
    theta_sat = []          # 相当外気温度[℃]
    theta_as_ave = []       # 通気層の平均温度[℃]
    effective_emissivity = []    # 有効放射率[-]
    h_cv = []               # 通気層の対流熱伝達率[W/(m2・K)]
    h_rv = []               # 通気層の放射熱伝達率[W/(m2・K)]
    u_o = []                # 室外側から通気層までの熱貫流率[W/(m2・K)]
    u_i = []                # 室内側から通気層までの熱貫流率[W/(m2・K)]
    q_room_side = []        # 室内表面熱流[W/m2]

    # エラーログ出力用の設定
    log = Log()
    saved_handler = np.seterrcall(log)

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
        for row in df.itertuples():
            print(row[0])
            # パラメータを設定
            parms = (vws.Parameters(theta_e=row.theta_e,
                                   theta_r=row.theta_r,
                                   J_surf=row.j_surf,
                                   a_surf=row.a_surf,
                                   C_1=row.C_1,
                                   C_2=row.C_2,
                                   l_h=row.l_h,
                                   l_w=row.l_w,
                                   l_d=row.l_d,
                                   angle=row.angle,
                                   v_a=row.v_a,
                                   l_s=row.l_s,
                                   emissivity_1=row.emissivity_1,
                                   emissivity_2=row.emissivity_2))

            # 対流熱伝達率、放射熱伝達率を計算
            effective_emissivity_buf = htc.get_e(eps1=row.emissivity_1, eps2=row.emissivity_2)
            if parms.theta_r == 20.0:
                h_cv_buf = htc.get_h_cv(calc_mode="simplified_winter", v_a=row.v_a)
                h_rv_buf = htc.get_h_rv(calc_mode="simplified_winter", eps_eff=effective_emissivity_buf)
            else:
                h_cv_buf = htc.get_h_cv(calc_mode="simplified_summer", v_a=row.v_a)
                h_rv_buf = htc.get_h_rv(calc_mode="simplified_summer", eps_eff=effective_emissivity_buf)

            effective_emissivity.append(effective_emissivity_buf)
            h_cv.append(h_cv_buf)
            h_rv.append(h_rv_buf)

            # 通気層平均温度を取得
            theta_as_ave_buf, u_o_buf, u_i_buf = vws.get_vent_wall_temperature_by_simplified_calculation_no_02(parm=parms, h_out=h_out)
            theta_as_ave.append(theta_as_ave_buf)

            # 相当外気温度を計算
            theta_sat.append(epf.get_theta_SAT(row.theta_e, row.a_surf, row.j_surf, h_out))

            # 室外側から通気層までの熱貫流率、室内側から通気層までの熱貫流率
            u_o.append(u_o_buf)
            u_i.append(u_i_buf)

            # 室内側表面熱流を計算
            q_room_side.append(epf.get_heat_flow_room_side_by_vent_layer_heat_transfer_coeff(u_i=u_i_buf, theta_as_ave=theta_as_ave_buf, theta_r=row.theta_r))

    # 計算結果をDataFrameに追加
    df['theta_sat'] = theta_sat
    df['theta_as_ave'] = theta_as_ave
    df['effective_emissivity'] = effective_emissivity
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['u_o'] = u_o
    df['u_i'] = u_i
    df['q_room_side'] = q_room_side

    return df


def get_wall_status_data_by_simplified_calculation_no_03() -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、簡易計算法案No.3（通気層を有する壁体の修正熱貫流率、修正日射熱取得率から
    室内表面熱流を求める）による計算結果を保有するDataFrameを作成する

    :param: なし
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    parameter_name = ['theta_e', 'theta_r', 'j_surf', 'a_surf', 'C_1', 'C_2', 'l_h', 'l_w', 'l_d', 'angle',
                      'v_a', 'l_s', 'emissivity_1', 'emissivity_2']
    df = get_parameter_data_frame()
                      
    # 固定値の設定
    h_out = global_number.get_h_out()

    # 計算結果格納用配列を用意
    theta_sat = []          # 相当外気温度[℃]
    h_cv = []               # 通気層の対流熱伝達率[W/(m2・K)]
    h_rv = []               # 通気層の放射熱伝達率[W/(m2・K)]
    u_dash = []             # 修正熱貫流率[W/(m2・K)]
    eta_dash = []           # 修正日射熱取得率[-]
    q_room_side = []        # 室内表面熱流[W/m2]

    # エラーログ出力用の設定
    log = Log()
    saved_handler = np.seterrcall(log)

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
        for row in df.itertuples():
            print(row[0])
            # パラメータを設定
            parms = (vws.Parameters(theta_e=row.theta_e,
                                   theta_r=row.theta_r,
                                   J_surf=row.j_surf,
                                   a_surf=row.a_surf,
                                   C_1=row.C_1,
                                   C_2=row.C_2,
                                   l_h=row.l_h,
                                   l_w=row.l_w,
                                   l_d=row.l_d,
                                   angle=row.angle,
                                   v_a=row.v_a,
                                   l_s=row.l_s,
                                   emissivity_1=row.emissivity_1,
                                   emissivity_2=row.emissivity_2))

            # 相当外気温度を計算
            theta_sat.append(epf.get_theta_SAT(row.theta_e, row.a_surf, row.j_surf, h_out))

            # 対流熱伝達率、放射熱伝達率、修正熱貫流率、修正日射熱取得率、室内側表面熱流を計算
            h_cv_buf, h_rv_buf, u_dash_buf, eta_dash_buf, q_room_side_buf \
                = vws.get_vent_wall_performance_factor_by_simplified_calculation_no_03(parm=parms, h_out=h_out)

            # 配列に格納
            h_cv.append(h_cv_buf)
            h_rv.append(h_rv_buf)
            u_dash.append(u_dash_buf)
            eta_dash.append(eta_dash_buf)
            q_room_side.append(q_room_side_buf)

    # 計算結果をDataFrameに追加
    df['theta_sat'] = theta_sat
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['u_dash'] = u_dash
    df['eta_dash'] = eta_dash
    df['q_room_side'] = q_room_side

    return df


def get_wall_status_data_by_simplified_calculation_no_04() -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、簡易計算法案No.4（簡易計算法案No.3をさらに簡略化）による計算結果を保有するDataFrameを作成する

    :param: なし
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    df = get_parameter_data_frame()
    
    # 固定値の設定
    h_out = global_number.get_h_out()

    # 計算結果格納用配列を用意
    theta_sat = []          # 相当外気温度[℃]
    h_cv = []               # 通気層の対流熱伝達率[W/(m2・K)]
    h_rv = []               # 通気層の放射熱伝達率[W/(m2・K)]
    u_dash = []             # 修正熱貫流率[W/(m2・K)]
    eta_dash = []           # 修正日射熱取得率[-]
    q_room_side = []        # 室内表面熱流[W/m2]

    # エラーログ出力用の設定
    log = Log()
    saved_handler = np.seterrcall(log)

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
        for row in df.itertuples():
            print(row[0])
            # パラメータを設定
            parms = (vws.Parameters(theta_e=row.theta_e,
                                   theta_r=row.theta_r,
                                   J_surf=row.j_surf,
                                   a_surf=row.a_surf,
                                   C_1=row.C_1,
                                   C_2=row.C_2,
                                   l_h=row.l_h,
                                   l_w=row.l_w,
                                   l_d=row.l_d,
                                   angle=row.angle,
                                   v_a=row.v_a,
                                   l_s=row.l_s,
                                   emissivity_1=row.emissivity_1,
                                   emissivity_2=row.emissivity_2))

            # 相当外気温度を計算
            theta_sat.append(epf.get_theta_SAT(row.theta_e, row.a_surf, row.j_surf, h_out))

            # 対流熱伝達率、放射熱伝達率、修正熱貫流率、修正日射熱取得率、室内側表面熱流を計算
            h_cv_buf, h_rv_buf, u_dash_buf, eta_dash_buf, q_room_side_buf \
                = vws.get_vent_wall_performance_factor_by_simplified_calculation_no_04(parm=parms, h_out=h_out)

            # 配列に格納
            h_cv.append(h_cv_buf)
            h_rv.append(h_rv_buf)
            u_dash.append(u_dash_buf)
            eta_dash.append(eta_dash_buf)
            q_room_side.append(q_room_side_buf)

    # 計算結果をDataFrameに追加
    df['theta_sat'] = theta_sat
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['u_dash'] = u_dash
    df['eta_dash'] = eta_dash
    df['q_room_side'] = q_room_side

    return df


def dump_csv_all_case_result():
    # 総当たりのパラメータと計算結果を取得し、CSVに出力

    # 詳細計算
    print("Detailed Calculation")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation("detailed", "detailed"))
    df.to_csv("wall_status_data_frame_detailed.csv")

    # 放射熱伝達率の検証： 冬期条件の簡易計算
    print("Simplified Calculation: h_rv_winter")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_winter"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_winter.csv")

    # 放射熱伝達率の検証： 夏期条件の簡易計算
    print("Simplified Calculation: h_rv_summer")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_summer"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_summer.csv")

    # 放射熱伝達率の検証： 放射熱伝達率ゼロ
    print("Simplified Calculation: h_rv_zero")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_zero"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_zero.csv")

    # 放射熱伝達率の検証：　通年の簡易計算
    print("Simplified Calculation: h_rv_all_season")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_all_season"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_all_season.csv")

    # 対流熱伝達率の検証： 冬期条件の簡易計算
    print("Simplified Calculation: h_cv_winter")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="simplified_winter", calc_mode_h_rv="detailed"))
    df.to_csv("wall_status_data_frame_h_cv_simplified_winter.csv")

    # 対流熱伝達率の検証： 夏期条件の簡易計算
    print("Simplified Calculation: h_cv_summer")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="simplified_summer", calc_mode_h_rv="detailed"))
    df.to_csv("wall_status_data_frame_h_cv_simplified_summer.csv")

    # 対流熱伝達率の検証：　通年の簡易計算
    print("Simplified Calculation: h_cv_all_season")
    df = pd.DataFrame(get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="simplified_all_season", calc_mode_h_rv="detailed"))
    df.to_csv("wall_status_data_frame_h_cv_simplified_all_season.csv")

    # 簡易計算法案No.1（簡易版の行列式）による計算
    print("Simplified Calculation No.1")
    df = pd.DataFrame(get_wall_status_data_by_simplified_calculation_no_01())
    df.to_csv("wall_status_data_frame_simplified_calculation_no01.csv")

    # 簡易計算法案No.2（簡易式）による計算
    print("Simplified Calculation No.2")
    df = pd.DataFrame(get_wall_status_data_by_simplified_calculation_no_02())
    df.to_csv("wall_status_data_frame_simplified_calculation_no02.csv")

    # 簡易計算法案No.3（通気層を有する壁体の修正熱貫流率、修正日射熱取得率から室内表面熱流を求める）による計算
    print("Simplified Calculation No.3")
    df = pd.DataFrame(get_wall_status_data_by_simplified_calculation_no_03())
    df.to_csv("wall_status_data_frame_simplified_calculation_no03.csv")

    # 簡易計算法案No.4（簡易計算法案No.3をさらに簡略化）による計算
    print("Simplified Calculation No.4")
    df = pd.DataFrame(get_wall_status_data_by_simplified_calculation_no_04())
    df.to_csv("wall_status_data_frame_simplified_calculation_no04.csv")


if __name__ == '__main__':

    dump_csv_all_case_result()


# デバッグ用
# dump_csv_all_case_result()
# print(np.median([-20,40]))
