import pandas as pd
import numpy as np

import ventilation_layer.global_number as global_number
import ventilation_layer.ventilation_wall as vw
import ventilation_layer.envelope_performance_factors as epf
import ventilation_layer.heat_transfer_coefficient as htc
import ventilation_layer.parameters as pm
from ventilation_layer.heat_transfer_coefficient import HCVMode, HCVModeDetail, HCVModeSimple, HRVMode, HRVModeDetail, HRVModeSimple


def get_wall_status_data_by_detailed_calculation(pms: pm.ParametersOptions, hcv_mode: HCVMode, hrv_mode: HRVMode) -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、各ケースの計算結果を保有するDataFrameを作成する

    :param calc_mode_h_cv: 対流熱伝達率の計算モード
    :param calc_mode_h_rv: 放射熱伝達率の計算モード
    :return: DataFrame
    """

    df = pms.get_df()

    # 固定値の設定
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
            
        # 通気層の状態値を取得
        results = [
            vw.get_wall_status_values(
                index=i,
                theta_e=pms.theta_e,
                theta_r=pms.theta_r,
                j_surf=pms.j_surf,
                a_surf=pms.a_surf,
                c_1=pms.C_1,
                c_2=pms.C_2,
                l_h=pms.l_h,
                l_w=pms.l_w,
                l_d=pms.l_d,
                angle=pms.angle,
                v_a=pms.v_a,
                eps_1=pms.emissivity_1,
                eps_2=pms.emissivity_2,
                hcv_mode=hcv_mode,
                hrv_mode=hrv_mode,
                h_out=h_out,
                h_in=h_in
            ) for i, pms in enumerate(pms.get_parameters_list())
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
        h_cv = np.vectorize(htc.get_h_cv)(hcv_mode=hcv_mode, v_a=v_a, theta_1=theta_1, theta_2=theta_2, angle=angle, l_h=l_h, l_d=l_d)

        # the radiative heat transfer coefficient, W/m2K
        h_rv = np.vectorize(htc.get_h_rv)(hrv_mode=hrv_mode, eps_eff=eps_eff, theta_1=theta_1, theta_2=theta_2)

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

        df["delta_theta_as_e-theta_r"] = df["theta_as_e"] - df["theta_r"]
        df["delta_theta_sat-theta_r"] = df["theta_sat"] - df["theta_r"]
    

    return df
