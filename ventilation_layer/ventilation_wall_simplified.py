import math
import numpy as np
from dataclasses import dataclass


import ventilation_layer.heat_transfer_coefficient as htc
import ventilation_layer.ventilation_wall as vw
import ventilation_layer.envelope_performance_factors as epf
from ventilation_layer.global_number import get_c_air, get_rho_air, get_h_out
import ventilation_layer.parameters as pm
from ventilation_layer.heat_transfer_coefficient import HCVMode, HCVModeDetail, HCVModeSimple, HRVMode, HRVModeDetail, HRVModeSimple


def get_datas(pmo: pm.ParametersOptions, hcv_mode: HCVMode, hrv_mode: HRVMode):

    # 計算結果格納用配列を用意
    h_cv = []               # 通気層の対流熱伝達率[W/(m2・K)]
    h_rv = []               # 通気層の放射熱伝達率[W/(m2・K)]
    u_dash = []             # 修正熱貫流率[W/(m2・K)]
    eta_dash = []           # 修正日射熱取得率[-]
    q_room_side = []        # 室内表面熱流[W/m2]

    for parms in pmo.get_parameters_list():

        # 対流熱伝達率、放射熱伝達率、修正熱貫流率、修正日射熱取得率、室内側表面熱流を計算
        h_cv_buf, h_rv_buf, u_dash_buf, eta_dash_buf, q_room_side_buf = calc(
            hcv_mode=hcv_mode,
            hrv_mode=hrv_mode,
            theta_e=parms.theta_e,
            theta_r=parms.theta_r,
            j_surf=parms.j_surf,
            a_surf=parms.a_surf,
            C_1=parms.C_1,
            C_2=parms.C_2,
            l_h=parms.l_h,
            l_w=parms.l_w,
            l_d=parms.l_d,
            v_a=parms.v_a,
            emissivity_1=parms.emissivity_1,
            emissivity_2=parms.emissivity_2
        )

        # 配列に格納
        h_cv.append(h_cv_buf)
        h_rv.append(h_rv_buf)
        u_dash.append(u_dash_buf)
        eta_dash.append(eta_dash_buf)
        q_room_side.append(q_room_side_buf)

    df = pmo.get_df()

    # 計算結果をDataFrameに追加
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['u_dash'] = u_dash
    df['eta_dash'] = eta_dash
    df['q_room_side'] = q_room_side

    return df


def calc(
        hcv_mode: HCVMode,
        hrv_mode: HRVMode,
        theta_e: float,
        theta_r: float,
        j_surf: float,
        a_surf: float,
        C_1: float,
        C_2: float,
        l_h: float,
        l_w: float,
        l_d: float,
        v_a: float,
        emissivity_1: float,
        emissivity_2: float
    ):
    """修正U値と修正η値を求める。(対流・放射熱伝達率、室内側熱流も計算する。)
    """

    # 有効放射率の計算
    e_eff = htc.get_e(emissivity_1, emissivity_2)

    # 対流熱伝達率、放射熱伝達率の計算
    h_cv = htc.get_h_cv(hcv_mode=hcv_mode, v_a=v_a)
    h_rv = htc.get_h_rv(hrv_mode=hrv_mode, eps_eff=e_eff)

    u_dash, eta_dash = get_u_eta(
        hcv_mode=hcv_mode,
        hrv_mode=hrv_mode,
        theta_e=theta_e,
        a_surf=a_surf,
        C_1=C_1,
        C_2=C_2,
        l_h=l_h,
        l_w=l_w,
        l_d=l_d,
        v_a=v_a,
        emissivity_1=emissivity_1,
        emissivity_2=emissivity_2
    )

    # 室内表面熱流を計算
    q_room_side = u_dash * (theta_e - theta_r) + eta_dash * j_surf

    return h_cv, h_rv, u_dash, eta_dash, q_room_side


def get_u_eta(
        hcv_mode: HCVMode,
        hrv_mode: HRVMode,
        theta_e: float,
        a_surf: float,
        C_1: float,
        C_2: float,
        l_h: float,
        l_w: float,
        l_d: float,
        v_a: float,
        emissivity_1: float,
        emissivity_2: float
    ) -> tuple[float, float]:
    """修正U値と修正η値を求める。

    """

    # 室外側熱伝達率, W/(m2 K)
    h_out = get_h_out()

    # 有効放射率の計算
    e_eff = htc.get_e(emissivity_1, emissivity_2)

    # 対流熱伝達率、放射熱伝達率の計算
    h_cv = htc.get_h_cv(hcv_mode=hcv_mode, v_a=v_a)
    h_rv = htc.get_h_rv(hrv_mode=hrv_mode, eps_eff=e_eff)

    # 熱伝達率の計算
    h_v = 2.0 * h_rv + h_cv

    # 通気風量の計算
    v_vent = v_a * l_d * l_w

    # 通気層の熱抵抗の値を設定
    if v_a > 0.0:
        beta = (2 * h_cv * l_w) / (get_c_air() * get_rho_air(theta_e) * v_vent)
        epc_s = 1.0 / l_h * 1.0 / beta * (math.exp(-beta * l_h) - 1.0)
        epc_s_dash = - ((2.0 * h_cv) * epc_s) / (1.0 + epc_s)
        r_e = 1.0 / epc_s_dash + h_rv / (h_v * h_cv)
        h_e = 1 / r_e
    else:
        h_e = 0.0
    
    # 熱抵抗を設定
    u_o_s = 1.0 / epf.get_r_o(C_1)
    u_i_s = 1.0 / epf.get_r_i(C_2)

    # R_sat
    r_sat = 1 / u_o_s + 1 / h_v

    p_1 = 1.0 / (1 / r_sat + h_e) + 1.0 / h_v

    p_2 = 1.0 / (r_sat * h_e + 1)
    
    # 修正U値を計算
    u_dash = 1.0 / (p_1 + 1.0 / u_i_s)

    # 修正η値を計算
    eta_dash = p_2 * u_dash * a_surf / h_out

    return u_dash, eta_dash
