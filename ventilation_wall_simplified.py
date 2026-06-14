import math
import numpy as np
from dataclasses import dataclass


import ventilation_layer.heat_transfer_coefficient as htc
import ventilation_layer.ventilation_wall as vw
import envelope_performance_factors as epf
from ventilation_layer.global_number import get_c_air, get_rho_air, get_h_out
import ventilation_layer.parameters as pm


def calc(parm: pm.Parameters):
    """
    簡易計算法案No.3：通気層を有する壁体の修正熱貫流率、修正日射熱取得率、室内表面熱流を求める関数

    :param parm:    計算条件パラメータ群
    :param h_out:   室外側総合熱伝達率[W/(m2・K)]
    :return:        修正熱貫流率[W/(m2・K)], 修正日射熱取得率[-], 室内表面熱流[W/m2]
    """

    # 室外側熱伝達率, W/(m2 K)
    h_out = get_h_out()

    # 有効放射率の計算
    effective_emissivity = htc.get_e(parm.emissivity_1, parm.emissivity_2)

    # 対流熱伝達率、放射熱伝達率の計算
    if parm.theta_r == 20.0:
        h_cv = htc.get_h_cv(calc_mode="simplified_winter", v_a=parm.v_a)
        h_rv = htc.get_h_rv(calc_mode="simplified_winter", eps_eff=effective_emissivity)
    else:
        h_cv = htc.get_h_cv(calc_mode="simplified_winter", v_a=parm.v_a)
        h_rv = htc.get_h_rv(calc_mode="simplified_winter", eps_eff=effective_emissivity)

    # 熱伝達率の計算
    h_v = 2.0 * h_rv + h_cv

    # 通気風量の計算
    v_vent = parm.v_a * parm.l_d * parm.l_w

    # 通気層の熱抵抗の値を設定
    if parm.v_a > 0.0:
        beta = (2 * h_cv * parm.l_w) / (get_c_air() * get_rho_air(parm.theta_e) * v_vent)
        epc_s = 1.0 / parm.l_h * 1.0 / beta * (math.exp(-beta * parm.l_h) - 1.0)
        epc_s_dash = - ((2.0 * h_cv) * epc_s) / (1.0 + epc_s)
        r_e = 1.0 / epc_s_dash + h_rv / (h_v * h_cv)
        h_e = 1 / r_e
    else:
        h_e = 0.0
    
    # 熱抵抗を設定
    u_o_s = 1.0 / epf.get_r_o(parm.C_1)
    u_i_s = 1.0 / epf.get_r_i(parm.C_2)

    # R_sat
    r_sat = 1 / u_o_s + 1 / h_v

    p_1 = 1.0 / (1 / r_sat + h_e) + 1.0 / h_v

    p_2 = 1.0 / (r_sat * h_e + 1)
    
    # 修正U値を計算
    u_dash = 1.0 / (p_1 + 1.0 / u_i_s)

    # 修正η値を計算
    eta_dash = p_2 * u_dash * parm.a_surf / h_out

    # 室内表面熱流を計算
    q_room_side = u_dash * (parm.theta_e - parm.theta_r) + eta_dash * parm.j_surf

    return h_cv, h_rv, u_dash, eta_dash, q_room_side

