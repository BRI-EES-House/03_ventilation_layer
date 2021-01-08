import math
import numpy as np
import heat_transfer_coefficient as htc
import ventilation_wall as vw
import envelope_performance_factors as epf
from global_number import get_c_air, get_rho_air


def get_vent_wall_temperature_by_simplified_matrix(parm: vw.Parameters, h_out: float, h_in: float) -> np.zeros(3):
    """
    簡易版の行列式により各部位の温度を求める関数

    :param parm:    計算条件パラメータ群
    :param h_out:   室外側総合熱伝達率[W/(m2・K)]
    :param h_in:    室内側総合熱伝達率[W/(m2・K)]
    :return:        各部位の温度(3×1の行列)[degC], 対流熱伝達率[W/(m2・K)], 放射熱伝達率[W/(m2・K)], 室内側から通気層表面までの熱抵抗[(m2・K)/W]
    """

    # 相当外気温度を計算
    theta_SAT = epf.get_theta_SAT(theta_e=parm.theta_e, a_surf=parm.a_surf, j_surf=parm.J_surf, h_out=h_out)

    # 行列の初期化
    matrix_coeff = np.zeros(shape=(3, 3))
    matrix_const = np.zeros(3)

    # 対流熱伝達率の計算
    h_cv = htc.convective_heat_transfer_coefficient_simplified_all_season(v_a=parm.v_a)

    # 有効放射率の計算
    effective_emissivity = htc.effective_emissivity_parallel(parm.emissivity_1, parm.emissivity_2)

    # 放射熱伝達率の計算
    h_rv = htc.radiative_heat_transfer_coefficient_simplified_all_season(effective_emissivity=effective_emissivity)

    # 通気風量の計算
    v_vent = parm.v_a * parm.l_d * parm.l_w

    # 通気層の平均空気温度の計算用の値を設定
    beta = 0.0
    epc_s = 0.0
    if parm.v_a > 0.0:
        beta = (2 * h_cv * parm.l_w) / (get_c_air(parm.theta_e) * get_rho_air(parm.theta_e) * v_vent)
        epc_s = 1.0 / parm.l_h * 1.0 / beta * (math.exp(-beta * parm.l_h) - 1)

    # 熱抵抗を設定
    R_o = epf.get_r_o(parm.C_1)
    R_i = epf.get_r_i(parm.C_2)

    # 行列に値を設定
    matrix_coeff[0][0] = 1.0/R_o + h_cv + h_rv
    matrix_coeff[0][1] = -h_cv
    matrix_coeff[0][2] = -h_rv
    matrix_coeff[1][0] = (1.0 + epc_s)/2.0
    matrix_coeff[1][1] = -1.0
    matrix_coeff[1][2] = (1.0 + epc_s)/2.0
    matrix_coeff[2][0] = -h_rv
    matrix_coeff[2][1] = -h_cv
    matrix_coeff[2][2] = 1.0/R_i + h_cv + h_rv

    matrix_const[0] = (1.0/R_o) * theta_SAT
    matrix_const[1] = epc_s * parm.theta_e
    matrix_const[2] = (1.0/R_i) * parm.theta_r

    # 逆行列を計算
    matrix_coeff_inv = np.linalg.inv(matrix_coeff)

    # 各部位の温度を計算
    matrix_temp = np.matmul(matrix_coeff_inv, matrix_const)

    return matrix_temp, h_cv, h_rv, R_i


def get_vent_wall_temperature_by_simplified_equation(parm: vw.Parameters, h_out: float):
    """
    簡易式により通気層の平均温度を求める関数

    :param parm:    計算条件パラメータ群
    :param h_out:   室外側総合熱伝達率[W/(m2・K)]
    :return:        通気層の平均温度[degC], 室外側から通気層までの熱貫流率[W/(m2・K)], 室内側から通気層までの熱貫流率[W/(m2・K)]
    """
    # 相当外気温度を計算
    theta_sat = epf.get_theta_SAT(theta_e=parm.theta_e, a_surf=parm.a_surf, j_surf=parm.J_surf, h_out=h_out)

    # 対流熱伝達率の計算
    h_cv = htc.convective_heat_transfer_coefficient_simplified_all_season(v_a=parm.v_a)

    # 有効放射率の計算
    effective_emissivity = htc.effective_emissivity_parallel(parm.emissivity_1, parm.emissivity_2)

    # 放射熱伝達率の計算
    h_rv = htc.radiative_heat_transfer_coefficient_simplified_all_season(effective_emissivity=effective_emissivity)

    # 室外側から通気層までの熱貫流率、室内側から通気層までの熱貫流率を計算
    u_o = epf.get_u_o(parm.C_1, h_cv, h_rv)
    u_i = epf.get_u_i(parm.C_2, h_cv, h_rv)

    # theta_weを計算
    theta_we = (u_o * theta_sat + u_i * parm.theta_r) / (u_o + u_i)

    # 通気風量の計算
    v_vent = parm.v_a * parm.l_d * parm.l_w

    # 通気層の平均空気温度の計算用の値を設定
    if parm.v_a > 0.0:
        w_h = (u_o + u_i) / (get_c_air(parm.theta_e) * get_rho_air(parm.theta_e) * v_vent)
        epc = 1.0 - math.exp(-w_h * parm.l_h)
        x = 1.0 - epc / (w_h * parm.l_h)
    else:
        w_h = 0.0
        epc = 0.0
        x = 1.0

    theta_as_ave = (1.0 - x) * parm.theta_e + x * theta_we

    return theta_as_ave, u_o, u_i


# デバッグ用
# parm_1: vw.Parameters = vw.Parameters(10, 20, 500, 1.0, 50.25, 2.55, 3.0, 0.05, 0.05, 45.0, 0.5, 0.45, 0.9, 0.9)
# temps = get_vent_wall_temperature(parm_1, h_out=25.0, h_in=9.0)
# print(temps)