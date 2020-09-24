import math
import numpy as np
import heat_transfer_coefficient
from dataclasses import dataclass


@dataclass
class Parameters:

    # 外気温度, degree C
    theta_e: float

    # 室内温度,　degree C
    theta_r: float

    # 外気側表面に入射する日射量, W/m2
    J_surf: float

    # 外気側表面日射吸収率
    a_surf: float

    # 外気側部材の熱コンダクタンス,W/(m2・K)
    C_1: float

    # 室内側部材の熱コンダクタンス, W/(m2・K)
    C_2: float

    # 通気層の長さ, m
    l_h: float

    # 通気層の幅, m
    l_w: float

    # 通気層の厚さ, m
    l_d: float

    # 通気層の傾斜角, °
    angle: float

    # 通気層の平均風速, m/s
    # Note: 通気層の風速は計算により求める方法もあるが、ひとまず与条件とする
    v_a: float

    # 通気胴縁または垂木の間隔, m
    l_s: float

    # 通気層に面する面1の放射率, -
    emissivity_1: float

    # 通気層に面する面2の放射率, -
    emissivity_2: float


# 通気層内の各層の温度の計算
def calculate_layer_temperatures(parm):

    # 固定値を設定
    # Note: パラメータごとにこれらの値を計算する方法もあるが、ひとまず固定値とする
    h_out = 25.0    # 室外側総合熱伝達率, W/(m2・K)
    h_in = 9.0      # 室内側総合熱伝達率, W/(m2・K)
    c_a = 1006     # 空気の定圧比熱, J/(kg・K)
    # Note: 仕様書に記載間違いあり（誤：1.006　正：1006)
    rho_a = 1.2     # 空気の密度, kg/m3

    # 相当外気温度を計算
    theta_SAT = parm.theta_e + (parm.a_surf * parm.J_surf) / h_out

    # 行列の初期化
    matrix_coeff = np.zeros(shape=(5, 5))
    matrix_temp = np.zeros(shape=(5, 1))
    matrix_temp_prev = np.zeros(shape=(5, 1))
    matrix_const = np.zeros(shape=(5, 1))

    # 通気層内の各点の温度の初期値を設定
    matrix_temp_prev[0][0] = parm.theta_e
    matrix_temp_prev[1][0] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 3)
    matrix_temp_prev[2][0] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 2)
    matrix_temp_prev[3][0] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 1)
    matrix_temp_prev[4][0] = (matrix_temp_prev[1][0] + matrix_temp_prev[2][0]) / 2

    # 収束計算により変化しない要素を設定
    matrix_coeff[0][0] = h_out + parm.C_1
    matrix_coeff[0][1] = -parm.C_1
    matrix_coeff[1][0] = parm.C_1
    matrix_coeff[2][3] = parm.C_2
    matrix_coeff[3][2] = parm.C_2
    matrix_coeff[3][3] = -(h_in + parm.C_2)

    matrix_const[0][0] = h_out * theta_SAT
    # Note: 仕様書に記載間違いあり（誤：theta_SAT　正：h_out * theta_SAT)
    matrix_const[3][0] = -h_in * parm.theta_r
    matrix_const[4][0] = parm.theta_e   # Note: この処理は不要

    # 収束計算の開始
    # Note: 100回で収束しないときどうするか？
    for i in range(100):

        # 通気層内の表面温度を設定
        theta_1 = matrix_temp_prev[1][0]
        theta_2 = matrix_temp_prev[2][0]

        # 対流熱伝達率の計算
        h_cv = heat_transfer_coefficient.convective_heat_transfer_coefficient(parm.v_a, theta_1, theta_2, parm.angle, parm.l_h, parm.l_d)

        # 有効放射率の計算
        effective_emissivity = heat_transfer_coefficient.effective_emissivity_parallel(parm.emissivity_1, parm.emissivity_2)

        # 放射熱伝達率の計算
        h_rv = heat_transfer_coefficient.radiative_heat_transfer_coefficient(theta_1, theta_2, effective_emissivity)

        # 通気風量の計算
        v_vent = parm.v_a * parm.l_d * parm.l_w

        # 通気層の平均空気温度の計算用の値を設定
        beta = (2 * h_cv * parm.l_w) / (c_a * rho_a * v_vent)

        # 行列に値を設定
        matrix_coeff[1][1] = -(h_cv + h_rv + parm.C_1)
        matrix_coeff[1][2] = h_rv
        matrix_coeff[1][4] = h_cv
        matrix_coeff[2][1] = h_rv
        matrix_coeff[2][2] = -(h_cv + h_rv + parm.C_2)
        matrix_coeff[2][4] = h_cv
        matrix_coeff[4][1] = (beta * parm.l_h + math.exp(-beta * parm.l_h) - 1) / 2
        matrix_coeff[4][2] = (beta * parm.l_h + math.exp(-beta * parm.l_h) - 1) / 2
        matrix_coeff[4][4] = -beta * parm.l_h
        matrix_const[4][0] = (math.exp(-beta * parm.l_h) - 1) * parm.theta_e

        # matrix_coeffの逆行列を計算
        inv_matrix_coeff = np.linalg.inv(matrix_coeff)

        # matrix_temp (通気層内の各点の温度）の計算
        matrix_temp = np.matmul(inv_matrix_coeff, matrix_const)

        # matrix_temp のすべての要素が近似的に近い値か判定する
        if np.allclose(matrix_temp, matrix_temp_prev):
            print("Done!")
            print("θ_e:", parm.theta_e)
            print("θ_sat:", theta_SAT)
            print("θ_out,surf:", matrix_temp[0][0])
            print("θ_1,surf:", matrix_temp[1][0])
            print("θ_2,surf:", matrix_temp[2][0])
            print("θ_in,surf:", matrix_temp[3][0])
            print("θ_as:", matrix_temp[4][0])
            print("θ_r:", parm.theta_r)
            break
        else:
            # 収束しなかったときは各要素を平均値で置き換える
            for j in range(5):
                matrix_temp_prev[j][0] = (matrix_temp[j][0] + matrix_temp_prev[j][0]) / 2

    return matrix_temp, h_rv, h_cv


# 通気層を有する壁体の熱貫流率の計算（W/(m2・K)）
# Note: 日射熱取得率も同時に計算するように変更したため、関数名を変更する必要がある
def overall_heat_transfer_coefficient(parm):

    # 固定値を設定
    # Note:固定値はグローバル変数にするのがよいか
    h_out = 25.0    # 室外側総合熱伝達率, W/(m2・K)
    h_in = 9.0      # 室内側総合熱伝達率, W/(m2・K)
    r_s_e = 0.11    # 室外側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    r_s_r = 0.11    # 室内側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    # Note: 省エネ基準の規定による表面熱伝達抵抗は、屋根、外壁で値が異なる。どのように判定するかが課題。

    # 相当外気温度を計算
    theta_SAT = parm.theta_e + (parm.a_surf * parm.J_surf) / h_out

    # 各層の温度を計算
    matrix_temp, h_rv, h_cv = calculate_layer_temperatures(parm)

    # 省エネ基準でのU値を計算
    u_s = 1/(r_s_e + 1/parm.C_2 + r_s_r)

    # 温度、風速依存の熱伝達率を使用したU値に修正
    u_s_dash = 1/(1/u_s - r_s_e + 1/(h_rv + h_cv))

    # 通気層を有する壁体の熱貫流率(W/(m2・K))を計算
    u_e = u_s_dash * (matrix_temp[4][0] - parm.theta_r) / (theta_SAT - parm.theta_r)

    # 通気層を有する壁体の日射熱取得率(-)を計算
    eta_e = u_s_dash * parm.a_surf/h_out * (matrix_temp[4][0] - parm.theta_r) / (theta_SAT - parm.theta_r)

    return u_e, eta_e


# 通気層を有する壁体の日射熱取得率(-)の計算
# Note: この関数は不要
def solar_heat_gain_coefficient(parm):

    # 固定値を設定
    # Note:固定値はグローバル変数にするのがよいか
    h_in = 9.0  # 室内側総合熱伝達率, W/(m2・K)

    # 内外温度差を0とする（室内温度 = 外気温度とする）
    parm.theta_r = parm.theta_e

    # 各層の温度を計算
    matrix_temp, h_rv, h_cv = calculate_layer_temperatures(parm)

    # 通気層を有する壁体の日射熱取得率(-)を計算
    eta_e = h_in * (matrix_temp[3][0] - parm.theta_r) / parm.J_surf

    return eta_e


# デバッグ用
# parm_1 = Parameters(20, 25, 500, 0.9, 10, 0.5, 6.0, 0.45, 0.018, 90, 0.2, 0.45, 0.9, 0.9)
# u_e, eta_e = overall_heat_transfer_coefficient(parm_1)
# print("通気層を有する壁体の熱貫流率U_e:", u_e)
# print("通気層を有する壁体の日射熱取得率η_e):", eta_e)
