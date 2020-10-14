import math
from scipy import optimize
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


@dataclass
class WallStatusValues:

    # 通気層内の各点の温度, degree C
    matrix_temp: np.zeros(shape=(5, 1))

    # 対流熱伝達率, W/(m2・K)
    h_cv: float

    # 放射熱伝達率, W/(m2・K)
    h_rv: float


# 熱収支式を解く関数
def get_heat_balance(matrix_temp: np.zeros(shape=(5, 1)), parm: Parameters) -> np.zeros(shape=(5, 1)):

    # 固定値を設定
    # Note: パラメータごとにこれらの値を計算する方法もあるが、ひとまず固定値とする
    h_out = 25.0  # 室外側総合熱伝達率, W/(m2・K)
    h_in = 9.0  # 室内側総合熱伝達率, W/(m2・K)
    c_a = 1006  # 空気の定圧比熱, J/(kg・K)
    rho_a = 1.2  # 空気の密度, kg/m3

    # 相当外気温度を計算
    theta_SAT = parm.theta_e + (parm.a_surf * parm.J_surf) / h_out

    # 行列の初期化
    matrix_coeff = np.zeros(shape=(5, 5))
    matrix_const = np.zeros(shape=(5, 1))

    # 通気層内の表面温度を設定
    theta_1 = matrix_temp[1][0]
    theta_2 = matrix_temp[2][0]

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
    matrix_coeff[0][0] = h_out + parm.C_1
    matrix_coeff[0][1] = -parm.C_1
    matrix_coeff[1][0] = parm.C_1
    matrix_coeff[1][1] = -(h_cv + h_rv + parm.C_1)
    matrix_coeff[1][2] = h_rv
    matrix_coeff[1][4] = h_cv
    matrix_coeff[2][1] = h_rv
    matrix_coeff[2][2] = -(h_cv + h_rv + parm.C_2)
    matrix_coeff[2][3] = parm.C_2
    matrix_coeff[2][4] = h_cv
    matrix_coeff[3][2] = parm.C_2
    matrix_coeff[3][3] = -(h_in + parm.C_2)
    matrix_coeff[4][1] = (beta * parm.l_h + math.exp(-beta * parm.l_h) - 1) / 2
    matrix_coeff[4][2] = (beta * parm.l_h + math.exp(-beta * parm.l_h) - 1) / 2
    matrix_coeff[4][4] = -beta * parm.l_h
    matrix_const[0][0] = h_out * theta_SAT
    matrix_const[3][0] = -h_in * parm.theta_r
    matrix_const[4][0] = (math.exp(-beta * parm.l_h) - 1) * parm.theta_e

    # 熱収支を計算
    q_balance = np.matmul(matrix_coeff, matrix_temp) - matrix_const

    return q_balance


# 通気層の状態値を取得する
def get_wall_status_values(parm: Parameters) -> WallStatusValues:

    # 通気層内の各点の温度の初期値を設定
    matrix_temp = np.zeros(shape=(5, 1))
    matrix_temp[0][0] = parm.theta_e
    matrix_temp[1][0] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 3)
    matrix_temp[2][0] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 2)
    matrix_temp[3][0] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 1)
    matrix_temp[4][0] = (matrix_temp[1][0] + matrix_temp[2][0]) / 2

    # 通気層内の各点の熱収支式が成り立つときの各点の温度を取得
    answer_T = optimize.root(fun=get_heat_balance, x0=matrix_temp, args=parm, method='broyden1')
    matrix_temp_fixed = answer_T.x

    # 対流熱伝達率の計算
    h_cv = heat_transfer_coefficient.convective_heat_transfer_coefficient(parm.v_a, matrix_temp_fixed[1][0], matrix_temp_fixed[2][0], parm.angle,
                                                                          parm.l_h, parm.l_d)

    # 有効放射率の計算
    effective_emissivity = heat_transfer_coefficient.effective_emissivity_parallel(parm.emissivity_1, parm.emissivity_2)

    # 放射熱伝達率の計算
    h_rv = heat_transfer_coefficient.radiative_heat_transfer_coefficient(matrix_temp_fixed[1][0], matrix_temp_fixed[2][0], effective_emissivity)

    return WallStatusValues(matrix_temp = matrix_temp_fixed, h_cv = h_cv, h_rv = h_rv)


# デバッグ用
# parm_1: Parameters = Parameters(20, 25, 500, 0.9, 10, 0.5, 6.0, 0.45, 0.018, 90, 0.2, 0.45, 0.9, 0.9)
