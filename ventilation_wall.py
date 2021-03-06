import math
from scipy import optimize
import numpy as np
import heat_transfer_coefficient
from dataclasses import dataclass
from global_number import get_c_air, get_rho_air


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

    # 通気層の傾斜角, degree
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

    # 各層の熱収支
    matrix_heat_balance: np.zeros(shape=(5, 1))

    # 対流熱伝達率, W/(m2・K)
    h_cv: float

    # 放射熱伝達率, W/(m2・K)
    h_rv: float

    # 最適化が正常に終了したかどうか
    is_optimize_succeed: bool

    # 最適化の終了ステータス
    optimize_status: int

    # 最適化の終了メッセージ
    optimize_message: str


def get_heat_balance(matrix_temp: np.zeros(5), parm: Parameters, calc_mode_h_cv: str, calc_mode_h_rv: str,
                     h_out: float, h_in: float) -> np.zeros(5):
    """
    熱収支式を解く関数

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param parm:        計算条件パラメータ群
    :param calc_mode_h_cv:   対流熱伝達率の計算モード
    :param calc_mode_h_rv:   放射熱伝達率の計算モード
    :param h_out:       室外側総合熱伝達率, W/(m2・K)
    :param h_in:        室内側総合熱伝達率, W/(m2・K)
    :return: 　         各層の熱収支, W/m2
    """

    # 相当外気温度を計算
    theta_SAT = parm.theta_e + (parm.a_surf * parm.J_surf) / h_out

    # 行列の初期化
    matrix_coeff = np.zeros(shape=(5, 5))
    matrix_const = np.zeros(5)

    # 通気層内の表面温度を設定
    theta_1 = matrix_temp[1]
    theta_2 = matrix_temp[2]

    # 対流熱伝達率の計算
    h_cv = heat_transfer_coefficient.get_convective_heat_transfer_coefficient(calc_mode_h_cv, parm.v_a, theta_1, theta_2, parm.angle, parm.l_h, parm.l_d)

    # 有効放射率の計算
    effective_emissivity = heat_transfer_coefficient.effective_emissivity_parallel(parm.emissivity_1, parm.emissivity_2)

    # 放射熱伝達率の計算
    h_rv = heat_transfer_coefficient.get_radiative_heat_transfer_coefficient(calc_mode_h_rv, theta_1, theta_2, effective_emissivity)

    # 通気風量の計算
    v_vent = parm.v_a * parm.l_d * parm.l_w

    # 通気層の平均空気温度の計算用の値を設定
    beta = 0.0
    if parm.v_a > 0.0:
        beta = (2 * h_cv * parm.l_w) / (get_c_air(matrix_temp[4]) * get_rho_air(matrix_temp[4]) * v_vent)

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
    matrix_coeff[4][4] = -1.0
    matrix_const[0] = h_out * theta_SAT
    matrix_const[3] = -h_in * parm.theta_r

    if parm.v_a > 0.0:
        matrix_coeff[4][1] = (1.0 + 1.0 / parm.l_h * 1.0 / beta * (math.exp(-beta * parm.l_h) - 1)) / 2
        matrix_coeff[4][2] = (1.0 + 1.0 / parm.l_h * 1.0 / beta * (math.exp(-beta * parm.l_h) - 1)) / 2
        matrix_const[4] = 1.0 / parm.l_h * 1.0 / beta * (math.exp(-beta * parm.l_h) - 1) * parm.theta_e
    else:
        matrix_coeff[4][1] = 0.5
        matrix_coeff[4][2] = 0.5
        matrix_const[4] = 0.0

    # 熱収支を計算
    q_balance = np.matmul(matrix_coeff, matrix_temp) - matrix_const

    return q_balance


def get_wall_status_values(parm: Parameters, calc_mode_h_cv: str, calc_mode_h_rv: str,
                           h_out: float, h_in: float) -> WallStatusValues:
    """
    通気層の状態値を取得する

    :param parm: 計算条件パラメータ群
    :param calc_mode_h_cv:   対流熱伝達率の計算モード
    :param calc_mode_h_rv:   放射熱伝達率の計算モード
    :param h_out: 室外側総合熱伝達率, W/(m2・K)
    :param h_in:  室内側総合熱伝達率, W/(m2・K)
    :return: 通気層の状態値（通気層の各層の温度、各層の熱収支、対流熱伝達率、放射熱伝達率、最適化の終了ステータス、終了メッセージ）
    """

    # 通気層内の各点の温度の初期値を設定
    matrix_temp = np.zeros(5)
    matrix_temp[0] = parm.theta_e
    matrix_temp[1] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 3)
    matrix_temp[2] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 2)
    matrix_temp[3] = parm.theta_e + (parm.theta_r - parm.theta_e) / (4 * 1)
    matrix_temp[4] = (matrix_temp[1] + matrix_temp[2]) / 2

    # 通気層内の各層の熱収支式の最適解を収束計算で求める
    optimize_result = optimize.root(fun=get_heat_balance, x0=matrix_temp, args=(parm, calc_mode_h_cv, calc_mode_h_rv, h_out, h_in), method='lm')

    # 収束した場合は各層の状態値を設定、収束しなかった場合はすべて無効（Nan）とする
    if optimize_result.success:

        # 熱収支式が成り立つときの各層の温度を取得
        matrix_temp_fixed = optimize_result.x

        # 各層の熱収支を計算
        heat_balance = get_heat_balance(matrix_temp_fixed, parm, calc_mode_h_cv, calc_mode_h_rv, h_out, h_in)

        # 対流熱伝達率の計算
        h_cv = heat_transfer_coefficient.get_convective_heat_transfer_coefficient(calc_mode=calc_mode_h_cv, v_a=parm.v_a, theta_1=matrix_temp_fixed[1], theta_2=matrix_temp_fixed[2], angle=parm.angle, l_h=parm.l_h, l_d=parm.l_d)

        # 有効放射率の計算
        effective_emissivity = heat_transfer_coefficient.effective_emissivity_parallel(emissivity_1=parm.emissivity_1, emissivity_2=parm.emissivity_2)

        # 放射熱伝達率の計算
        h_rv = heat_transfer_coefficient.get_radiative_heat_transfer_coefficient(calc_mode=calc_mode_h_rv, theta_1=matrix_temp_fixed[1], theta_2=matrix_temp_fixed[2], effective_emissivity=effective_emissivity)

    else:
        matrix_temp_fixed = np.full(5, np.nan)
        heat_balance = np.full(5, np.nan)
        h_cv = np.nan
        h_rv = np.nan

    return WallStatusValues(matrix_temp=matrix_temp_fixed, matrix_heat_balance=heat_balance, h_cv=h_cv, h_rv=h_rv,
                            is_optimize_succeed=optimize_result.success, optimize_status=optimize_result.status,
                            optimize_message=optimize_result.message
                            )


def get_heat_flow_0(matrix_temp: np.ndarray, param: Parameters, h_out: float) -> float:
    """
    各部温度から屋外側表面熱流を計算する

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param param:       計算条件パラメータ群
    :param h_out:       室外側総合熱伝達率, W/(m2・K)
    :return:            屋外側表面熱流, W/m2
    """

    # 相当外気温度を計算
    theta_sat = param.theta_e + (param.a_surf * param.J_surf) / h_out

    return h_out * (theta_sat - matrix_temp[0])


def get_heat_flow_1(matrix_temp: np.ndarray, param: Parameters) -> float:
    """
    各部温度から外装材伝導熱量を計算する

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param param:       計算条件パラメータ群
    :return:            外装材伝導熱量, W/m2
    """

    return param.C_1 * (matrix_temp[0] - matrix_temp[1])


def get_heat_flow_exhaust(matrix_temp: np.ndarray, param: Parameters, theta_as_in: float, h_cv: float) -> float:
    """
    通気層からの排気熱量の計算

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param param:       計算条件パラメータ群
    :param theta_as_in: 通気層への流入温度=外気温度, degC
    :param h_cv:        通気層の対流熱伝達率, W/m2K
    :return:            通気層の排気熱量, W/m2
    """

    if param.v_a > 0.0:

        # 通気風量の計算
        v_vent = param.v_a * param.l_d * param.l_w

        ec = math.exp(- 2.0 * h_cv * param.l_w * param.l_h / (get_c_air(matrix_temp[4]) * get_rho_air(matrix_temp[4]) * v_vent))

        # 出口温度の計算
        theta_out = (1.0 - ec) * (matrix_temp[1] + matrix_temp[2]) / 2.0 + ec * theta_as_in

        # 通気層の排気熱量
        return get_c_air(matrix_temp[4]) * get_rho_air(matrix_temp[4]) * v_vent * (theta_out - theta_as_in) / (param.l_w * param.l_h)

    else:
        return 0.0


def get_heat_flow_convect_vent_layer(matrix_temp: np.ndarray, param: Parameters, h_cv: float) -> float:
    """
    通気層内表面から通気層空気への対流熱量

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param param:       計算条件パラメータ群
    :param h_cv:        通気層対流熱伝達率, W/m2K
    :return:            通気層内表面から通気層空気への対流熱量、W/m2
    """

    return 2.0 * h_cv * ((matrix_temp[1] + matrix_temp[2]) / 2.0 - matrix_temp[4])


def get_heat_flow_2(matrix_temp: np.ndarray, h_cv: float, h_rv: float) -> tuple:
    """
    通気層熱伝達量の計算

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param h_cv:        通気層対流熱伝達率, W/m2K
    :param h_rv:        通気層放射熱伝達率, W/m2K
    :return:            通気層熱伝達量, W/m2
    """

    # 熱伝達量
    return (h_cv * (matrix_temp[1] - matrix_temp[2]), h_rv * (matrix_temp[1] - matrix_temp[2]))


def get_heat_flow_3(matrix_temp: np.ndarray, param: Parameters) -> float:
    """
    各部温度から断熱材+内装材伝導熱量を計算する

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param param:       計算条件パラメータ群
    :return:            断熱材+内装材伝導熱量, W/m2
    """

    return param.C_2 * (matrix_temp[2] - matrix_temp[3])


def get_heat_flow_4(matrix_temp: np.ndarray, param: Parameters, h_in: float) -> float:
    """
    各部温度から室内表面熱流を計算する

    :param matrix_temp: 各部温度計算結果 (5,1), degC
    :param param:       計算条件パラメータ群
    :param h_in:        室内側総合熱伝達率, W/(m2・K)
    :return:            断熱材+内装材伝導熱量, W/m2
    """

    return h_in * (matrix_temp[3] - param.theta_r)


# デバッグ用
# parm_1: Parameters = Parameters(-20, 20, 500, 1.0, 50.25, 2.55, 3.0, 0.05, 0.05, 45.0, 0.5, 0.45, 0.9, 0.9)
# status = get_wall_status_values(parm_1, h_out=25.0, h_in=9.0)
# print(status.is_optimize_succeed)
# print(status.optimize_status)
# print(status.optimize_message)
# print(status.matrix_heat_balance)
