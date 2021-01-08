import numpy as np
import ventilation_wall
import global_number


def overall_heat_transfer_coefficient(theta_e: float, theta_r: float, a_surf: float, j_surf: float, C_2: float,
                                      angle: float, theta_as_e: float, h_cv: float, h_rv: float) -> float:
    """
    通気層を考慮した壁体の相当熱貫流率[W/(m2・K)]を計算する

    :param theta_e: 外気温度[degC]
    :param theta_r: 室内温度[degC]
    :param a_surf:  外気側表面日射吸収率[-]
    :param j_surf:  外気側表面に入射する日射量[W/m2]
    :param C_2:     室内側部材の熱コンダクタンス[W/(m2・K)]
    :param angle:   通気層の傾斜角[°]
    :param theta_as_e: 通気層の等価温度[degC]
    :param h_cv:    通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:    通気層の放射熱伝達率[W/(m2・K)]
    :return:        通気層を考慮した壁体の相当熱貫流率[W/(m2・K)]
    """

    # 固定値を設定
    h_out = global_number.get_h_out()   # 室外側総合熱伝達率, W/(m2・K)

    # 相当外気温度を計算
    theta_SAT = get_theta_SAT(theta_e, a_surf, j_surf, h_out)

    # 温度、風速依存の熱伝達率を使用したU値に修正
    u_s_dash = get_u_s_dash(angle, C_2, h_cv, h_rv)

    # 通気層を有する壁体の相当熱貫流率を求めるための補正係数を計算
    k_e = get_k_e(theta_as_e, theta_r, theta_SAT)

    # 通気層を有する壁体の熱貫流率(W/(m2・K))を計算（k_eがnanのときはnanとする）
    if np.isnan(k_e):
        u_e = np.nan
    else:
        u_e = u_s_dash * k_e

    return u_e


def get_k_e(theta_as_e: float, theta_r: float, theta_SAT: float) -> float:
    """
    通気層を有する壁体の相当熱貫流率を求めるための補正係数（k_e）[-]を計算する

    :param theta_as_e:    通気層の等価温度[degC]
    :param theta_r:       室内温度[degC]
    :param theta_SAT:     相当外気温度[degC]
    :return:              通気層を有する壁体の相当熱貫流率を求めるための補正係数[-]
    """

    if abs(theta_SAT - theta_r) < 0.001:
        # 相当外気温度と室温が同じ値の場合、貫流熱がゼロになるので、Nanとする
        k_e = np.nan
    else:
        k_e = (theta_as_e - theta_r) / (theta_SAT - theta_r)

    return k_e


# 通気層を有する壁体の日射熱取得率(-)の計算
# TODO: parmを個別の変数に書き下す、固定値はglobal_numberから取得する、通気層の角度（angle）によって表面熱伝達抵抗を変更する
def solar_heat_gain_coefficient(parm: ventilation_wall.Parameters, theta_as_ave: float, h_cv: float, h_rv: float) -> float:
    """

    :param parm:
    :param theta_as_ave:
    :param h_cv:
    :param h_rv:
    :return:
    """
    # 固定値を設定
    # Note:固定値はグローバル変数にするのがよいか
    h_out = 25.0    # 室外側総合熱伝達率, W/(m2・K)
    r_s_e = 0.11  # 室外側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    r_s_r = 0.11  # 室内側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    # Note: 省エネ基準の規定による表面熱伝達抵抗は、屋根、外壁で値が異なる。どのように判定するかが課題。

    # 相当外気温度を計算
    theta_SAT = get_theta_SAT(parm.theta_e, parm.a_surf, parm.J_surf, h_out)

    # 温度、風速依存の熱伝達率を使用したU値を計算
    u_s_dash = get_u_s_dash(r_s_e, r_s_r, parm.C_2, h_cv, h_rv)

    # 通気層を有する壁体の日射熱取得率(-)を計算
    eta_e = u_s_dash * parm.a_surf/h_out * (theta_as_ave - parm.theta_r) / (theta_SAT - parm.theta_r)

    return eta_e


def get_theta_SAT(theta_e: float, a_surf: float, j_surf: float, h_out: float) -> float:
    """
    相当外気温度[℃]を計算する

    :param theta_e: 外気温度[degC]
    :param a_surf:  外気側表面日射吸収率[-]
    :param j_surf:  外気側表面に入射する日射量[W/m2]
    :param h_out:   室外側総合熱伝達率[W/(m2・K)]
    :return:        相当外気温度[℃]
    """
    return theta_e + (a_surf * j_surf) / h_out


def get_u_s_dash(angle: float, C_2: float, h_cv: float, h_rv: float) -> float:
    """
    通気層から室内までの熱貫流率（温度、風速依存の熱伝達率を使用した熱貫流率)U'_s[W/(m2・K)]を計算する

    :param angle:   通気層の傾斜角[°]
    :param C_2:     室内側部材の熱コンダクタンス[W/(m2・K)]
    :param h_cv:    通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:    通気層の放射熱伝達率[W/(m2・K)]
    :return:        通気層から室内までの熱貫流率（温度、風速依存の熱伝達率を使用）[W/(m2・K)]
    """

    # 表面熱伝達抵抗を設定
    # 室外側表面熱伝達抵抗は、省エネ基準に規定の値とする（通気層の角度（angle）が90°の場合は外壁の値（0.11）、90°以外の場合は屋根の値（0.09）を使用）
    if angle == 90.0:
        r_s_e = 0.11  # 室外側表面熱伝達抵抗, (m2・K)/W
        r_s_r = 0.11  # 室内側表面熱伝達抵抗, (m2・K)/W
    else:
        r_s_e = 0.09  # 室外側表面熱伝達抵抗, (m2・K)/W
        r_s_r = 0.09  # 室内側表面熱伝達抵抗, (m2・K)/W

    # 省エネ基準でのU値を計算
    u_s = 1 / (r_s_e + 1 / C_2 + r_s_r)

    if abs(h_rv + h_cv) < 0.001:
        u_s_dash = np.nan
    else:
        u_s_dash = 1.0 / (1.0/u_s - r_s_e + 1.0/(h_rv + h_cv))

    return u_s_dash


def get_u_o(C_1: float, h_cv: float, h_rv: float) -> float:
    """
    室外側から通気層までの熱貫流率を計算する

    :param C_1:     外気側部材の熱コンダクタンス[W/(m2・K)]
    :param h_cv:    通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:    通気層の放射熱伝達率[W/(m2・K)]
    :return:        室外側から通気層までの熱貫流率[W/(m2・K)]
    """

    # 固定値を設定
    h_out = global_number.get_h_out()  # 室外側総合熱伝達率, W/(m2・K)

    return 1.0 / (1.0/h_out + 1.0/C_1 + 1.0/(h_cv + h_rv))


def get_u_i(C_2: float, h_cv: float, h_rv: float) -> float:
    """
    室内側から通気層までの熱貫流率を計算する

    :param C_2:     室内側部材の熱コンダクタンス[W/(m2・K)]
    :param h_cv:    通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:    通気層の放射熱伝達率[W/(m2・K)]
    :return:        室内側から通気層までの熱貫流率[W/(m2・K)]
    """

    # 固定値を設定
    h_in = global_number.get_h_in()  # 室内側総合熱伝達率, W/(m2・K)

    return 1.0 / (1.0/(h_cv + h_rv) + 1.0/C_2 + 1.0/h_in)


def get_r_o(C_1: float) -> float:
    """
    室外側から通気層表面までの熱抵抗を計算する

    :param C_1:     外気側部材の熱コンダクタンス[W/(m2・K)]
    :return:        室外側から通気層表面までの熱抵抗[(m2・K)/W]
    """

    # 固定値を設定
    h_out = global_number.get_h_out()  # 室外側総合熱伝達率, W/(m2・K)

    return 1.0 / h_out + 1.0 / C_1


def get_r_i(C_2: float) -> float:
    """
    室内側から通気層表面までの熱抵抗を計算する

    :param C_2:     室内側部材の熱コンダクタンス[W/(m2・K)]
    :return:        室内側から通気層表面までの熱抵抗[(m2・K)/W]
    """

    # 固定値を設定
    h_in = global_number.get_h_in()  # 室内側総合熱伝達率, W/(m2・K)

    return 1.0 / h_in + 1.0 / C_2


def get_theata_as_e(theta_as_ave: float, theta_1_surf: float, h_cv: float, h_rv: float) -> float:
    """
    通気層の等価温度[℃]を計算する

    :param theta_as_ave:    通気層の平均温度[degC]
    :param theta_1_surf:    通気層内の面1の温度[degC]
    :param h_cv:            通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:            通気層の放射熱伝達率[W/(m2・K)]
    :return:                通気層の等価温度[degC]
    """

    if abs(h_rv + h_cv) < 0.001:
        theata_as_e = np.nan
    else:
        theata_as_e = (theta_as_ave * h_cv + theta_1_surf * h_rv) / (h_cv + h_rv)

    return theata_as_e

def get_heat_flow_room_side(angle: float, C_2: float, h_cv: float, h_rv: float,
                            theta_as_e: float, theta_r: float) -> float:
    """
    壁体の室内側表面熱流[W/m2]を計算する

    :param angle:   通気層の傾斜角[°]
    :param C_2:     室内側部材の熱コンダクタンス[W/(m2・K)]
    :param h_cv:    通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:    通気層の放射熱伝達率[W/(m2・K)]
    :param theta_as_e:  通気層の平均温度[℃]
    :param theta_r:     室内温度[℃]
    :return:            通気層の等価温度[℃]
    """

    # 通気層から室内までの熱貫流率（温度、風速依存の熱伝達率を使用した熱貫流率)[W / (m2・K)]を取得
    u_s_dash = get_u_s_dash(angle, C_2, h_cv, h_rv)

    return u_s_dash * (theta_as_e - theta_r)
