import ventilation_wall
import global_number


def overall_heat_transfer_coefficient(theta_e: float, theta_r: float, a_surf: float, j_surf: float, C_2: float,
                                      angle: float, theta_as_ave: float, h_cv: float, h_rv: float) -> float:
    """
    通気層を有する壁体の熱貫流率[W/(m2・K)]を計算する

    :param theta_e: 外気温度[℃]
    :param theta_r: 室内温度[℃]
    :param a_surf:  外気側表面日射吸収率[-]
    :param j_surf:  外気側表面に入射する日射量[W/m2]
    :param C_2:     室内側部材の熱コンダクタンス[W/(m2・K)]
    :param angle:   通気層の傾斜角[°]
    :param theta_as_ave: 通気層の平均温度[℃]
    :param h_cv:    通気層の対流熱伝達率[W/(m2・K)]
    :param h_rv:    通気層の放射熱伝達率[W/(m2・K)]
    :return:        通気層の等価温度[℃]
    """

    # 固定値を設定
    h_out = global_number.get_h_out()   # 室外側総合熱伝達率, W/(m2・K)

    # 表面熱伝達抵抗を設定
    # TODO: 通気層の角度（angle）が90°の場合は外壁の値（0.11）、90°以外の場合は屋根の値（0.09）を使用する
    r_s_e = 0.11    # 室外側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    r_s_r = 0.11    # 室内側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値

    # 相当外気温度を計算
    theta_SAT = get_theta_SAT(theta_e, a_surf, j_surf, h_out)

    # 温度、風速依存の熱伝達率を使用したU値に修正
    u_s_dash = get_u_s_dash(r_s_e, r_s_r, C_2, h_cv, h_rv)

    # 通気層を有する壁体の相当熱貫流率を求めるための補正係数を計算
    k_e = get_weight_factor_of_u_s_dash(theta_as_ave, theta_r, theta_SAT)

    # 通気層を有する壁体の熱貫流率(W/(m2・K))を計算
    u_e = u_s_dash * k_e

    return u_e


# 通気層を有する壁体の日射熱取得率(-)の計算
def solar_heat_gain_coefficient(parm: ventilation_wall.Parameters, theta_as_ave: float, h_cv: float, h_rv: float):

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


# 相当外気温度を計算
def get_theta_SAT(theta_e: float, a_surf: float, J_surf: float, h_out: float):
    return theta_e + (a_surf * J_surf) / h_out


# 温度、風速依存の熱伝達率を使用したU値を計算
def get_u_s_dash(r_s_e: float, r_s_r: float, C_2: float, h_cv: float, h_rv: float):
    # 省エネ基準でのU値を計算
    u_s = 1 / (r_s_e + 1 / C_2 + r_s_r)
    return 1/(1/u_s - r_s_e + 1/(h_rv + h_cv))
