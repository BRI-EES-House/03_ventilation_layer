import ventilation_wall


# 通気層を有する壁体の熱貫流率の計算（W/(m2・K)）
def overall_heat_transfer_coefficient(parm: ventilation_wall.Parameters, theta_as_ave: float, h_cv: float, h_rv: float):

    # 固定値を設定
    # Note:固定値はグローバル変数にするのがよいか
    h_out = 25.0    # 室外側総合熱伝達率, W/(m2・K)
    r_s_e = 0.11    # 室外側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    r_s_r = 0.11    # 室内側表面熱伝達抵抗, (m2・K)/W, 省エネ基準の規定値
    # Note: 省エネ基準の規定による表面熱伝達抵抗は、屋根、外壁で値が異なる。どのように判定するかが課題。

    # 相当外気温度を計算
    theta_SAT = parm.theta_e + (parm.a_surf * parm.J_surf) / h_out


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
def solar_heat_gain_coefficient(parm: ventilation_wall.Parameters, theta_as_ave: float, h_cv: float, h_rv: float):

    # 固定値を設定
    # Note:固定値はグローバル変数にするのがよいか
    h_in = 9.0  # 室内側総合熱伝達率, W/(m2・K)



    # 通気層を有する壁体の日射熱取得率(-)を計算
    eta_e = h_in * (matrix_temp[3][0] - parm.theta_r) / parm.J_surf

    return eta_e

