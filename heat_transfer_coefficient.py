import math


# 有効放射率の計算（無限の平行面の場合）
def effective_emissivity_parallel(emissivity_1, emissivity_2):
    effective_emissivity = 1 / ((1 / emissivity_1) + (1 / emissivity_2) - 1)
    return effective_emissivity


# 有効放射率の計算（二次元空間の場合）
def effective_emissivity_two_dimension(emissivity_1, emissivity_2, l_d, l_s):
    effective_emissivity = 1.0 / (1.0 / emissivity_1 + 1.0 / emissivity_2 - 2.0 + 1.0 / (1.0 / 2.0 * (1.0 + math.sqrt(1.0 + l_d ** 2.0 / l_s ** 2.0) - l_d / l_s)))
    return effective_emissivity


# 放射熱伝達率[W/(m2・K)]の計算
def radiative_heat_transfer_coefficient(theta_1, theta_2, effective_emissivity):
    t_m = (theta_1 + 273.15 + theta_2 + 273.15) / 2
    sigma = 5.67 * (10 ** (-8))  # ステファン・ボルツマン定数（W/(m^2・K^4)）
    h_r = 4 * sigma * effective_emissivity * (t_m ** 3)
    return h_r


# 対流熱伝達率の計算[W/(m2・K)]
def convective_heat_transfer_coefficient(v_a, theta_1, theta_2, angle, l_h, l_d):

    # 固定値を設定
    lambda_a = 0.026  # 空気の熱伝導率, W/(m・K)

    if theta_1 == theta_2:
        # 両表面の温度（theta_1とtheta_2）が同じ値のときはh_c = 1とする
        h_c = 1
    else:
        # ヌセルト数を計算
        nusselt_number = get_nusselt_number(theta_1, theta_2, angle, l_h, l_d)

        # 密閉空気層の自然対流熱伝達率を計算
        h_base = nusselt_number * lambda_a / l_d

        # 通気層の対流熱伝達率の計算
        h_c = 2 * h_base * 4 * v_a

    return h_c


# ヌセルト数の計算
def get_nusselt_number(theta_1, theta_2, angle, l_h, l_d):

    # 固定値を設定
    g = 9.8                     # 重力加速度,m/s^2
    beta_a = 3.7 * (10 ** -3)    # 空気の体膨張率, 1/K
    lambda_a = 0.026            # 空気の熱伝導率, W/(m・K)
    c_a = 1006                 # 空気の定圧比熱, J/(kg・K)
    rho_a = 1.2                 # 空気の密度, kg/m3
    mu_a = 1.8 * (10 ** -5)     # 空気の粘性率, Pa・s

    # レーリー数の計算
    rayleigh_number = (g * beta_a * abs((theta_1 + 273.15) - (theta_2 + 273.15)) * (l_d ** 3) * (rho_a ** 2) * c_a) / (mu_a * lambda_a)

    # ヌセルト数の計算
    nusselt_number = 0
    nu_ct = (1 + ((0.104 * rayleigh_number ** 0.293) / (1 + (6310 / rayleigh_number) ** 1.36)) ** 3) ** (1 / 3)
    nu_u1 = 0.242 * (rayleigh_number * l_d / l_h) ** 0.273
    nu_ut = 0.0605 * rayleigh_number ** (1 / 3)

    # 傾斜角が0°（水平）のとき
    if angle == 0:
        if rayleigh_number > 5830:
            nusselt_number = 1.44 * (1 - 1708/rayleigh_number) + (rayleigh_number/5830) ** (1/3)
        elif 1708 < rayleigh_number <= 5830:
            nusselt_number = 1 + 1.44 * (1 - 1708/rayleigh_number)
        elif rayleigh_number <= 1708:
            nusselt_number = 1

    # 傾斜角が90°（鉛直）のとき
    elif angle == 90:
        # nu_ct = (1 + ((0.104 * rayleigh_number ** 0.293) / (1 + (6310 / rayleigh_number) ** 1.36)) ** 3) ** (1 / 3)
        # nu_u1 = 0.242 * (rayleigh_number * l_d / l_h) ** 0.273
        # nu_ut = 0.0605 * rayleigh_number ** (1/3)
        nusselt_number = max(nu_ct, nu_u1, nu_ut)

    # 傾斜角が0°<γ≤60°のとき
    elif 0 < angle <= 60:
        buff = rayleigh_number * math.cos(angle)
        if buff >= 5830:
            nusselt_number = 1.44 * (1 - 1708/buff) * (1 - (1708 * (math.sin(1.8 * angle) ** 1.6))/buff) + (buff/5830) ** (1/3)
        elif 1708 <= buff < 5830:
            nusselt_number = 1.44 * (1 - 1708 / buff) * (1 - (1708 * (math.sin(1.8 * angle) ** 1.6)) / buff)
        elif buff < 1708:
            nusselt_number = 1.0

    # 傾斜角が60°<γ<90°のとき
    elif 60 < angle < 90:
        buff_g = 0.5 / (1 + (rayleigh_number/3165) ** 20.6) ** 0.1
        nu_60_1 = (1 + ((0.0936 * rayleigh_number ** 0.314) ** 7) / (1 + buff_g)) ** (1 / 7)
        nu_60_2 = (0.1044 + 0.1759 * l_d/l_h) * rayleigh_number ** 0.283
        nu_60 = max(nu_60_1, nu_60_2)
        nu_v = max(nu_ct, nu_u1, nu_ut)
        nusselt_number = nu_60 * (90 - angle)/30 + nu_v * (angle - 60)/30
       
    else:
        raise ValueError("指定された傾斜角は計算対象外です")

    return nusselt_number
