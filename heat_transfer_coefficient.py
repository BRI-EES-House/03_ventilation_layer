import math
from global_number import get_abs_temp, get_sgm, get_g, get_lambda_air, get_beta_air, get_mu_air, get_pr_air, get_c_air, get_rho_air


def effective_emissivity_parallel(emissivity_1: float, emissivity_2: float) -> float:
    """
    有効放射率の計算（無限の平行面の場合）

    :param emissivity_1:    面1の放射率, -
    :param emissivity_2:    面2の放射率, -
    :return:                有効放射率, -
    """

    effective_emissivity = 1 / ((1 / emissivity_1) + (1 / emissivity_2) - 1)
    return effective_emissivity


def effective_emissivity_two_dimension(emissivity_1: float, emissivity_2: float, l_d: float, l_s: float) -> float:
    """
    有効放射率の計算（二次元空間の場合）

    :param emissivity_1:    面1の放射率, -
    :param emissivity_2:    面2の放射率, -
    :param l_d:             通気層の厚さ, m
    :param l_s:             通気胴縁または垂木の間隔, m
    :return:                有効放射率, -
    """
    effective_emissivity = 1.0 / (1.0 / emissivity_1 + 1.0 / emissivity_2 - 2.0 + 1.0 /
                                  (1.0 / 2.0 * (1.0 + math.sqrt(1.0 + l_d ** 2.0 / l_s ** 2.0) - l_d / l_s)))
    return effective_emissivity


def get_radiative_heat_transfer_coefficient(calc_mode: str, theta_1: float, theta_2: float, effective_emissivity: float) -> float:
    """
    計算モードに応じた放射熱伝達率を計算する

    :param calc_mode:   計算モード
    :param theta_1:     通気層に面する面1の表面温度, degC
    :param theta_2:     通気層に面する面2の表面温度, degC
    :param effective_emissivity: 有効放射率, -
    :return:            放射熱伝達率, W/(m2・K)
    """
    if calc_mode == "detailed":
        h_rv = radiative_heat_transfer_coefficient_detailed(theta_1, theta_2, effective_emissivity)
    elif calc_mode == "simplified_winter":
        h_rv = radiative_heat_transfer_coefficient_simplified_winter(effective_emissivity)
    elif calc_mode == "simplified_summer":
        h_rv = radiative_heat_transfer_coefficient_simplified_summer(effective_emissivity)
    elif calc_mode == "simplified_all_season":
        h_rv = radiative_heat_transfer_coefficient_simplified_all_season(effective_emissivity)
    elif calc_mode == "simplified_zero":
        h_rv = 0.0
    else:
        raise ValueError("指定された計算モードは対象外です")

    return h_rv


def radiative_heat_transfer_coefficient_simplified_winter(effective_emissivity: float) -> float:
    """
    放射熱伝達率[W/(m2・K)]の計算（簡易計算、冬期条件）

    :param effective_emissivity: 有効放射率, -
    :return:                     放射熱伝達率, W/(m2・K)
    """
    return 5.054 * effective_emissivity


def radiative_heat_transfer_coefficient_simplified_summer(effective_emissivity: float) -> float:
    """
    放射熱伝達率[W/(m2・K)]の計算（簡易計算、夏期条件）

    :param effective_emissivity: 有効放射率, -
    :return:                     放射熱伝達率, W/(m2・K)
    """
    return 6.615 * effective_emissivity


def radiative_heat_transfer_coefficient_simplified_all_season(effective_emissivity: float) -> float:
    """
    放射熱伝達率[W/(m2・K)]の計算（簡易計算、通年）

    :param effective_emissivity: 有効放射率, -
    :return:                     放射熱伝達率, W/(m2・K)
    """
    return 5.88 * effective_emissivity


def radiative_heat_transfer_coefficient_detailed(theta_1: float, theta_2: float, effective_emissivity: float) -> float:
    """
    放射熱伝達率[W/(m2・K)]の計算（詳細計算）

    :param theta_1:     通気層に面する面1の表面温度, degC
    :param theta_2:     通気層に面する面2の表面温度, degC
    :param effective_emissivity: 有効放射率, -
    :return:            放射熱伝達率, W/(m2・K)
    """
    t_m = (theta_1 + get_abs_temp() + theta_2 + get_abs_temp()) / 2
    h_rv = 4 * get_sgm() * effective_emissivity * (t_m ** 3)
    return h_rv


def get_convective_heat_transfer_coefficient(calc_mode: str, v_a: float, theta_1: float, theta_2: float, angle: float, l_h: float, l_d: float) -> float:
    """
    計算モードに応じた対流熱伝達率を計算する

    :param calc_mode:   計算モード
    :param v_a:         通気層の平均風速, m/s
    :param theta_1:     通気層に面する面1の表面温度, degC
    :param theta_2:     通気層に面する面2の表面温度, degC
    :param angle:       通気層の傾斜角, degree
    :param l_h:         通気層の長さ, m
    :param l_d:         通気層の厚さ, m
    :return:            対流熱伝達率, W/(m2・K)
    """
    if calc_mode == "detailed":
        h_cv = convective_heat_transfer_coefficient_detailed(v_a, theta_1, theta_2, angle, l_h, l_d)
    elif calc_mode == "simplified_winter":
        h_cv = convective_heat_transfer_coefficient_simplified_winter(v_a)
    elif calc_mode == "simplified_summer":
        h_cv = convective_heat_transfer_coefficient_simplified_summer(v_a)
    elif calc_mode == "simplified_all_season":
        h_cv = convective_heat_transfer_coefficient_simplified_all_season(v_a)
    else:
        raise ValueError("指定された計算モードは対象外です")

    return h_cv


def convective_heat_transfer_coefficient_simplified_winter(v_a: float) -> float:
    """
    対流熱伝達率[W/(m2・K)]の計算（簡易計算、冬期条件）

    :param v_a: 通気層の平均風速, m/s
    :return:    対流熱伝達率, W/(m2・K)
    """
    return 4.077 * v_a + 2.302


def convective_heat_transfer_coefficient_simplified_summer(v_a: float) -> float:
    """
    対流熱伝達率[W/(m2・K)]の計算（簡易計算、夏期条件）

    :param v_a: 通気層の平均風速, m/s
    :return:    対流熱伝達率, W/(m2・K)
    """
    return 4.113 * v_a + 1.844


def convective_heat_transfer_coefficient_simplified_all_season(v_a: float) -> float:
    """
    対流熱伝達率[W/(m2・K)]の計算（簡易計算、通年）

    :param v_a: 通気層の平均風速, m/s
    :return:    対流熱伝達率, W/(m2・K)
    """
    return 4.096 * v_a + 2.06



    theta_ave = (theta_1 + theta_2) / 2.0

    if theta_1 == theta_2:
        # 両表面の温度（theta_1とtheta_2）が同じ値のときはh_c = 0.0とする
        h_c = 0.0
    else:
        # ヌセルト数を計算
        nusselt_number = get_nusselt_number(theta_1, theta_2, angle, l_h, l_d)

        # 密閉空気層の自然対流熱伝達率を計算
        h_base = nusselt_number * get_lambda_air(theta_ave) / l_d

        # 通気層の対流熱伝達率の計算
        if v_a > 0.0:
            h_c = 2 * h_base + 4 * v_a
        else:
            h_c = h_base

    return h_c


# ヌセルト数の計算
def get_nusselt_number(theta_1: float, theta_2: float, angle: float, l_h: float, l_d: float) -> float:
    """
    ヌセルト数の計算

    :param theta_1:     通気層に面する面1の表面温度, degC
    :param theta_2:     通気層に面する面2の表面温度, degC
    :param angle:       通気層の傾斜角, degree
    :param l_h:         通気層の長さ, m
    :param l_d:         通気層の厚さ, m
    :return:            ヌセルト数
    """

    # 表面温度の平均値
    theta_ave = (theta_1 + theta_2) / 2.0

    # プラントル数の計算
    pr = get_pr_air(theta_ave)

    # レーリー数の計算
    rayleigh_number = (get_g() * get_beta_air(theta_ave) * abs(theta_1 - theta_2) * (l_d ** 3) * (get_rho_air(theta_ave) ** 2) * get_c_air(theta_ave)) / (get_mu_air(theta_ave) * get_lambda_air(theta_ave))

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
        buff = rayleigh_number * math.cos(math.radians(angle))
        if buff >= 5830:
            nusselt_number = 1.44 * (1 - 1708/buff) * (1 - (1708 * (math.sin(1.8 * math.radians(angle)) ** 1.6))/buff) + (buff/5830) ** (1/3)
        elif 1708 <= buff < 5830:
            nusselt_number = 1.44 * (1 - 1708 / buff) * (1 - (1708 * (math.sin(1.8 * math.radians(angle)) ** 1.6)) / buff)
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
