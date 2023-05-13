import math
from typing import Optional

from ventilation_layer import global_number as gn


def get_e(eps1: float, eps2: float, l_d: Optional[float] = None, l_s: Optional[float] = None, method: Optional[str] = "parallel") -> float:

    if method == "parallel":
        return _get_e_parallel(eps1=eps1, eps2=eps2)
    elif method == "two_dimension":
        return _get_e_two_dimension(eps1=eps1, eps2=eps2, l_d=l_d, l_s=l_s)
    else:
        raise ValueError()
        

def _get_e_parallel(eps1: float, eps2: float) -> float:
    """有効放射率の計算（無限の平行面の場合）

    Args:
        eps1: 面1の放射率, -
        eps2: 面2の放射率, -
    Returns:
        有効放射率, -
    """

    return 1 / ((1 / eps1) + (1 / eps2) - 1)


def _get_e_two_dimension(eps1: float, eps2: float, l_d: float, l_s: float) -> float:
    """有効放射率の計算（二次元空間の場合）

    Args:
        eps1: 面1の放射率, -
        eps2: 面2の放射率, -
        l_d: 通気層の厚さ, m
        l_s: 通気胴縁または垂木の間隔, m
    Returns:
        有効放射率, -
    """

    return 1.0 / (1.0 / eps1 + 1.0 / eps2 - 2.0 + 1.0 / (1.0 / 2.0 * (1.0 + math.sqrt(1.0 + l_d**2.0 / l_s**2.0) - l_d / l_s)))


def get_h_rv(eps_eff: float, calc_mode: Optional[str] = "detailed", theta_1: Optional[float] = None, theta_2: Optional[float] = None) -> float:
    """計算モードに応じた放射熱伝達率を計算する

    Args:
        calc_mode: 計算モード
        theta_1: 通気層に面する面1の表面温度, degrees
        theta_2: 通気層に面する面2の表面温度, degrees
        eps_eff: 有効放射率, -
    Returns:
        放射熱伝達率, W/(m2・K)
    """

    if calc_mode == "detailed":
        h_rv = _get_h_rv_detailed(theta_1, theta_2, eps_eff)
    elif calc_mode == "simplified_winter":
        h_rv = _get_h_rv_simplified_winter(eps_eff)
    elif calc_mode == "simplified_summer":
        h_rv = _get_h_rv_simplified_summer(eps_eff)
    elif calc_mode == "simplified_all_season":
        h_rv = _get_h_rv_simplified_all_season(eps_eff)
    elif calc_mode == "simplified_zero":
        h_rv = 0.0
    else:
        raise ValueError("指定された計算モードは対象外です")

    return h_rv


def _get_h_rv_simplified_winter(eps_eff: float) -> float:
    """放射熱伝達率[W/(m2・K)]の計算（簡易計算、冬期条件）

    Args:
        effective_emissivity: 有効放射率, -
    Returns:
        放射熱伝達率, W/(m2・K)
    """

    return 5.054 * eps_eff


def _get_h_rv_simplified_summer(eps_eff: float) -> float:
    """放射熱伝達率[W/(m2・K)]の計算（簡易計算、夏期条件）

    Args:
        eps_eff: 有効放射率, -
    Returns:
        放射熱伝達率, W/(m2・K)
    """

    return 6.615 * eps_eff


def _get_h_rv_simplified_all_season(eps_eff: float) -> float:
    """放射熱伝達率[W/(m2・K)]の計算（簡易計算、通年）

    Args:
        eps_eff: 有効放射率, -
    Returns:
        放射熱伝達率, W/(m2・K)
    """

    return 5.88 * eps_eff


def _get_h_rv_detailed(theta_1: float, theta_2: float, eps_eff: float) -> float:
    """放射熱伝達率[W/(m2・K)]の計算（詳細計算）

    Args:
        theta_1: 通気層に面する面1の表面温度, degrees
        theta_2: 通気層に面する面2の表面温度, degrees
        eps_eff: 有効放射率, -
    Returns:
        放射熱伝達率, W/(m2・K)
    """

    t_m = (theta_1 + gn.get_abs_temp() + theta_2 + gn.get_abs_temp()) / 2
    
    h_rv = 4 * gn.get_sgm() * eps_eff * (t_m ** 3)
    
    return h_rv


def get_h_cv(
        calc_mode: str,
        v_a: float,
        theta_1: Optional[float] = None,
        theta_2: Optional[float] = None,
        angle: Optional[float] = None,
        l_h: Optional[float] = None,
        l_d: Optional[float] = None
    ) -> float:
    """計算モードに応じた対流熱伝達率を計算する

    Args:
        calc_mode: 計算モード
        v_a: 通気層の平均風速, m/s
        theta_1: 通気層に面する面1の表面温度, degrees
        theta_2: 通気層に面する面2の表面温度, degrees
        angle: 通気層の傾斜角, degrees
        l_h: 通気層の長さ, m
        l_d: 通気層の厚さ, m
    Returns:
        対流熱伝達率, W/(m2・K)
    """
    if calc_mode == "detailed":
        return _get_h_cv_detailed(v_a, theta_1, theta_2, angle, l_h, l_d)
    elif calc_mode == "simplified_winter":
        return _get_h_cv_simplified_winter(v_a)
    elif calc_mode == "simplified_summer":
        return _get_h_cv_simplified_summer(v_a)
    elif calc_mode == "simplified_all_season":
        return _get_h_cv_simplified_all_season(v_a)
    else:
        raise ValueError("指定された計算モードは対象外です")


def _get_h_cv_simplified_winter(v_a: float) -> float:
    """対流熱伝達率[W/(m2・K)]の計算（簡易計算、冬期条件）

    Args:
        v_a: 通気層の平均風速, m/s
    Returns:
        対流熱伝達率, W/(m2・K)
    """

    return 4.077 * v_a + 2.302


def _get_h_cv_simplified_summer(v_a: float) -> float:
    """対流熱伝達率[W/(m2・K)]の計算（簡易計算、夏期条件）

    Args:
        v_a: 通気層の平均風速, m/s
    Returns:
        対流熱伝達率, W/(m2・K)
    """

    return 4.113 * v_a + 1.844


def _get_h_cv_simplified_all_season(v_a: float) -> float:
    """対流熱伝達率[W/(m2・K)]の計算（簡易計算、通年）

    Args:
        v_a: 通気層の平均風速, m/s
    Returns:
        対流熱伝達率, W/(m2・K)
    """

    return 4.096 * v_a + 2.06


def _get_h_cv_detailed(v_a: float, theta_1: float, theta_2: float, angle: float, l_h: float, l_d: float) -> float:
    """対流熱伝達率[W/(m2・K)]の計算（詳細計算）

    Args:
        v_a: 通気層の平均風速, m/s
        theta_1: 通気層に面する面1の表面温度, degC
        theta_2: 通気層に面する面2の表面温度, degC
        angle: 通気層の傾斜角, degree
        l_h: 通気層の長さ, m
        l_d: 通気層の厚さ, m
    Returns
        対流熱伝達率, W/(m2・K)
    """

    theta_ave = (theta_1 + theta_2) / 2.0

    if theta_1 == theta_2:
        # 両表面の温度（theta_1とtheta_2）が同じ値のときはh_c = 0.0とする
        h_cv = 0.0
    else:
        # 風速係数, (W/s)(m3 K)
        c_v = 4

        # ヌセルト数を計算
        nusselt_number = _get_n_u(theta_1, theta_2, angle, l_h, l_d)

        # 密閉空気層の自然対流熱伝達率を計算
        h_base = nusselt_number * gn.get_lambda_air(theta_ave) / l_d

        # 通気層の対流熱伝達率の計算
        h_cv = 2 * h_base + c_v * v_a

    return h_cv


def _get_n_u(theta_1: float, theta_2: float, angle: float, l_h: float, l_d: float) -> float:
    """ヌセルト数の計算

    Args:
        theta_1: 通気層に面する面1の表面温度, degrees
        theta_2: 通気層に面する面2の表面温度, degrees
        angle: 通気層の傾斜角, degrees
        l_h: 通気層の長さ, m
        l_d: 通気層の厚さ, m
    Returns:
        ヌセルト数
    """

    # 表面温度の平均値
    theta_ave = (theta_1 + theta_2) / 2.0

    # レーリー数の計算
    r_a = (gn.get_g() * gn.get_beta_air(theta_ave) * abs(theta_1 - theta_2) * (l_d ** 3) * (gn.get_rho_air(theta_ave) ** 2) * gn.get_c_air()) / (gn.get_mu_air(theta_ave) * gn.get_lambda_air(theta_ave))

    # ヌセルト数の計算
    nusselt_number = 0
    nu_ct = (1.0 + ((0.104 * r_a ** 0.293) / (1.0 + (6310.0 / r_a) ** 1.36)) ** 3) ** (1 / 3)
    nu_u1 = 0.242 * (r_a * l_d / l_h) ** 0.273
    nu_ut = 0.0605 * r_a ** (1 / 3)

    # 傾斜角が0°（水平）のとき
    if angle == 0.0:
        if r_a > 5830.0:
            nusselt_number = 1.44 * (1.0 - 1708.0/r_a) + (r_a/5830.0) ** (1/3)
        elif r_a > 1708.0:
            nusselt_number = 1.0 + 1.44 * (1.0 - 1708.0/r_a)
        else:
            nusselt_number = 1.0

    # 傾斜角が90°（鉛直）のとき
    elif angle == 90.0:
        nusselt_number = max(nu_ct, nu_u1, nu_ut)

    # 傾斜角が0°<γ≤60°のとき
    elif 0.0 < angle <= 60.0:
        buff = r_a * math.cos(math.radians(angle))
        if buff >= 5830.0:
            nusselt_number = 1.44 * (1.0 - 1708.0/buff) * (1.0 - (1708.0 * (math.sin(1.8 * math.radians(angle)) ** 1.6))/buff) + (buff/5830.0) ** (1/3)
        elif 1708.0 <= buff < 5830.0:
            nusselt_number = 1.44 * (1.0 - 1708.0/buff) * (1.0 - (1708.0 * (math.sin(1.8 * math.radians(angle)) ** 1.6))/buff)
        elif buff < 1708.0:
            nusselt_number = 1.0
        else:
            raise Exception()

    # 傾斜角が60°<γ<90°のとき
    elif 60.0 < angle < 90.0:
        buff_g = 0.5 / (1.0 + (r_a/3165.0) ** 20.6) ** 0.1
        nu_60_1 = (1.0 + ((0.0936 * r_a ** 0.314) ** 7) / (1.0 + buff_g)) ** (1 / 7)
        nu_60_2 = (0.1044 + 0.1759 * l_d/l_h) * r_a ** 0.283
        nu_60 = max(nu_60_1, nu_60_2)
        nu_v = max(nu_ct, nu_u1, nu_ut)
        nusselt_number = nu_60 * (90.0 - angle)/30.0 + nu_v * (angle - 60.0)/30.0
       
    else:
        raise ValueError("指定された傾斜角は計算対象外です")

    return nusselt_number


def get_h_cv_Jurges(v_a: float) -> float:
    """Get convective heat transfer coefficient base on Jurges equation.

    Args:
        v_a: air speed, m/s

    Returns:
        forced convective heat transfer coefficient, W/m2K
    """

    if v_a <= 4.9:
        return 5.6 + 3.9 * v_a
    else:
        return 7.2 * v_a**0.78

