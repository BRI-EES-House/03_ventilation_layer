
def get_c_air(T: float) -> float:
    """
    Returns:
        空気の定圧比熱, J/(kg・K)
    """

    return 1006.0


def get_rho_air(t: float) -> float:
    """
    Returns:
        空気の密度, kg/m3
    """

    return 1.293 / (1.0 + t / get_abs_temp())


def get_g() -> float:
    """
    Returns:
        重力加速度, m/s**2
    """

    return 9.8


def get_lambda_air(t: float) -> float:
    """
    Returns:
        空気の熱伝導率, W/(m・K)
    """
    return 0.0241 + 7.7e-5 * t


def get_beta_air(t: float) -> float:
    """
    Returns:
        空気の体膨張率, 1/K
    """
    return 1.0 / (t + get_abs_temp())


def get_mu_air(t: float) -> float:
    """
    Returns:
        空気の粘性率, Pa・s
    """
    return (0.0074237 / (t + 390.15)) * ((t + get_abs_temp()) / 293.15) ** 1.5


def get_new_air(t: float) -> float:
    """
    :param t: 空気温度
    :return: 空気の動粘性係数[m2/s]
    """
    return get_mu_air(t) / get_rho_air(t)


def get_a_air(t: float) -> float:
    """
    :param t: 空気温度
    :return: 空気の熱拡散率
    """
    return get_lambda_air(t) / get_c_air(t) / get_rho_air(t)


def get_pr_air(t: float) -> float:
    """
    :param t: 空気温度
    :return: プラントル数
    """
    return get_new_air(t) / get_a_air(t)


def get_gr_air(tw: float, tf: float, d: float) -> float:
    """
    :param tw: 表面温度, C
    :param tf: 流体温度, C
    :param d: 代表長さ, m
    :return: グラスホフ数
    """
    # 膜温度の計算
    t_ave = (tw + tf) / 2.0
    return get_g() * get_beta_air(tf) * abs(tw - tf) * d ** 3.0 / get_new_air(t_ave) ** 2.0


def get_sgm() -> float:
    """
    Returns:
        ステファンボルツマン定数
    """
    return 5.67e-8


def get_abs_temp() -> float:
    """
    Returns:
        絶対温度, K
    """
    return 273.15


def get_h_out() -> float:
    """
    Returns:
        室外側総合熱伝達率, W/(m2・K)
    """
    return 25.0


def get_h_in() -> float:
    """
    Returns:
        室内側総合熱伝達率, W/(m2・K)
    """
    return 9.0


def get_surface_albedo() -> float:
    """
    :return: 地面の日射反射率（アルベド）, -
    """
    return 0.2
