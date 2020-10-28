
def get_c_air() -> float:
    """
    Returns:
        空気の定圧比熱, J/(kg・K)
    """

    return 1006.0


def get_rho_air() -> float:
    """
    Returns:
        空気の密度, kg/m3
    """

    return 1.2


def get_g() -> float:
    """
    Returns:
        重力加速度, m/s**2
    """

    return 9.8


def get_lambda_air() -> float:
    """
    Returns:
        空気の熱伝導率, W/(m・K)
    """
    return 0.026


def get_beta_air() -> float:
    """
    Returns:
        空気の体膨張率, 1/K
    """
    return 3.7e-3


def get_mu_air() -> float:
    """
    Returns:
        空気の粘性率, Pa・s
    """
    return 1.8e-5


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
