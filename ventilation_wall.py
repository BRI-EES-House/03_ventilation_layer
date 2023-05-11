import math
from scipy import optimize
import numpy as np
from dataclasses import dataclass

from ventilation_layer import heat_transfer_coefficient as htc
from ventilation_layer import global_number as gn


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

    # the heat flow from the outside(inclued solar irradiance) to the exterior surface, W/m2
    q_flow_out: float

    # the heat flow from the exterior surface to the surface of exterior side facing the ventilation layer, W/m2
    q_flow_em: float

    # the exhausted heat flow from the ventilation layer, W/m2
    q_flow_exhaust: float

    # the heat flow to the inside, W/m2
    q_flow_in: float

    # 最適化が正常に終了したかどうか
    is_optimize_succeed: bool

    # 最適化の終了ステータス
    optimize_status: int

    # 最適化の終了メッセージ
    optimize_message: str


def _get_heat_balance(
        matrix_temp: np.zeros(5),
        theta_e: float,
        theta_r: float,
        j_surf: float,
        a_surf: float,
        c_1: float,
        c_2: float,
        l_h: float,
        l_w: float,
        l_d: float,
        angle: float,
        v_a: float,
        eps1: float,
        eps2: float,
        calc_mode_h_cv: str,
        calc_mode_h_rv: str,
        h_out: float,
        h_in: float
    ) -> np.zeros(5):
    """
    熱収支式を解く関数

    Args:
        matrix_temp: the temperature of the points below, degrees
            0: the temperature on the exterior surface
            1: the temperature on the exterior side surface facing the ventilation layer
            2: the temperature on the interior side surface facing the ventilation layer
            3: the temperature on the interior surface
            4: the temperature of the air in the ventilation layer
        theta_e: the outdoor temperature, degrees
        theta_r: the indoor temperature, degrees
        j_surf: the solar irradiance on the exterior surface, W/m2
        a_surf: the solar absorption ratio on the external surface, -
        c_1: the thermal conductance of the outside material, W/m2K
        c_2: the thermal conductance of the inside material, W/m2K
        l_h: the length of the ventilation layer, m
        l_w: the width of the ventilation layer, m
        l_d: the thickness of the ventilation layer, m
        angle: the angle of the ventilation layer, degrees
        v_a: the mean air velocity of the ventilation layer, m/s
        eps1: the emissivity of the surface 1 facing the ventilation layer
        eps2: the emissivity of the surface 2 facing the ventilation layer
        calc_mode_h_cv: 対流熱伝達率の計算モード
        calc_mode_h_rv: 放射熱伝達率の計算モード
        h_out: 室外側総合熱伝達率, W/(m2・K)
        h_in: 室内側総合熱伝達率, W/(m2・K)
    Returns:
        各層の熱収支, W/m2
    """

    # SAT temp, degrees
    theta_sat = theta_e + (a_surf * j_surf) / h_out

    # the temperature of the surfaces facing the ventilation layer, degrees
    theta_1 = matrix_temp[1]
    theta_2 = matrix_temp[2]
    theta_4 = matrix_temp[4]

    # the convective heat transfer coefficient, W/m2K
    h_cv = htc.get_h_cv(calc_mode_h_cv, v_a, theta_1, theta_2, angle, l_h, l_d)

    # the effective emissivity, -
    eps_eff = htc.get_e(eps1=eps1, eps2=eps2)

    # the radiative heat transfer coefficient, W/m2K
    h_rv = htc.get_h_rv(eps_eff, calc_mode_h_rv, theta_1, theta_2)

    # the ventilation air volume, m3/s
    v_vent = v_a * l_d * l_w

    if v_a > 0.0:
        beta = (2 * h_cv * l_w) / (gn.get_c_air() * gn.get_rho_air(theta_4) * v_vent)
        a41 = (1.0 + 1.0 / l_h * 1.0 / beta * (math.exp(-beta * l_h) - 1)) / 2
        a42 = (1.0 + 1.0 / l_h * 1.0 / beta * (math.exp(-beta * l_h) - 1)) / 2
        b4 = 1.0 / l_h * 1.0 / beta * (math.exp(-beta * l_h) - 1) * theta_e
    else:
        a41 = 0.5
        a42 = 0.5
        b4 = 0.0

    # the matrix equation
    # A \theta = B

    # A
    a = np.array([
        [h_out + c_1, -c_1, 0.0, 0.0, 0.0],
        [c_1, -(h_cv + h_rv + c_1), h_rv, 0.0, h_cv],
        [0.0, h_rv, -(h_cv + h_rv + c_2), c_2, h_cv],
        [0.0, 0.0, c_2, -(h_in + c_2), 0.0],
        [0.0, a41, a42, 0.0, -1.0]
    ])

    # B
    b = np.array([h_out * theta_sat, 0.0, 0.0, -h_in * theta_r, b4])

    # 熱収支を計算
    q_balance = np.matmul(a, matrix_temp) - b

    return q_balance


def get_wall_status_values(
        theta_e: float,
        theta_r: float,
        j_surf: float,
        a_surf: float,
        c_1: float,
        c_2: float,
        l_h: float,
        l_w: float,
        l_d: float,
        angle: float,
        v_a: float,
        eps_1: float,
        eps_2: float,
        calc_mode_h_cv: str,
        calc_mode_h_rv: str,
        h_out: float,
        h_in: float
    ) -> WallStatusValues:
    """Calculation of the status(temperature, heat flow, etc.) of the ventilation layer.

    Args:
        theta_e: the outdoor temperature, degrees
        theta_r: the indoor temperature, degrees
        j_surf: the solar irradiance on the exterior surface, W/m2
        a_surf: the solar absorption ratio on the exterior surface, -
        c_1: the thermal conductance of the outside material, W/m2K
        c_2: the thermal conductance of the inside material, W/m2K
        l_h: the length of the ventilation layer, m
        l_w: the width of the ventilation layer, m
        l_d: the thickness of the ventilation layer, m
        angle: the angle of the ventilation layer, degrees
        v_a: the mean air velocity of the ventilation layer, m/s
        emissivity_1: the emissivity of the surface 1 facing the ventilation layer, -
        emissivity_2: the emissivity of the surface 2 facing the ventilation layer, -
        calc_mode_h_cv: the calculation mode for the convective heat transfer coefficient
        calc_mode_h_rv: the calculation mode for the radiative heat transfer coefficient
        h_out: the overall heat transfer coefficient of the inside, W/m2K
        h_in: the overall heat transfer coefficient of the outside, W/m2K
    Returns:
        the status(temperature, heat flow, etc.) of the ventilation layer
        通気層の状態値（通気層の各層の温度、各層の熱収支、対流熱伝達率、放射熱伝達率、最適化の終了ステータス、終了メッセージ）
    """

    # the initial temperature of the points in the ventilation layer
    t0 = theta_e
    t1 = theta_e + (theta_r - theta_e) * 1 / 4
    t2 = theta_e + (theta_r - theta_e) * 2 / 4
    t3 = theta_e + (theta_r - theta_e) * 3 / 4
    t4 = (t1 + t2) / 2
    
    x0 = np.array([t0, t1, t2, t3, t4])

    def f(matrix_temp):
        return _get_heat_balance(
            matrix_temp=matrix_temp, theta_e=theta_e, theta_r=theta_r, j_surf=j_surf, a_surf=a_surf,
            c_1=c_1, c_2=c_2, l_h=l_h, l_w=l_w, l_d=l_d,
            angle=angle, v_a=v_a, eps1=eps_1, eps2=eps_2,
            calc_mode_h_cv=calc_mode_h_cv, calc_mode_h_rv=calc_mode_h_rv, h_out=h_out, h_in=h_in)


    # 通気層内の各層の熱収支式の最適解を収束計算で求める
    optimize_result = optimize.root(fun=f, x0=x0, method='lm')

    if optimize_result.success:

        # the temperatures, degrees
        matrix_temp_fixed = optimize_result.x

        # the heat balance, W/m2
        heat_balance = f(matrix_temp=matrix_temp_fixed)

        # the convective heat transfer coefficient, W/m2K
        h_cv = htc.get_h_cv(calc_mode=calc_mode_h_cv, v_a=v_a, theta_1=matrix_temp_fixed[1], theta_2=matrix_temp_fixed[2], angle=angle, l_h=l_h, l_d=l_d)

        # the effective emissivity
        eps_eff = htc.get_e(eps1=eps_1, eps2=eps_2)

        # the radiative heat transfer coefficient, W/m2K
        h_rv = htc.get_h_rv(eps_eff=eps_eff, calc_mode=calc_mode_h_rv, theta_1=matrix_temp_fixed[1], theta_2=matrix_temp_fixed[2])

        # the heat flow from the outside(inclued solar irradiance) to the exterior surface, W/m2
        q_flow_out = _get_q_flow_out(theta_e=theta_e, a_surf=a_surf, j_surf=j_surf, theta_0=matrix_temp_fixed[0], h_out=h_out)

        # the heat flow from the exterior surface to the surface of exterior side facing the ventilation layer, W/m2
        q_flow_em = _get_q_flow_em(theta_0=matrix_temp_fixed[0], theta_1=matrix_temp_fixed[1], c_1=c_1)

        # the exhausted heat flow from the ventilation layer, W/m2
        q_flow_exhaust = _get_q_flow_exhaust(v_a=v_a, l_d=l_d, l_w=l_w, l_h=l_h, theta_1=matrix_temp_fixed[1], theta_2=matrix_temp_fixed[2], theta_4=matrix_temp_fixed[4], theta_as_in=theta_e, h_cv=h_cv)

        # the heat flow to the inside, W/m2
        q_flow_in = _get_q_flow_in(theta_3=matrix_temp_fixed[3], theta_r=theta_r, h_in=h_in)

        return WallStatusValues(
            matrix_temp=matrix_temp_fixed,
            matrix_heat_balance=heat_balance,
            h_cv=h_cv,
            h_rv=h_rv,
            q_flow_out=q_flow_out,
            q_flow_em=q_flow_em,
            q_flow_exhaust=q_flow_exhaust,
            q_flow_in=q_flow_in,
            is_optimize_succeed=optimize_result.success,
            optimize_status=optimize_result.status,
            optimize_message=optimize_result.message
        )

    # If the optimized result is false, the return values are set to be np.nan.
    else:

        return WallStatusValues(
            matrix_temp=np.full(5, np.nan),
            matrix_heat_balance=np.full(5, np.nan),
            h_cv=np.nan,
            h_rv=np.nan,
            q_flow_out=np.nan,
            q_flow_em=np.nan,
            q_flow_exhaust=np.nan,
            q_flow_in=np.nan,
            is_optimize_succeed=optimize_result.success,
            optimize_status=optimize_result.status,
            optimize_message=optimize_result.message
        )


def _get_q_flow_out(theta_e: float, a_surf: float, j_surf: float, theta_0: float, h_out: float) -> float:
    """Calculate the heat flow from the outside to the exterior surface

    Args:
        theta_e: the outdoor temperature, degrees
        a_surf: the solar absorption ratio on the exterior surface, -
        j_surf: the solar irradiance on the exterior surface, W/m2
        theta_0: the temperature on the exterior surface, degrees
        h_out: the overall heat transfer coefficient, W/m2K
    Returns
        the heat flow from the outside to the exterior surface, W/m2
    """

    # SAT temperature, degrees
    theta_sat = theta_e + (a_surf * j_surf) / h_out

    return h_out * (theta_sat - theta_0)


def _get_q_flow_em(theta_0: float, theta_1: float, c_1: float) -> float:
    """Calculate the heat flow from the exterior surface to the surface of exterior side facing the ventilation layer

    Args:
        theta_0: the temperature on the exterior surface, degrees
        theta_1: the temperature on the surface of the exterior side facing the ventilation layer, degrees
        c_1: the thermal conductance of the outside material, W/m2K
    Returns:
        the heat flow from the exterior surface to the surface of exterior side facing the ventilation layer, W/m2
    """

    return c_1 * (theta_0 - theta_1)


def _get_q_flow_exhaust(v_a: float, l_d: float, l_w: float, l_h: float, theta_1: float, theta_2: float, theta_4: float, theta_as_in: float, h_cv: float) -> float:
    """Calculate the exhausted heat flow from the ventilation layer.

    Args:
        v_a: the mean air velocity of the ventilation layer, m/s
        l_d: the thickness of the ventilation layer, m
        l_w: the width of the ventilation layer, m
        l_h: the length of the ventilation layer, m
        theta_1: the temperature on the surface of the exterior side facing the ventilation layer, degrees
        theta_2: the temperature on the surface of the interior side facing the ventilation layer, degrees
        theta_4: the air temperature of the ventilation layer, degrees
        theta_as_in: the inlet temperature of the ventilation layer, degrees
            this temperature is equal to the outdoor air temperature, degrees
        h_cv: the convective heat transfer coefficient, W/m2K
    Returns:
        the exhausted heat flow from the ventilation layer, W/m2
    """

    if v_a > 0.0:

        # the air flow of the ventilation layer, m3/s
        v_vent = v_a * l_d * l_w

        ec = math.exp(- 2.0 * h_cv * l_w * l_h / (gn.get_c_air() * gn.get_rho_air(theta_4) * v_vent))

        # the air temperature at the outlet of the ventilation layer, degrees
        theta_out = (1.0 - ec) * (theta_1 + theta_2) / 2.0 + ec * theta_as_in

        return gn.get_c_air() * gn.get_rho_air(theta_4) * v_vent * (theta_out - theta_as_in) / (l_w * l_h)

    else:

        return 0.0


def _get_q_flow_in(theta_3: float, theta_r: float, h_in: float) -> float:
    """Calculate the heat flow to the inside.

    Args:
        theta_3: the surface temperature of the inside(room), degrees
        theta_r: the room temperature, degrees
        h_in: the overall heat transfer coefficient, W/m2K
    Returns:
        the heat flow to the inside., W/m2
    """

    return h_in * (theta_3 - theta_r)

