import math
import global_number as gn


def calc_takeda_method(delta_t: float, direction: str) -> float:

    '''
    武田仁他：最新建築環境工学　第４版
    自然対流熱伝達率の計算
    :param delta_t:温度差[C]
    :param direction:向き　'v':垂直、'u':上向き、'd':下向き
    :return:熱伝達率, W/(m2 K)
    '''

    if direction == 'v':        # 垂直
        return 1.98 * delta_t ** 0.25
    elif direction == 'u':      # 水平下向き
        return 2.67 * delta_t ** 0.25
    elif direction == 'd':      # 水平上向き
        return (0.64 + 0.87) / 2.0 * delta_t ** 0.25


def calc_shase_handbook_method(tw: float, tf: float, d: float, direction: str) -> float:

    '''
    空気調和・衛生工学便覧　第14班
    自然対流熱伝達率の計算
    :param tw: 表面温度[C]
    :param tf: 竜体温度[C]
    :param d: 代表長さ[m]
    :param direction: 向き　'v':垂直、'u':上向き、'd':下向き
    :return: 熱伝達率, W/(m2 K)
    '''

    # 膜温度
    t_ave = (tw + tf) / 2.0

    # プラントル数の計算
    pr = calc_prandtl_number(tw=tw, tf=tf)
    # グラスホフ数の計算
    gr = calc_grashof_number(tw=tw, tf=tf, d=d)

    # ヌセルト数の初期化
    nu = 0.0

    if direction == 'v':        # 垂直
        if gr * pr < 1.0e9:
            f_pr = 3.0 / 4.0 * (pr / (2.4 + 4.9 * math.sqrt(pr) + 5 * pr)) ** 0.25
            nu = 4.0 / 3.0 * f_pr * (gr * pr) ** 0.25
        elif gr * pr >= 1.0e9:
            nu = 0.13 * gr ** (1.0 / 3.0)
    elif direction == 'u':      # 水平上向き
        if 3.0e5 < gr * pr < 3.0e10:
            nu = 0.27 * (gr * pr) ** 0.25
    elif direction == 'd':      # 水平下向き
        if 1.0e5 < gr * pr < 2.0e7:
            nu = 0.54 * (gr * pr) ** 0.25
        elif 2.0e7 < gr * pr < 3.0e10:
            nu = 0.14 * (gr * pr) ** (1.0 / 3.0)

    return gn.get_lambda_air() * nu / d


def calc_kimura_method(tw: float, tf: float, d: float, direction: str) -> float:

    '''
    木村建一　建築設備基礎
    自然対流熱伝達率の計算
    :param tw: 表面温度[C]
    :param tf: 竜体温度[C]
    :param d: 代表長さ[m]
    :param direction: 向き　'v':垂直、'u':上向き、'd':下向き
    :return: 熱伝達率, W/(m2 K)
    '''

    # 膜温度
    t_ave = (tw + tf) / 2.0

    # プラントル数の計算
    pr = calc_prandtl_number(tw=tw, tf=tf)
    # グラスホフ数の計算
    gr = calc_grashof_number(tw=tw, tf=tf, d=d)

    if direction == 'v':        # 垂直
        if gr * pr < 1.0e9:
            return 1.42 * (abs(tw - tf) / d) ** 0.25
        else:
            return 1.31 * abs(tw - tf) ** (1.0 / 3.0)
    elif direction == 'u':      # 水平上向き
        if gr * pr < 2.e7:
            return 2.64 * (abs(tw - tf) / d) ** 0.25
        else:
            return 0.966 * abs(tw - tf) ** (1.0 / 3.0)
    elif direction == 'd':      # 水平下向き
        return 1.31 * (abs(tw - tf) / d) ** 0.25


def calc_ashrae_handbook_method(tw: float, tf: float, L: float, direction: str) -> float:

    '''
    ASHRAE Handbook Fundamentals 2013
    自然対流熱伝達率の計算
    :param tw: 表面温度[C]
    :param tf: 竜体温度[C]
    :param L: 代表長さ[m]
    :param direction: 向き　'v':垂直、'u':上向き、'd':下向き
    :return: 熱伝達率, W/(m2 K)
    '''

    # 膜温度
    t_ave = (tw + tf) / 2.0

    # プラントル数の計算
    pr = calc_prandtl_number(tw=tw, tf=tf)
    # グラスホフ数の計算
    gr = calc_grashof_number(tw=tw, tf=tf, d=d)
    # レーリー数の計算
    ra = pr * gr

    if direction == 'v':        # 垂直
        if 1.0e-1 < ra < 1.0e9:
            nu = 0.68 + 0.67 * ra ** 0.25 / (1.0 + (0.492 / pr) ** (9.0 / 16.0)) ** (4.0 / 9.0)
        elif 1.0e9 < ra < 1.0e12:
            nu = (0.825 + 0.387 * ra ** (1.0 / 6.0) / (1.0 + (0.492 / pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)) ** 2.0
    elif direction == 'u':      # 水平上向き
        if 1.0 < ra < 200.0:
            nu = 0.96 * ra ** (1.0 / 6.0)
        elif 200.0 < ra < 1.0e4:
            nu = 0.59 * ra ** 0.25
        elif 2.2e4 < ra < 8.0e6:
            nu = 0.54 * ra ** 0.25
        elif 8.0e6 < ra < 1.5e9:
            nu = 0.15 * ra ** (1.0 / 3.0)
    elif direction == 'd':      # 水平下向き
        if 1.0e5 < ra < 1.0e10:
            nu = 0.27 * ra ** 0.25

    return nu * gn.get_lambda_air() / L


def calc_grashof_number(tw: float, tf: float, d: float) -> float:

    '''
    グラスホフ数の計算
    :param tw: 表面温度[C]
    :param tf: 竜体温度[C]
    :param d: 代表長さ[m]
    :return: グラスホフ数
    '''

    return gn.get_sgm() * gn.get_beta_air() * abs(tw - tf) * d ** 3.0 \
           / (gn.get_mu_air() / gn.get_rho_air()) ** 2.0


def calc_prandtl_number(tw: float, tf: float) -> float:

    '''
    プラントル数の計算
    :param tw: 表面温度[C]
    :param tf: 竜体温度[C]
    :return: プラントル数
    '''

    # 動粘性係数
    new_air = gn.get_mu_air() / gn.get_rho_air()
    # 温度拡散率
    a_air = gn.get_lambda_air() / (gn.get_rho_air() * gn.get_c_air())

    return new_air / a_air
