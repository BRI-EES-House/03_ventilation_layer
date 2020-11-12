import math
from global_number import get_gr_air, get_pr_air, get_lambda_air


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
    pr = get_pr_air(t=t_ave)
    # グラスホフ数の計算
    gr = get_gr_air(tw=tw, tf=tf, d=d)

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

    return get_lambda_air(t=t_ave) * nu / d


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
    pr = get_pr_air(t=t_ave)
    # グラスホフ数の計算
    gr = get_gr_air(tw=tw, tf=tf, d=d)

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
    :param L: 代表長さ[m]（垂直：高さ、水平：A/L_peri）
    :param direction: 向き　'v':垂直、'u':上向き、'd':下向き
    :return: 熱伝達率, W/(m2 K)
    '''

    # 膜温度
    t_ave = (tw + tf) / 2.0

    # プラントル数の計算
    pr = get_pr_air(t=t_ave)
    # グラスホフ数の計算
    gr = get_gr_air(tw=tw, tf=tf, d=L)
    # レーリー数の計算
    ra = pr * gr

    nu = 0.0
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

    return nu * get_lambda_air(t=t_ave) / L


def calc_jsme_databook(tw: float, tf: float, L: float, b: float) -> float:

    '''
    日本機械学会　伝熱工学資料改訂第5版
    鉛直平行平板間の自然対流熱伝達率の計算
    :param tw: 表面温度[C]
    :param tf: 竜体温度[C]
    :param L: 代表長さ[m]（垂直：高さ）
    :param b: 平行平板の間隔[m]
    :return: 熱伝達率, W/(m2 K)
    '''

    # 膜温度
    t_ave = (tw + tf) / 2.0

    # プラントル数の計算
    pr = get_pr_air(t=t_ave)
    # print('Pr=', pr)
    # グラスホフ数の計算
    gr = get_gr_air(tw=tw, tf=tf, d=b)
    # print('Gr=', gr)
    # 修正レーリー数の計算
    ra = pr * gr * b / L
    # print('Rab=', ra)
    # print('log(ra)=', math.log(ra))

    nu = (1.0 / 24.0) * ra \
         / (1.0 + (1.0 / 18.0) ** 1.5 * (1.0 + (0.492 / pr) ** (9.0 / 16.0)) ** (2.0 / 3.0) * ra ** (9.0 / 8.0)) ** (2.0 / 3.0)
    # print('nu1=', nu)
    # nu = 1.0 / 18.8 * ra ** 0.85 * (1.0 - math.exp(-52.0 / ra)) ** 0.6
    # print('nu2=', nu)

    return nu * get_lambda_air(t=t_ave) / L


if __name__ == '__main__':
    print(calc_jsme_databook(tw=5, tf=18, L=3.0, b=0.03))
