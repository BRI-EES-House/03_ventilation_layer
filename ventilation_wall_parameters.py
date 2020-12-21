import itertools
import pandas as pd
import numpy as np
import global_number
import ventilation_wall as vw
import envelope_performance_factors as epf


class Log:
    def write(self, msg):
        """
        np.seterrのエラーログ記録用の処理（エラー発生時はコンソールにメッセージを出力する）
        :param msg:　エラーメッセージ
        :return: なし
        """
        print("LOG: %s" % msg)


def get_parameter_list() -> object:
    """
    複数のパラメータの総当たりの組み合わせ（直積）のリストを作成する
    :param なし
    :return: 総当たりのパラメータリスト
    """

    # Note: 各パラメータには、一部を除いて想定される上下限値と中央値の3点を与える
    theta_e = np.array([-20.0, np.median([-20.0, 40.0]), 40.0], dtype=float)        # 外気温度, degree C
    theta_r = np.array([20.0, np.median([20.0, 27.0]), 27.0], dtype=float)          # 室内温度,　degree C
    j_surf = np.array([0.0, np.median([0.0, 1000.0]), 1000.0], dtype=float)         # 外気側表面に入射する日射量, W/m2
    a_surf = np.array([0.0, np.median([0.0, 1.0]), 1.0], dtype=float)               # 外気側表面日射吸収率
    C_1 = np.array([0.5, np.median([0.5, 100.0]), 100.0], dtype=float)              # 外気側部材の熱コンダクタンス,W/(m2・K)
    C_2 = np.array([0.1, np.median([0.1, 5.0]), 5.0], dtype=float)                  # 室内側部材の熱コンダクタンス, W/(m2・K)
    l_h = np.array([3.0, np.median([3.0, 12.0]), 12.0], dtype=float)                # 通気層の長さ, m
    l_w = np.array([0.05, np.median([0.05, 10.0]), 10.0], dtype=float)              # 通気層の幅, m
    l_d = np.array([0.05, np.median([0.05, 0.3]), 0.3], dtype=float)                # 通気層の厚さ, m
    angle = np.array([0.0, np.median([0.0, 90.0]), 90.0], dtype=float)              # 通気層の傾斜角, °
    v_a = np.array([0.0, np.median([0.0, 1.0]), 1.0], dtype=float)                  # 通気層の平均風速, m/s
    l_s = [0.45]                                                                    # 通気胴縁または垂木の間隔, m
    emissivity_1 = [0.9]                                                            # 通気層に面する面1の放射率, -
    emissivity_2 = np.array([0.1, np.median([0.1, 0.9]), 0.9], dtype=float)         # 通気層に面する面2の放射率, -

    parameter_list = list(itertools.product(theta_e, theta_r, j_surf, a_surf, C_1, C_2, l_h, l_w, l_d, angle,
                                            v_a, l_s, emissivity_1, emissivity_2))

    return parameter_list


def get_wall_status_data_frame() -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、各ケースの計算結果を保有するDataFrameを作成する
    :param なし
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    parameter_name = ['theta_e', 'theta_r', 'j_surf', 'a_surf', 'C_1', 'C_2', 'l_h', 'l_w', 'l_d', 'angle',
                      'v_a', 'l_s', 'emissivity_1', 'emissivity_2']
    df = pd.DataFrame(get_parameter_list(), columns=parameter_name)

    # 固定値の設定
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    # 計算結果格納用配列を用意
    theta_out_surf = []     # 外気側表面温度[℃]
    theta_1_surf = []       # 通気層に面する面1の表面温度[℃]
    theta_2_surf = []       # 通気層に面する面1の表面温度[℃]
    theta_in_surf = []      # 室内側表面温度[℃]
    theta_as_ave = []       # 通気層の平均温度[℃]
    h_cv = []               # 通気層の対流熱伝達率[W/(m2・K)]
    h_rv = []               # 通気層の放射熱伝達率[W/(m2・K)]
    theta_as_e = []         # 通気層の等価温度[℃]
    k_e = []                # 通気層を有する壁体の相当熱貫流率を求めるための補正係数[-]

    # エラーログ出力用の設定
    log = Log()
    saved_handler = np.seterrcall(log)

    with np.errstate(all='log'):  # withスコープ内でエラーが出た場合、Logを出力する
        for row in df.itertuples():
            print(row[0])
            # パラメータを設定
            parms = (vw.Parameters(theta_e=row.theta_e,
                                   theta_r=row.theta_r,
                                   J_surf=row.j_surf,
                                   a_surf=row.a_surf,
                                   C_1=row.C_1,
                                   C_2=row.C_2,
                                   l_h=row.l_h,
                                   l_w=row.l_w,
                                   l_d=row.l_d,
                                   angle=row.angle,
                                   v_a=row.v_a,
                                   l_s=row.l_s,
                                   emissivity_1=row.emissivity_1,
                                   emissivity_2=row.emissivity_2))

            # 通気層の状態値を取得
            status = vw.get_wall_status_values(parms, h_out, h_in)
            theta_out_surf.append(status.matrix_temp[0])
            theta_1_surf.append(status.matrix_temp[1])
            theta_2_surf.append(status.matrix_temp[2])
            theta_in_surf.append(status.matrix_temp[3])
            theta_as_ave.append(status.matrix_temp[4])
            h_cv.append(status.h_cv)
            h_rv.append(status.h_rv)

            # 通気層の等価温度を取得
            theta_as_e.append(epf.get_theata_as_e(status.matrix_temp[4], status.matrix_temp[1],
                                                  status.h_cv, status.h_rv))

            # 相当外気温度を計算
            theta_SAT = epf.get_theta_SAT(row.theta_e, row.a_surf, row.j_surf, h_out)

            # 通気層を有する壁体の相当熱貫流率を求めるための補正係数を取得
            k_e.append(epf.get_weight_factor_of_u_s_dash(status.matrix_temp[4], row.theta_r, theta_SAT))

    # 計算結果をDataFrameに追加
    df['theta_out_surf'] = theta_out_surf
    df['theta_1_surf'] = theta_1_surf
    df['theta_2_surf'] = theta_2_surf
    df['theta_in_surf'] = theta_in_surf
    df['theta_as_ave'] = theta_as_ave
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['theta_as_e'] = theta_as_e
    df['k_e'] = k_e

    return df


def dump_csv_all_case_result():
    # 総当たりのパラメータと計算結果を取得
    df = pd.DataFrame(get_wall_status_data_frame())
    # CSV出力
    df.to_csv("wall_status_data_frame.csv")


if __name__ == '__main__':

    dump_csv_all_case_result()


# デバッグ用
# dump_csv_all_case_result()
# print(np.median([-20,40]))
