import itertools
import pandas as pd
import global_number
import ventilation_wall as vw


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
    :arg:
    :return: 総当たりのパラメータリスト
    """

    # TODO: テスト用にパラメータ数を制限しているため、最終的には想定したパラメータに合わせて変更すること
    theta_e = [round(x * 1.0, 1) for x in range(-20, 50, 10)]    # 外気温度, degree C
    theta_r = [25.0]                                            # 室内温度,　degree C
    j_surf = [round(x * 1.0, 1) for x in range(0, 1100, 200)]    # 外気側表面に入射する日射量, W/m2
    a_surf = [round(x * 0.1, 1) for x in range(0, 11, 2)]       # 外気側表面日射吸収率
    C_1 = [round(x * 0.1, 1) for x in range(5, 20, 5)]          # 外気側部材の熱コンダクタンス,W/(m2・K)
    C_2 = [round(x * 0.1, 1) for x in range(5, 20, 5)]          # 室内側部材の熱コンダクタンス, W/(m2・K)
    l_h = [3.0]                 # 通気層の長さ, m
    l_w = [3.0]                 # 通気層の幅, m
    l_d = [0.3]                 # 通気層の厚さ, m
    angle = [0.0, 45.0, 90.0]   # 通気層の傾斜角, °
    v_a = [1.0]                 # 通気層の平均風速, m/s
    l_s = [0.45]                # 通気胴縁または垂木の間隔, m
    emissivity_1 = [1.0]        # 通気層に面する面1の放射率, -
    emissivity_2 = [1.0]        # 通気層に面する面2の放射率, -

    parameter_list = list(itertools.product(theta_e, theta_r, j_surf, a_surf, C_1, C_2, l_h, l_w, l_d, angle,
                                            v_a, l_s, emissivity_1, emissivity_2))

    return parameter_list


def get_wall_status_data_frame() -> pd.DataFrame:
    """
    通気層を有する壁体の総当たりパラメータを取得し、各ケースの計算結果を保有するDataFrameを作成する
    :arg:
    :return: DataFrame
    """

    # パラメータの総当たりリストを作成する
    parameter_name = ['theta_e', 'theta_r', 'j_surf', 'a_surf', 'C_1', 'C_2', 'l_h', 'l_w', 'l_d', 'angle',
                      'v_a', 'l_s', 'emissivity_1', 'emissivity_2']
    df = pd.DataFrame(get_parameter_list(), columns=parameter_name)

    # 固定値の設定
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    # 計算結果格納用配列
    # TODO: for文を使わずに処理ができないか検討する(for文だと処理が遅い）
    # parms = []
    matrix_temp = []
    h_cv = []
    h_rv = []
    heat_flow_0 = []            # 屋外側表面熱流, W/m2
    heat_flow_exhaust = []      # 通気層からの排気熱量, W/m2
    heat_flow_4 = []            # 室内表面熱流, W/m2

    for row in df.itertuples():
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
        matrix_temp.append(status.matrix_temp)
        h_cv.append(status.h_cv)
        h_rv.append(status.h_rv)

        # 屋外側表面熱流、通気層からの排気熱量、室内表面熱流を取得
        heat_flow_0.append(vw.get_heat_flow_0(matrix_temp=status.matrix_temp, param=parms, h_out=h_out))
        heat_flow_exhaust.append(vw.get_heat_flow_exhaust(matrix_temp=status.matrix_temp, param=parms,
                                                          theta_as_in=parms.theta_e, h_cv=status.h_cv))
        heat_flow_4.append(vw.get_heat_flow_4(matrix_temp=status.matrix_temp, param=parms, h_in=h_in))

    # 計算結果をDataFrameに追加
    # TODO: matrix_tempは配列のままではなく個別の値としてDataFrameに入れる
    df['matrix_temp'] = matrix_temp
    df['h_cv'] = h_cv
    df['h_rv'] = h_rv
    df['heat_flow_0'] = heat_flow_0
    df['heat_flow_exhaust'] = heat_flow_exhaust
    df['heat_flow_4'] = heat_flow_4
    # df['parms'] = parms
    # print(df['parms'])
    # df['matrix_temp'] = vw.get_wall_status_values(df['parms'], h_out, h_in).matrix_temp

    return df


#デバッグ用
#get_wall_status_data_frame()
