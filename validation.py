import ventilation_wall as vw
import envelope_performance_factors as ep


def validation():
    '''
    通気層を有する壁体の熱貫流率の検証
    :return:
    '''

    # パラメータの設定
    parms = vw.Parameters(theta_e=0.0,
                          theta_r=20.0,
                          J_surf=500.0,
                          a_surf=0.8,
                          C_1=5.0,
                          C_2=0.5,
                          l_h=3.5,
                          l_w=0.455,
                          l_d=0.05,
                          angle=90.0,
                          v_a=0.5,
                          l_s=100.,
                          emissivity_1=0.9,
                          emissivity_2=0.9)

    # 通気層の状態値を取得
    status = vw.get_wall_status_values(parms)

    # 温度状態値の出力
    print('各部温度')
    print(status.matrix_temp)

    # 熱流の出力
    print('屋外表面熱流')
    print(vw.get_heat_flow_0(matrix_temp=status.matrix_temp, param=parms))

    # 外装材伝導熱流の出力
    print('外装材伝導熱流')
    print(vw.get_heat_flow_1(matrix_temp=status.matrix_temp, param=parms))

    # 通気層熱伝達
    print('通気層熱伝達')
    print(vw.get_heat_flow_2(matrix_temp=status.matrix_temp, param=parms, h_cv=status.h_cv, h_rv=status.h_rv))

    # 通気層排気熱量
    print('通気層からの排気熱量')
    print(vw.get_heat_flow_exhaust(matrix_temp=status.matrix_temp, param=parms, theta_as_in=parms.theta_e, h_cv=status.h_cv))

    # 通気層内表面から通気層空気への熱流
    print('通気層内表面から通気層空気への熱流')
    print(vw.get_heat_flow_convect_vent_layer(matrix_temp=status.matrix_temp, param=parms, h_cv=status.h_cv))

    # 断熱材＋内装材の伝導熱量
    print('断熱材＋内装材の伝導熱量')
    print(vw.get_heat_flow_3(matrix_temp=status.matrix_temp, param=parms))

    # 室内表面熱流の計算
    print('室内表面熱流')
    print(vw.get_heat_flow_4(matrix_temp=status.matrix_temp, param=parms))


if __name__ == '__main__':

    validation()
