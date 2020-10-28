import global_number
import ventilation_wall as vw
import envelope_performance_factors as ep


def validation0():
    '''
    通気層を有する壁体の熱貫流率の検証
    通気風量0での検証
    :return:
    '''

    # 固定値の設定
    c_a = global_number.get_c_air()
    rho_a = global_number.get_rho_air()
    h_out = global_number.get_h_out()
    h_in = global_number.get_h_in()

    # パラメータの設定
    theta_e = 0.0
    theta_r = 20.0
    J_surf = 500.0
    a_surf = 0.8
    C_1 = 5.0
    C_2 = 0.5

    parms = vw.Parameters(theta_e=theta_e,
                          theta_r=theta_r,
                          J_surf=J_surf,
                          a_surf=a_surf,
                          C_1=C_1,
                          C_2=C_2,
                          l_h=3.5,
                          l_w=0.455,
                          l_d=0.05,
                          angle=90.0,
                          v_a=0.0,
                          l_s=100.,
                          emissivity_1=0.9,
                          emissivity_2=0.9)

    # 通気層の状態値を取得
    status = vw.get_wall_status_values(parms, c_a, rho_a, h_out, h_in)

    # 温度状態値の出力
    print('各部温度')
    print(status.matrix_temp)
    print('h_cv= ', status.h_cv)
    print('h_rv= ', status.h_rv)

    # 熱流の出力
    print('屋外表面熱流')
    print(vw.get_heat_flow_0(matrix_temp=status.matrix_temp, param=parms, h_out=h_out))

    # 外装材伝導熱流の出力
    print('外装材伝導熱流')
    print(vw.get_heat_flow_1(matrix_temp=status.matrix_temp, param=parms))

    # 通気層熱伝達
    print('通気層熱伝達')
    print(vw.get_heat_flow_2(matrix_temp=status.matrix_temp, param=parms, h_cv=status.h_cv, h_rv=status.h_rv))

    # 通気層排気熱量
    print('通気層からの排気熱量')
    print(vw.get_heat_flow_exhaust(matrix_temp=status.matrix_temp, param=parms, theta_as_in=parms.theta_e, h_cv=status.h_cv, c_a=c_a, rho_a=rho_a))

    # 通気層内表面から通気層空気への熱流
    print('通気層内表面から通気層空気への熱流')
    print(vw.get_heat_flow_convect_vent_layer(matrix_temp=status.matrix_temp, param=parms, h_cv=status.h_cv))

    # 断熱材＋内装材の伝導熱量
    print('断熱材＋内装材の伝導熱量')
    print(vw.get_heat_flow_3(matrix_temp=status.matrix_temp, param=parms))

    # 室内表面熱流の計算
    print('室内表面熱流')
    q = vw.get_heat_flow_4(matrix_temp=status.matrix_temp, param=parms, h_in=h_in)
    print(q)

    # 通気層を考慮した熱貫流率
    h_o = 25.0
    theta_sat = theta_e + a_surf * J_surf / h_o
    Ue = q / (theta_sat - theta_r)
    print('Ue= ', Ue)

    # 通気層から室内までの熱貫流率
    h_i = 9.0
    h_as_s = 9.1                    # 通気層を有する壁体の省エネ基準における屋外側熱伝達率
    Us = 1.0 / (1.0 / h_as_s + 1.0 / C_2 + 1.0 / h_i)
                                    # 通気層から室内までの熱貫流率（屋外側熱伝達率は省エネ基準の値）
    Usd = 1.0 / (1.0 / (status.h_cv + status.h_rv) + 1.0 / C_2 + 1.0 / h_i)
                                    # 通気層から室内までの熱貫流率（屋外側熱伝達率はここで計算した値）
    print('Us= ', Us, 'Usd= ', Usd)

    # Usdを用いて計算した室内側熱流
    theta_as_e = (status.matrix_temp[4][0] * status.h_cv + status.matrix_temp[1][0] * status.h_rv) \
                / (status.h_cv + status.h_rv)
    q_usd = Usd * (theta_as_e - theta_r)
    print('q_usd= ', q_usd)
    Ue_calc = Usd * (theta_as_e - theta_r) / (theta_sat - theta_r)
    print('UsからUeを計算', Ue_calc)
    eta_e = Ue_calc * a_surf / h_o
    print('eta_e= ', eta_e)
    print('Usから熱流を計算', Ue * (theta_e - theta_r) + eta_e * J_surf)

    # 層構成から計算した熱貫流率
    Ue_from_layer = 1.0 / (1.0 / h_o + 1.0 / h_i + 1.0 / C_1 + 1.0 / C_2 + 1.0 / (status.h_cv + status.h_rv))
    print('Ue_from_layer= ', Ue_from_layer)


if __name__ == '__main__':

    validation0()
