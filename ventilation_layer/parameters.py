from dataclasses import dataclass
import numpy as np
import pandas as pd
import itertools

@dataclass
class Parameters:

    # 外気温度, deg.C
    theta_e: float

    # 室温, deg.C
    theta_r: float
 
    # 日射量, W/m2
    j_surf: float

    # 日射吸収率, -
    a_surf: float

    # 外装材の熱コンダクタンス, W/m2K
    C_1: float
    
    # 室内側部材の熱コンダクタンス, W/m2K
    C_2: float
    
    # 通気層の縦方向の長さ, m
    l_h: float
    
    # 通気層の横方向の長さ, m
    l_w: float

    # 通気層の幅, m
    l_d: float
    
    # 通気層の傾き, deg.
    angle: float

    # 通気層の風速, m/s
    v_a: float
    
    # 胴縁と垂木の間の幅, m
    l_s: float
    
    # 通気層に面する面の長波長放射率（室外側）, -
    emissivity_1: float

    # 通気層に面する面の長波長放射率（室内側）, -
    emissivity_2: float


@dataclass
class ParametersOptions:

    # 外気温度, deg.C
    theta_e: list
 
    # 日射量, W/m2
    j_surf: list

    # 日射吸収率, -
    a_surf: list

    # 外装材の熱コンダクタンス, W/m2K
    C_1: list
    
    # 室内側部材の熱コンダクタンス, W/m2K
    C_2: list
    
    # 通気層の縦方向の長さ, m
    l_h: list
    
    # 通気層の横方向の長さ, m
    l_w: list

    # 通気層の幅, m
    l_d: list
    
    # 通気層の傾き, deg.
    angle: list

    # 通気層の風速, m/s
    v_a: list
    
    # 胴縁と垂木の間の幅, m
    l_s: list
    
    # 通気層に面する面の長波長放射率（室外側）, -
    emissivity_1: list

    # 通気層に面する面の長波長放射率（室内側）, -
    emissivity_2: list

    @classmethod
    def SetParameters(
        cls,
        theta_e=None,
        j_surf=None,
        a_surf=None,
        C_1=None,
        C_2=None,
        l_h=None,
        l_w=None,
        l_d=None,
        angle=None,
        v_a=None,
        l_s=None,
        emissivity_1=None,
        emissivity_2=None
    ):

        theta_e = [-10.0, 0.0, 10.0, 25.0, 30.0, 35.0] if theta_e is None else theta_e

        j_surf = [0.0, 500.0, 1000.0] if j_surf is None else j_surf

        a_surf = [0.0, 0.5, 1.0] if a_surf is None else a_surf
    
        C_1 = [0.5, 50.25, 100.0] if C_1 is None else C_1
    
        C_2 = [0.1, 2.55, 5.0] if C_2 is None else C_2
    
        l_h = [3.0, 7.5, 12.0] if l_h is None else l_h
    
        l_w = [0.05, 5.025, 10.0] if l_w is None else l_w

        l_d = [0.005, 0.0175, 0.03] if l_d is None else l_d
    
        angle = [0.0, 45.0, 90.0] if angle is None else angle

        v_a = [0.0, 0.5, 1.0] if v_a is None else v_a
    
        l_s = [0.45] if l_s is None else l_s
    
        emissivity_1 = [0.9] if emissivity_1 is None else emissivity_1

        emissivity_2 = [0.1, 0.5, 0.9] if emissivity_2 is None else emissivity_2

        return ParametersOptions(
            theta_e=theta_e,
            j_surf=j_surf,
            a_surf=a_surf,
            C_1=C_1,
            C_2=C_2,
            l_h=l_h,
            l_w=l_w,
            l_d=l_d,
            angle=angle,
            v_a=v_a,
            l_s=l_s,
            emissivity_1=emissivity_1,
            emissivity_2=emissivity_2
        )

    def get_df(self):

        parameter_list = list(
            itertools.product(
                self.theta_e,
                self.j_surf,
                self.a_surf,
                self.C_1,
                self.C_2,
                self.l_h,
                self.l_w,
                self.l_d,
                self.angle,
                self.v_a,
                self.l_s,
                self.emissivity_1,
                self.emissivity_2
            )
        )
    
        parameter_name = ['theta_e', 'j_surf', 'a_surf', 'C_1', 'C_2', 'l_h', 'l_w', 'l_d', 'angle', 'v_a', 'l_s', 'emissivity_1', 'emissivity_2']

        df = pd.DataFrame(parameter_list, columns=parameter_name)

        # Give the temperature of 20.0 degrees for winter and 27.0 degrees for summer as the indoor temperature.
        df['theta_r'] = np.where(df.theta_e > 20.0, 27.0, 20.0)

        return df

    def check_options(self):

        df = self.get_df()

        print('外気温度：' + str(df['theta_e'].unique().tolist()))
        print('外気側表面に入射する日射量：' + str(df['j_surf'].unique().tolist()))
        print('外装材の日射吸収率：' + str(df['a_surf'].unique().tolist()))
        print('外装材の熱コンダクタンス：' + str(df['C_1'].unique().tolist()))
        print('断熱層の熱コンダクタンス：' + str(df['C_2'].unique().tolist()))
        print('通気層の長さ：' + str(df['l_h'].unique().tolist()))
        print('通気層の幅：' + str(df['l_w'].unique().tolist()))
        print('通気層の厚さ：' + str(df['l_d'].unique().tolist()))
        print('通気層の傾斜角：' + str(df['angle'].unique().tolist()))
        print('通気層の平均風速：' + str(df['v_a'].unique().tolist()))
        print('胴縁と垂木との間の幅：' + str(df['l_s'].unique().tolist()))
        print('通気層に面する面１の放射率：' + str(df['emissivity_1'].unique().tolist()))
        print('通気層に面する面２の放射率：' + str(df['emissivity_2'].unique().tolist()))

    def get_parameters_list(self):

        df = self.get_df()

        return [
            Parameters(
                theta_e=row.theta_e,
                theta_r=row.theta_r,
                j_surf=row.j_surf,
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
                emissivity_2=row.emissivity_2
            )
            for row in df.itertuples()
        ]
