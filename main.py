import pandas as pd

import ventilation_wall_parameters as vwp


if __name__ == '__main__':

    # 詳細計算
    print("Detailed Calculation")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation("detailed", "detailed"))
    df.to_csv("wall_status_data_frame_detailed.csv")

    # 放射熱伝達率の検証： 冬期条件の簡易計算
    print("Simplified Calculation: h_rv_winter")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_winter"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_winter.csv")

    # 放射熱伝達率の検証： 夏期条件の簡易計算
    print("Simplified Calculation: h_rv_summer")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_summer"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_summer.csv")

    # 放射熱伝達率の検証： 放射熱伝達率ゼロ
    print("Simplified Calculation: h_rv_zero")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_zero"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_zero.csv")

    # 放射熱伝達率の検証：　通年の簡易計算
    print("Simplified Calculation: h_rv_all_season")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="detailed", calc_mode_h_rv="simplified_all_season"))
    df.to_csv("wall_status_data_frame_h_rv_simplified_all_season.csv")

    # 対流熱伝達率の検証： 冬期条件の簡易計算
    print("Simplified Calculation: h_cv_winter")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="simplified_winter", calc_mode_h_rv="detailed"))
    df.to_csv("wall_status_data_frame_h_cv_simplified_winter.csv")

    # 対流熱伝達率の検証： 夏期条件の簡易計算
    print("Simplified Calculation: h_cv_summer")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="simplified_summer", calc_mode_h_rv="detailed"))
    df.to_csv("wall_status_data_frame_h_cv_simplified_summer.csv")

    # 対流熱伝達率の検証：　通年の簡易計算
    print("Simplified Calculation: h_cv_all_season")
    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation(calc_mode_h_cv="simplified_all_season", calc_mode_h_rv="detailed"))
    df.to_csv("wall_status_data_frame_h_cv_simplified_all_season.csv")
