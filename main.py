import pandas as pd

import ventilation_wall_parameters as vwp


if __name__ == '__main__':

    df = pd.DataFrame(vwp.get_wall_status_data_by_detailed_calculation("detailed", "detailed"))
    df.to_csv("wall_status_data_frame_detailed.csv")
