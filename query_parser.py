import pandas as pd


def csv(filename):
    csv = pd.read_csv(filename)
    csv.index += 1  # renumbers indices to start at 1

    labels = csv['Labels'].to_dict()

    pfams_raw = csv['Pfams'].to_dict()
    pfams_split = {k: v.split(';') for k, v in pfams_raw.items()}

    return labels, pfams_split
