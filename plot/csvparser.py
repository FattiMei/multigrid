import pandas as pd


def parse(csv_path):
    with open(csv_path, 'r') as file:
        header = file.readline().strip()

    df = pd.read_csv(csv_path, comment='#')
    columns = header.split(',')

    return (columns, df)
