import pandas as pd

def loaddata(datafile,header=0):
    dataset = pd.read_excel(datafile, header=header)
    print('Dataset loaded.')
    return dataset