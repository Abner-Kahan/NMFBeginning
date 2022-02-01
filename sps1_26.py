import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
data = pd.read_csv('nondaily_study_deidentified_7_16_2021.csv')
ind = 0
ListInd=[]
for n in data['SMOKING_FREQUENCY']:
    if (n == 1 or n == 2):
        ListInd.append(ind)
    ind+=1

print(len (ListInd))
dataFilt = data.loc[ListInd]
print(dataFilt)