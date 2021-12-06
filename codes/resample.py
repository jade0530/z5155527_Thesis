import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
from scipy import stats

file_name = sys.argv[1]
sheet = sys.argv[2]

deletion_length = 44
dup_length = 25

to_analyse_case = pd.read_excel(file_name, sheet_name=sheet)

#calculate case elements
for times in range(0,5):
    dup = to_analyse_case.loc[to_analyse_case['cnv type']=='Dup']
    deletion = to_analyse_case.loc[to_analyse_case['cnv type']=='Del']
    dup_resampled = dup.sample(n=dup_length)
    deletion_resampled = deletion.sample(n=deletion_length)
    resampled = pd.concat([dup_resampled,deletion_resampled])
    outfile = 'resample+' + 'number+' + str(times) + '+' + sheet + '+' + file_name 
    resampled.to_excel(outfile)
    


