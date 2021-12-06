import matplotlib.pyplot as plt
import sys
import pandas as pd
import numpy as np
from scipy import stats

file_name = sys.argv[1]
sheet_info = sys.argv[2]
sheet_name_case = 'case_' + sheet_info
sheet_name_control = 'control_' + sheet_info

to_analyse_case = pd.read_excel(file_name, sheet_name=sheet_name_case)
to_analyse_control = pd.read_excel(file_name, sheet_name=sheet_name_control)

res_mann_outgroup = pd.DataFrame()
res_noise = pd.DataFrame()
res_t_test = pd.DataFrame()

conc = pd.DataFrame()

# select not-0 elements 
to_analyse_case = to_analyse_case[to_analyse_case['count'] > 0]
to_analyse_case = to_analyse_control[to_analyse_control['count'] > 0]

# calculate there's how many dup in control, how many del
length = len(to_analyse_control.index)
dup_length = len(to_analyse_control.loc[to_analyse_control['cnv_type']=='Dup'].index)
deletion_length = len(to_analyse_control.loc[to_analyse_control['cnv_type']=='Del'].index)

# drop what's irrelevant 
to_analyse_case = to_analyse_case.drop(to_analyse_case.columns[0:8], axis=1, inplace=True)
to_analyse_control = to_analyse_control.drop(to_analyse_control.columns[0:8], axis=1, inplace=True)

# calculate control elements 
tissues_control = to_analyse_control.columns.values
 for tissue in tissues:
    # concat all out group values to a separate numpy array
    ingroup = to_analyse_control.loc[:, to_analyse_control.columns == tissue].values
    ingroup_list = np.concatenate(ingroup, axis=None)
    outgroup = to_analyse_control.loc[:, to_analyse_control.columns != tissue].values
    outgroup_list = np.concatenate(outgroup, axis=None)
    #print(len(outgroup_list))
    # perform t-test on: 1. one-sample t-test between background noise mean + ingroup;
    dev_from_back_noise = stats.ttest_1samp(ingroup_list, mean_background)
    conc.at['one_sample_noise' tissue] = dev_from_back_noise.pvalue
    # print(len(dev_from_back_noise.pvalue))
    # 2. independent t-test between outgroup values and ingroup values.
    dev_from_outgroup = stats.ttest_ind(ingroup_list, outgroup_list, equal_var=False)
    conc.at['ingroup_outgroup_ind_t', tissue] = dev_from_outgroup.pvalue
    # 3. possible independent t-test between the mean of ingroup and outgroup
    mann = stats.mannwhitneyu(ingroup_list, outgroup_list)
    conc.at['mann', tissue] = mann.pvalue

#calculate case elements
for times in range(0,10):
    dup = to_analyse_case.loc[to_analyse_case['cnv_type']=='Dup']
    deletion = to_analyse_case.loc[to_analyse_case['cnv_type']=='Del']
    dup_resampled = dup.sample(n=dup_length)
    deletion_resampled = deletion.sample(n=deletion_length)
    to_analyse = pd.concat([dup_resampled,deletion_resampled])
    to_analyse = to_analyse.drop(columns=['cnv_type'])
    # parse through the columns
    tissues = to_analyse.columns.values
    # concat all values to a separate numpy array
    background_noise = np.concatenate(to_analyse.values, axis=None)
    mean_background = np.mean(background_noise)

    for tissue in tissues:
        # concat all out group values to a separate numpy array
        ingroup = to_analyse.loc[:, to_analyse.columns == tissue].values
        ingroup_list = np.concatenate(ingroup, axis=None)
        outgroup = to_analyse.loc[:, to_analyse.columns != tissue].values
        outgroup_list = np.concatenate(outgroup, axis=None)
        #print(len(outgroup_list))
        # perform t-test on: 1. one-sample t-test between background noise mean + ingroup;
        dev_from_back_noise = stats.ttest_1samp(ingroup_list, mean_background)
        res_noise.at[times, tissue] = dev_from_back_noise.pvalue
        # print(len(dev_from_back_noise.pvalue))
        # 2. independent t-test between outgroup values and ingroup values.
        dev_from_outgroup = stats.ttest_ind(ingroup_list, outgroup_list, equal_var=False)
        res_t_test.at[times, tissue] = dev_from_outgroup.pvalue
        # 3. possible independent t-test between the mean of ingroup and outgroup
        mann = stats.mannwhitneyu(ingroup_list, outgroup_list)
        res_mann_outgroup.at[times, tissue] = mann.pvalue
    #    mean_t = np.mean(ingroup)
    #    mean_o = np.mean(outgroup, axis=1)
    #    dev_from_outgroup_v = stats.ttest_ind(mean_t, mean_o, equal_var=True)
    #    dev_from_outgroup = stats.ttest_ind(mean_t, mean_o, equal_var=False)

outfile1 = 'pvalues+one_sample+' + file_name + '+' + sheet_name_case
outfile2 = 'pvalues+ttest' + file_name + '+' + sheet_name_case
outfile3 = 'pvalues+mann' + file_name + '+' + sheet_name_case
outfile4 = 'pvalues+' file_name + '+' + sheet_name_control
res_noise.to_excel(outfile1)
res_t_test.to_excel(outfile2)
res_mann_outgroup.to_excel(outfile3)
conc.to_excel(outfile4)




