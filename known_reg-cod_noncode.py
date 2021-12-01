import pandas as pd
import numpy as np
import re
import sys
import openpyxl

threshold = sys.argv[1]
option = sys.argv[2] # 1 for count, 2 for coding gene count, 3 for noncoding gene count, 4 for tissue enrichment 
sheet = sys.argv[3]
print("threshold = " + threshold + ", option is " + option)
region_info_loc = 'case_control_specific_regions_filtered.xlsx'
gene_vs_tissue_loc = ''
pppd = ''
region_all = ''

if (option == '1'):
  gene_vs_tissue_loc = 'FANTOM5_genes_vs_traits_expession_avg_25Feb2021.txt'
  pppd = pd.read_csv(gene_vs_tissue_loc, delimiter='\t')
elif option == '2':
  gene_vs_tissue_loc = 'FANTOM5_genes_vs_traits_expession_avg_coding.csv'
  pppd = pd.read_csv(gene_vs_tissue_loc)
elif option == '3':
  gene_vs_tissue_loc = 'FANTOM5_gene_vs_expression_noncoding.csv'
  pppd = pd.read_csv(gene_vs_tissue_loc)
elif option == '4':
  gene_vs_tissue_loc = 'FantomCAT_robust_genes_vs_traits_pvalue_nonzero_04072021_new.txt'
  pppd = pd.read_csv(gene_vs_tissue_loc, delimiter='\t')

#coding: FANTOM5_genes_vs_traits_expession_avg_coding.csv FANTOM5_genes_vs_traits_expession_avg_25Feb2021
# noncoding FANTOM5_gene_vs_expression_noncoding.csv
# tissue enrichment: FantomCAT_robust_genes_vs_traits_pvalue_nonzero_04072021_new


if (sheet == 'case'):
  region_all = pd.read_excel(region_info_loc, sheet_name='case_specific_regions')

elif sheet == 'control':
  region_all = pd.read_excel(region_info_loc, sheet_name='control_specific_regions')

# data cleaning on FANTOM6 data
#t_rearranged = pppd.applymap(lambda x: 1 if x.dtype == 'float64' and x > int(threshold) else 0)
region = region_all.iloc[:,0:9]
chrs = pppd.loc[:,'chr']

occurred_gene = []

def get_chr(row):
    chr_num = re.search('X|Y|[0-9]*$', row).group(0)
    return chr_num

chrs_1 = chrs.apply(lambda x : get_chr(x))
pppd.loc[:, 'chr']= chrs_1

# select gene ids in range
# transfer relatove information
# calculate count
# append in table


def select_relative_gene_id(row):
  #global occurred_gene
  
  chro = str(row['chr'])
  r1 = row['start']
  r2 = row['end']
  #print(chr)
  # firstly selected all genes that overlaps with the region
  if (chro == '23'):
    chro = 'X'
  elif chro == '24':
    chro = 'Y'
  selected_gene_id = pppd.loc[pppd['chr']==chro]
  selected_gene_id = selected_gene_id.loc[((pppd['start'] < r2) & (pppd['end'] > r1))]
  
  
  
  # then cross compare with the existing list, start from empty []. select genes that are uniquely from the existing genes
  #print(len(occurred_gene))
  #if (len(occurred_gene) > 0):
  #  selected_gene_id = selected_gene_id[~selected_gene_id['geneID'].isin(occurred_gene)]
  # then add unique element of the list to the occurring list
  #occurred_gene = np.concatenate([selected_gene_id['geneID'].values, occurred_gene], axis=None)

  # count the existing list
  selected_genename = selected_gene_id['geneID']
  selected_gene_id = selected_gene_id.select_dtypes(include=['float64'])
  if (option == '4'):
    selected_gene_id = selected_gene_id.applymap(lambda x: 1 if x < float(threshold) else 0)
  else:
    selected_gene_id = selected_gene_id.applymap(lambda x: 1 if x > int(threshold) else 0)
  counted = selected_gene_id.sum()
  selected_gene_id['geneID']=selected_genename
  ##print(selected_gene_id.loc[selected_gene_id['nervous system']==1])
  no_of_v = selected_gene_id.count().values[0]
  count_series = pd.Series([no_of_v], index=['count'])
  #print(counted)
  to_fill = pd.concat([count_series, counted])
  #print(to_fill)
  region_all.iloc[row.name, 9:(len(region_all.columns))] = to_fill
  

region.apply(lambda x : select_relative_gene_id(x),axis=1)
if (option == '1'):
  out_name = 'count_patient_reg_' + sheet + '_' + threshold + '.xlsx'
elif option == '2':
  out_name = 'coding_patient_reg_' + sheet + '_' + threshold + '.xlsx'
elif option == '3':
  out_name = 'noncoding_patient_reg_' + sheet + '_' + threshold + '.xlsx'
elif option == '4':
  out_name = 't_en_patient_reg_' + sheet + '_' + threshold + '.xlsx'
region_all.to_excel(out_name)

