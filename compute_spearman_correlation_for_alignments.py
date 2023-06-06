import pandas as pd
from scipy.stats import spearmanr, pearsonr, kendalltau
import os
import numpy as np

base_path = 'output_alignment/alignment_checkpoint/'
dfs = []
for file in os.listdir(base_path):
    dfs.append(pd.read_csv(os.path.join(base_path, file)))
df1 = pd.concat(dfs)

# df1 = pd.read_csv('output_alignment/alignmnets_y0_com_bal_v1.csv')

df1_beta_star = (df1['E1'] / df1['Ea_hat']) ** 0.5
df1['FBeta*'] = ((1 + df1_beta_star ** 2) * df1['Ea']) / (df1['E1'] + df1_beta_star ** 2 * df1['Ea_hat'])
measures = ['NC', 'EC', 'ICS', 'S3', 'F1', 'ACC', 'BM', 'MK', 'MCC', 'F0.1', 'F0.33', 'F3', 'F10', 'FBeta*', 'BMK2', 'Ilia']
df1['BMK2'] = (df1['BM'] + df1['MK']) / 2.0
df1['Ilia'] = (df1['TP'] / df1['E1'] + df1['TN'] / (df1['omega'] - df1['E1'])) / 2.0
df2 = pd.DataFrame(columns=measures)
static_measures = ['NC', 'S3', 'F1', 'ACC', 'MCC']
for measure in static_measures:
    df2[measure] = df1[measure]

dynamic_measures = [('EC', 'ICS'), ('BM', 'MK'), ('F0.1', 'F10'), ('F0.33', 'F3')]
for m1, m2 in dynamic_measures:
    df2[m1] = df1[m2]
    df2[m2] = df2[m1]

df2_beta_star = 1 / ((df1['E1'] / df1['Ea_hat']) ** 0.5)
df2['FBeta*'] = ((1 + df2_beta_star ** 2) * df1['Ea']) / (df1['E1'] + df2_beta_star ** 2 * df1['Ea_hat'])
df2['BMK2'] = (df2['BM'] + df2['MK']) / 2.0
df = pd.concat([df1[measures], df2[measures]])
print('number of rows: ' + str(len(df)))

rho, p_value = spearmanr(df[measures])
pd.DataFrame(rho, columns=measures, index=measures).to_csv('output_alignment/spearman_correlation_ba_er_ws.csv')

# rho = np.corrcoef(df[measures])

pc = pd.DataFrame(columns=measures, index=['NC'])
for measure in measures:
    pc[measure]['NC'] = pearsonr(df['NC'], df[measure])[0]
pc.to_csv('output_alignment/pearson_correlation_ba_er_ws.csv')

kt = pd.DataFrame(columns=measures, index=['NC'])
for measure in measures:
    kt[measure]['NC'] = kendalltau(df['NC'], df[measure])[0]
kt.to_csv('output_alignment/kendalltau_correlation_ba_er_ws.csv')
