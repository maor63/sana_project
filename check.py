import pandas as pd
from scipy import stats

file_name = 'y0_y2_ec_ics_s3_fbeta_0.1_0.4_5min_iter.csv'
df = pd.read_excel('%s.xlsx' % file_name)
gb = df.groupby('method')['nc']
groups = gb.groups.keys()
df_results = pd.DataFrame(index=list(groups)[1:], columns=list(groups)[:-1])
for i, group1 in enumerate(list(groups)[:-1]):
    group1_nc = gb.get_group(group1)
    col = [] + [None] * i
    for j, group2 in enumerate(list(groups)[i + 1:]):
        group2_nc = gb.get_group(group2)

        statistic, pvalue = stats.ttest_ind(group1_nc, group2_nc)
        col.append(pvalue)
    df_results[group1] = col
df_results.T.to_csv('ttest_{}.csv'.format(file_name))
