import pandas as pd
import glob
import re

rgx = re.compile('(\S+)\.fve\.abundance\.tsv')
dfs = []
for f in glob.glob('*.fve.abundance.tsv'):
    full_df = pd.read_csv(f, sep='\t')
    m = rgx.search(f)
    if m:
        sample_name = m.group(1)
        sorted_df = pd.read_csv(f'{sample_name}.fve.sorted.abundance.tsv',
            sep='\t', skiprows=1, header=None, usecols=[0,3,4,5,6],
            names=['target_id', 'estimated_abundance', 'support', 'predicted_support', 'c_ratio'])
        final_df = full_df.merge(sorted_df, how='left', on='target_id')
        final_df['sample'] = sample_name
        dfs.append(final_df)

final_df = pd.concat(dfs)
final_df = final_df[['sample', 'target_id', 'length', 'eff_length', 'est_counts', 'tpm',
       'estimated_abundance', 'support', 'predicted_support', 'c_ratio']]
final_df.to_csv('gov.10m.vs.pelagiphages.fve.tsv.gz', sep='\t', index=False)
