import sys
import gzip
import pandas as pd

file = sys.argv[1]
norm_hc = sys.argv[2]

vcf = gzip.open(norm_hc, 'rb')
num = 0

for line in vcf:
    if line.decode().startswith('##'):
        num += 1
    elif line.decode().startswith('#'):
        sample_name = line.decode().split('\t')[9]
    else:
        break
vcf.close()

df = pd.read_feather(file)
dfHC = pd.read_table(norm_hc, compression='gzip', skiprows=num)

_tmp_chr = list(dfHC['#CHROM'])
dfHC['#CHROM'] = [_.lstrip('chr') for _ in _tmp_chr]
_tmp_gt = list(dfHC[sample_name])
for n, _gt in enumerate(_tmp_gt):
    gt = _gt.split(':', 1)[0]
    if gt in ['0/1', '0|1']:
        _tmp_gt[n] = '0/1'
    elif gt in ['1/1', '1|1']:
        _tmp_gt[n] = '1/1'
    elif gt in ['1/0', '1|0']:
        _tmp_gt[n] = '1/0'
    else:
        _tmp_gt[n] = gt
dfHC['GT'] = _tmp_gt

dfTg = pd.DataFrame(dfHC, columns=['#CHROM', 'POS', 'REF', 'ALT', "GT"])
dfTg['tag'] = 1
df2 = pd.merge(df, dfTg, how='left', on=['#CHROM', 'POS', 'REF', 'ALT', "GT"])
df2['tag'].fillna(0, inplace=True)
df2['tag'] = df2['tag'].astype(int)

df2.to_feather(file + '_withTag.feather')
