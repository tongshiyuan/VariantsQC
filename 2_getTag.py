import sys
import gzip
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', help='input file.')
parser.add_argument('-c', '--hc', help='high confident vcf of giab with bcftool norm.')
parser.add_argument('--xc', action="store_true", help='Whether to include X chromosome.')
args = parser.parse_args()
file = args.infile
norm_hc = args.hc


# identify hc vcf header number
def get_num(vcf_file):
    vcf = gzip.open(vcf_file, 'rb')
    num, sample_name = 0, ''
    for line in vcf:
        if line.decode().startswith('##'):
            num += 1
        elif line.decode().startswith('#'):
            sample_name = line.decode().strip().split('\t')[9]
        else:
            break
    vcf.close()
    return num, sample_name


num, sample_name = get_num(norm_hc)
dfHC = pd.read_table(norm_hc, compression='gzip', skiprows=num)

cl = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
      '21', '22']
if args.xc:
    cl.append('X')
# identify input file
if file.endswith('vcf.gz'):
    num2, sample_name2 = get_num(file)
    df = pd.read_table(file, compression='gzip', skiprows=num2)
    _tmp_chr = df['#CHROM']
    df['#CHROM'] = [_.lstrip('chr') for _ in _tmp_chr]
    sample_info_list = list(df[sample_name2])
    gt_list = []
    for i in sample_info_list:
        gt = i.split(':', 1)[0]
        if gt in ['0/1', '0|1']:
            gt_list.append('0/1')
        elif gt in ['1/1', '1|1']:
            gt_list.append('1/1')
        elif gt in ['1/0', '1|0']:
            gt_list.append('1/0')
        else:
            gt_list.append(gt)
    df['GT'] = gt_list
else:
    df = pd.read_feather(file)
df = df[df['#CHROM'].isin(cl)]
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
print(df2['tag'].value_counts())
df2.to_feather(file + '_withTag.feather')
