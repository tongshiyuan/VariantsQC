import os
import sys
import gzip
import datetime
import argparse
import pandas as pd


def get_args():
    args_dict = {}
    description = 'from vcf to annotated matrix.'
    parser = argparse.ArgumentParser()
    parser.description = description
    parser.add_argument('-i', '--vcf', help='input file of vcf')
    parser.add_argument('-o', '--annofile', help='output file for next annotation')
    parser.add_argument('-r', '--reference', help='reference with fasta format.')
    parser.add_argument('-d', '--humandb',
                        help='directory of annovar database. Must have the reference databases:'
                             'hgxx_refGene, hgxx_rmsk, hgxx_cpgIslandExt, hgxx_genomicSuperDups')
    parser.add_argument('-t', '--tmp', default='./tmp', help='temp directory, [./tmp/]')
    parser.add_argument('--dp', default=0, type=int, help='depth for filter variants. [0]')
    parser.add_argument('--refTag', default='hg38', choices=['hg19', 'hg38'],
                        help='reference of variants hg19/hg38. [hg38]')
    parser.add_argument('--thread', default=1, type=int, help='thread of annovar annotation. [1]')
    parser.add_argument('-j', '--job', default=1, type=int,
                        help='multiple running to get reference, if set -j > 1, please install parallel. [1]')
    parser.add_argument('-n', '--paranum', default=100000, type=int,
                        help='number of lines that get reference one job, please set -j > 1. Recommend 50k-1m [100000]')
    parser.add_argument('--reserved', action="store_true", help='Whether to reserved temp documents.')
    # 提取参数
    args = parser.parse_args()
    if args.vcf:
        args_dict['vcf'] = args.vcf
    else:
        sys.exit('Error: Please set input file.')
    if args.annofile:
        args_dict['outfile'] = args.annofile
    else:
        args_dict['outfile'] = 'result.matrix'
        print('Warn: Result will be output as < result.avinput >.')
    args_dict['tmp'] = args.tmp
    if os.path.exists(args_dict['tmp']):
        args_dict['tmp'] = args_dict['tmp'].rstrip('/') + '_tmp4v2a' + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        print('Warn: Temp directory is exits, the temp directory will change to %s .' % (args_dict['tmp']))
    os.makedirs(args_dict['tmp'])
    # 其他参数
    args_dict['depth'] = args.dp
    args_dict['refTag'] = args.refTag
    args_dict['thread'] = args.thread
    if not args.humandb:
        sys.exit('Error: Please set humandb directory.')
    args_dict['humandb'] = args.humandb
    if not args.reference:
        sys.exit('Error: Please give reference.')
    args_dict['reference'] = args.reference
    args_dict['job'] = args.job
    args_dict['lines'] = args.paranum
    args_dict['reserved'] = args.reserved

    return args_dict


def vcf2matrix(file, matrix):
    # 文件识别
    if file.endswith('vcf.gz'):
        vcf = gzip.open(file, "rb")
    elif file.endswith('vcf'):
        vcf = open(file)
    else:
        sys.exit('Error: Can not identify the format of < %s > .' % file)
    # 输出表头
    var_matrix = open(matrix, 'w')
    col_h = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
             'AC', 'AF', 'AN', 'BaseQRankSum', 'ExcessHet', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'MQRankSum',
             'ReadPosRankSum', 'DP_1',
             'GT', 'AD', 'DP', 'GQ', 'PL0', 'PL1', 'PL2',
             'start_expend', 'end_expend']
    col = ['NA'] * col_h.__len__()
    var_matrix.write('\t'.join(col_h) + '\n')
    for var in vcf:
        if not var.decode().startswith('#'):
            var_info_list = var.decode().strip().split('\t')
            # 去除 chr，保留 vcf 前8列
            col[0] = var_info_list[0].lstrip('chr')
            col[1:8] = var_info_list[1:8]
            # info
            info = var_info_list[7].split(';')
            info_dict = {}
            for i in info:
                info_dict[i.split('=')[0]] = i.split('=')[1]
            for i in range(8, 21):
                col[i] = info_dict.get(col_h[i], 'NA')
            col[21] = info_dict.get('DP', 'NA')
            # format
            format_ = var_info_list[8].split(':')
            format_dcit = {}
            for i in ['GT', 'AD', 'DP', 'GQ', 'PL']:
                format_dcit[i] = var_info_list[9].split(':')[format_.index(i)]
            for i in range(22, 26):
                col[i] = format_dcit.get(col_h[i], 'NA')
            col[26] = format_dcit['PL'].split(',')[0]
            col[27] = format_dcit['PL'].split(',')[1]
            col[28] = format_dcit['PL'].split(',')[2]
            # 用于碱基注释
            col[29] = str(eval(col[1]) - 250)
            col[30] = str(eval(col[1]) + 250)

            var_matrix.write('\t'.join(col) + '\n')
    #         print(var_info_list)
    #         print(info_dict)
    #         print(col)
    #         break
    var_matrix.close()
    print('Log: Get raw matrix done.')


def matrix2avinput(matrix, matrix_filter, dp_th):
    df = pd.read_table(matrix, low_memory=False)
    df2 = df[(df['DP'] > dp_th) & (df['ALT'] != '*') & df['#CHROM'].isin(
        ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
         '21', '22', 'X', 'Y'])]
    _df = df2['AD'].str.split(',', expand=True)
    _df.columns = ['AD_r', 'AD_f']
    df2 = df2.join(_df)
    df2.loc[((df2['REF'].isin(['A', 'T', 'C', 'G'])) & (df2['ALT'].isin(['A', 'T', 'C', 'G']))), 'vtype'] = 'SNV'
    df2.loc[(~df2['ALT'].isin(['A', 'T', 'C', 'G'])), 'vtype'] = 'INS'
    df2.loc[(~df2['REF'].isin(['A', 'T', 'C', 'G'])), 'vtype'] = 'DEL'
    df2 = df2[(df2['AD_f'] != '0') | (df2['AD_r'] != '0')]
    df2 = df2.apply(pd.to_numeric, errors='ignore')
    df2['VarBalance'] = df2.apply(lambda x: min(x['AD_r'], x['AD_f']) / max(x['AD_f'], x['AD_r']), axis=1)
    df2['Vaf'] = df2.apply(lambda x: x['AD_f'] / x['DP'], axis=1)
    df2['varLen'] = df2.apply(lambda x: max(len(x['REF']), len(x['ALT'])), axis=1)
    df2.index = pd.RangeIndex(len(df2.index))
    df2.to_csv(matrix_filter, sep='\t', index=False)
    print('Log: Filter raw matrix done.')
    return df2


def annotation(mat, matrix_filter, directory, script, ref, db, th):
    anno_vcf = directory + '/anno.vcf'
    anno_avinput = directory + '/var.avinput'
    anno_out = directory + '/var_tmp'
    raw_anno = directory + '/var_tmp.%s_multianno.txt' % ref  # var_tmp.hg38_multianno.txt
    anno_info = directory + '/var_tmp.info'
    final_anno = directory + '/var.anno.tsv'

    mat[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].to_csv(anno_vcf, sep='\t', index=False)
    v2a_cmd = 'perl %s/convert2annovar.pl -format vcf4 %s > %s' % (script, anno_vcf, anno_avinput)
    v2a_log = os.system(v2a_cmd)
    if v2a_log:
        sys.exit('Error: Something wrong with convert vcf to annovar format.')

    ann_cmd = 'perl %s/table_annovar.pl %s %s --buildver %s -out %s -remove ' \
              '-protocol refGene,rmsk,cpgIslandExt,genomicSuperDups -operation g,r,r,r ' \
              '-nastring . --thread %d' % (script, anno_avinput, db, ref, anno_out, th)
    ann_log = os.system(ann_cmd)
    if ann_log:
        sys.exit('Error: Something wrong with annovar annotation.')
    cmb_cmd = 'cut -f 6- %s > %s && paste %s %s > %s' % (raw_anno, anno_info, matrix_filter, anno_info, final_anno)
    cmb_log = os.system(cmb_cmd)
    if cmb_log:
        sys.exit('Error: Something wrong with merge annotation.')
    print('Log: Annotation done.')


def anno_reference(file, seq_matrix, ref, tmp_dir, script, para, num):
    # cut region
    region_file = tmp_dir + '/region.tsv'
    cut_cmd = 'cut -f 1,30,31 %s | tail -n+2 > %s' % (file, region_file)
    cut_log = os.system(cut_cmd)
    if cut_log:
        sys.exit('Error: Something wrong with get reference region.')
    # get ref
    if para > 1:
        os.makedirs(tmp_dir + '/tmp_split/')
        tmp_pr = tmp_dir + '/tmp_split/varSplit'
        split_cmd = 'split -a 2 -d -l %d %s %s' % (num, region_file, tmp_pr)
        split_log = os.system(split_cmd)
        if split_log:
            sys.exit('Error: Something wrong with split file.')
        anno_ref_cmd = 'ls %s | parallel -j %d "bash %s/getRef.sh %s/{} %s " ' % (
            tmp_dir + '/tmp_split/', para, script, tmp_dir + '/tmp_split/', ref)
        anno_ref_log = os.system(anno_ref_cmd)
        if anno_ref_log:
            sys.exit('Error: Something wrong with annotate reference.')
        merge_cmd = 'cat %s*.seq > %s' % (tmp_pr, region_file + '.seq')
        merge_log = os.system(merge_cmd)
        if merge_log:
            sys.exit('Error: Something wrong with merge split file.')
    else:
        anno_ref_cmd = 'bash %s/getRef.sh %s %s' % (script, region_file, ref)
        anno_ref_log = os.system(anno_ref_cmd)
        if anno_ref_log:
            sys.exit('Error: Something wrong with annotate reference.')
    # merge
    ref_out = open(region_file + '.format.seq', 'w')
    ref_out.write('chr\tstart\tend\tSequence\n')
    with open(region_file + '.seq') as f:
        i = f.readline()
        ref_out.write(i.upper().strip())
        for i in f:
            if i[0] in '123456789XY':
                ref_out.write('\n' + i.upper().strip())
            else:
                ref_out.write(i.upper().strip())
    ref_out.write('\n')
    ref_out.close()
    seq = tmp_dir + '/seq.txt'
    # seq_matrix = tmp_dir + '/annoWithSeq.tsv'
    mergeSeq_cmd = 'cut -f 4 %s > %s && paste %s %s > %s' % (region_file + '.format.seq', seq, file, seq, seq_matrix)
    mergeSeq_log = os.system(mergeSeq_cmd)
    if mergeSeq_log:
        sys.exit('Error: Something wrong with merge sequence.')

    print('Log: Get reference done.')


def remove_tmp(tmp_dir, remain):
    if remain:
        return
    else:
        rm_tmp_cmd = 'rm -rf %s' % tmp_dir
        rm_tmp_log = os.system(rm_tmp_cmd)
        if rm_tmp_log:
            print('Error: Something wrong with delete temp file.')
    print('Log: Remove temp file done.')


def main():
    args = get_args()
    vcf2matrix(args['vcf'], args['tmp'] + '/matrix_raw.tsv')
    mat_df = matrix2avinput(args['tmp'] + '/matrix_raw.tsv', args['tmp'] + '/matrix.vcf', args['depth'])
    script = os.path.split(os.path.realpath(__file__))[0] + '/bin'
    annotation(mat_df, args['tmp'] + '/matrix.vcf', args['tmp'], script, args['refTag'], args['humandb'],
               args['thread'])
    anno_reference(args['tmp'] + '/var.anno.tsv', args['outfile'], args['reference'], args['tmp'], script, args['job'],
                   args['lines'])
    remove_tmp(args['tmp'], args['reserved'])
    print('Log: Finish!')


if __name__ == '__main__':
    main()
