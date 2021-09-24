import sys
import numpy as np
import pandas as pd
import tensorflow as tf


def snv_encode(seq, w):
    code = {
        'N': [0, 0, 0, 0],
        'A': [1, 0, 0, 0],
        'T': [0, 1, 0, 0],
        'C': [0, 0, 1, 0],
        'G': [0, 0, 0, 1],
    }
    _mat = []
    if len(seq) != w:
        seq += 'N' * (w - len(seq))
    for r in seq:
        _tmp = [code[r], ]
        _mat.append(_tmp)
    return np.array(_mat, dtype=np.float32).transpose((2, 0, 1))  # hwc: (4,w,11)


def toTFRecord(writer, mat, fetureList, tag):
    features = {
        'matrix': tf.train.Feature(bytes_list=tf.train.BytesList(value=[mat.tobytes()])),
        # 'shape': tf.train.Feature(int64_list=tf.train.Int64List(value=mat.shape)),
        'feature': tf.train.Feature(float_list=tf.train.FloatList(value=fetureList)),
        'tags': tf.train.Feature(int64_list=tf.train.Int64List(value=[tag]))
    }
    tf_features = tf.train.Features(feature=features)
    tf_example = tf.train.Example(features=tf_features)
    tf_serialized = tf_example.SerializeToString()
    # 写入一个序列化的样本
    writer.write(tf_serialized)


file = sys.argv[1]
out_tfrecode = sys.argv[2]
w = 501
writer = tf.io.TFRecordWriter(out_tfrecode)
df = pd.read_feather(file)
col_list = df.columns.to_list()
seq_idx = col_list.index('Sequence') + 1
tag_idx = col_list.index('tag') + 1
features = ['rmsk', 'cpgIslandExt', 'genomicSuperDups', 'VarBalance', 'Vaf', 'DEL', 'INS', 'SNV', 'Het', 'Hom',
            'QUAL_norm', 'FS_norm', 'MQ_norm', 'QD_norm', 'SOR_norm', 'DP_rate', 'GQ_norm', 'PL0_norm', 'PL1_norm',
            'PL2_norm', 'varLen_norm', 'BaseQRankSum_norm', 'MQRankSum_norm', 'ReadPosRankSum_norm']
f_idx_list = []
for f in features:
    f_idx_list.append(col_list.index(f) + 1)

for tup in df.itertuples():
    seq = tup[seq_idx]
    seq_matrix = snv_encode(seq, w)
    tag = tup[tag_idx]
    feature_list = []
    for i in f_idx_list:
        feature_list.append(tup[i])
    toTFRecord(writer, seq_matrix, feature_list, tag)

writer.close()
