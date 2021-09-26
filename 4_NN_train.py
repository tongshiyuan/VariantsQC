import os
import sys
import math
import datetime
from tqdm import tqdm
import tensorflow as tf

NUM_CLASSES = 1
bs = 1024
initial_lr = 0.0001
epochs = 100


# os.environ["CUDA_VISIBLE_DEVICES"]="1"

class BasicBlock(tf.keras.layers.Layer):
    def __init__(self, filter_num, stride=1):
        super(BasicBlock, self).__init__()
        self.conv1 = tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(1, 1), strides=stride, padding="same",
                                            kernel_initializer='he_uniform')
        self.bn1 = tf.keras.layers.BatchNormalization()
        self.conv2 = tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(1, 1), strides=1, padding="same",
                                            kernel_initializer='he_uniform')
        self.bn2 = tf.keras.layers.BatchNormalization()
        if stride != 1:
            self.downsample = tf.keras.Sequential()
            self.downsample.add(tf.keras.layers.Conv2D(filters=filter_num, kernel_size=(1, 1), strides=stride,
                                                       kernel_initializer='he_uniform'))
            self.downsample.add(tf.keras.layers.BatchNormalization())
        else:
            self.downsample = lambda x: x

    def call(self, inputs, training=None, **kwargs):
        residual = self.downsample(inputs)

        x = self.conv1(inputs)
        x = self.bn1(x, training=training)
        x = tf.nn.relu(x)
        x = self.conv2(x)
        x = self.bn2(x, training=training)

        output = tf.nn.relu(tf.keras.layers.add([residual, x]))

        return output


def make_basic_block_layer(filter_num, blocks, stride=1):
    res_block = tf.keras.Sequential()
    res_block.add(BasicBlock(filter_num, stride=stride))

    for _ in range(1, blocks):
        res_block.add(BasicBlock(filter_num, stride=1))

    return res_block


cov_input = tf.keras.layers.Input(shape=(4, 501, 1))
var_input = tf.keras.layers.Input(shape=(24,))

# cov
conv1 = tf.keras.layers.Conv2D(filters=64,
                               kernel_size=(3, 3), strides=2,
                               activation='relu', padding="same",
                               kernel_initializer='he_uniform')(cov_input)
conv1 = tf.keras.layers.BatchNormalization()(conv1)
conv1 = make_basic_block_layer(filter_num=64, blocks=3)(conv1)
conv1 = make_basic_block_layer(filter_num=128, blocks=4, stride=2)(conv1)
conv1 = make_basic_block_layer(filter_num=256, blocks=6, stride=2)(conv1)
conv1 = make_basic_block_layer(filter_num=512, blocks=3, stride=2)(conv1)

conv1 = tf.keras.layers.GlobalAveragePooling2D()(conv1)
conv1 = tf.keras.layers.Dense(units=18, activation=tf.keras.activations.sigmoid, kernel_initializer='he_uniform')(conv1)
# nocov
# var1 = tf.keras.layers.GlobalAveragePooling2D()(var_input)
var1 = tf.keras.layers.Dense(36, activation='relu', kernel_initializer='he_uniform')(var_input)
var1 = tf.keras.layers.Dense(18, activation='relu', kernel_initializer='he_uniform')(var1)
var1 = tf.keras.layers.Dense(18, activation='relu', kernel_initializer='he_uniform')(var1)
# sum
concat = tf.keras.layers.concatenate([conv1, var1])
out = tf.keras.layers.Dense(36, activation='relu', kernel_initializer='he_uniform')(concat)
out = tf.keras.layers.Dense(36, activation='relu', kernel_initializer='he_uniform')(out)
out = tf.keras.layers.Dense(18, activation='relu', kernel_initializer='he_uniform')(out)
out = tf.keras.layers.Dense(1, activation='sigmoid')(out)

sum_model = tf.keras.models.Model(inputs=[cov_input, var_input], outputs=out)
sum_model.summary()


@tf.function
def train_step(train_conv, train_var, train_labels):
    with tf.GradientTape() as tape:
        output = sum_model([train_conv, train_var], training=True)
        loss = loss_object(train_labels, output)
    gradients = tape.gradient(loss, sum_model.trainable_variables)
    optimizer.apply_gradients(zip(gradients, sum_model.trainable_variables))

    train_loss(loss)
    train_accuracy(train_labels, output)
    train_precision(train_labels, output)
    train_recall(train_labels, output)
    train_auc(train_labels, output)


@tf.function
def val_step(val_conv, val_var, val_labels):
    output = sum_model([val_conv, val_var], training=False)
    loss = loss_object(val_labels, output)

    val_loss(loss)
    val_accuracy(val_labels, output)
    val_precision(val_labels, output)
    val_recall(val_labels, output)
    val_auc(val_labels, output)


def scheduler(now_epoch):
    end_lr_rate = 0.01  # end_lr = initial_lr * end_lr_rate
    rate = ((1 + math.cos(now_epoch * math.pi / epochs)) / 2) * (1 - end_lr_rate) + end_lr_rate  # cosine
    new_lr = rate * initial_lr
    # writing lr into tensorboard
    with train_writer.as_default():
        tf.summary.scalar('learning rate', data=new_lr, step=epoch)

    return new_lr


# using keras low level api for training
# loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
loss_object = tf.keras.losses.BinaryCrossentropy(from_logits=False)
# optimizer = tf.keras.optimizers.SGD(learning_rate=initial_lr, momentum=0.9)
optimizer = tf.keras.optimizers.Adam(learning_rate=initial_lr)

train_loss = tf.keras.metrics.Mean(name='train_loss')
# train_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='train_accuracy')
train_accuracy = tf.keras.metrics.BinaryAccuracy(name='train_accuracy')
train_precision = tf.keras.metrics.Precision(name='train_precision')
train_recall = tf.keras.metrics.Recall(name='train_recall')
train_auc = tf.keras.metrics.AUC(name='train_auc')

val_loss = tf.keras.metrics.Mean(name='val_loss')
# val_accuracy = tf.keras.metrics.SparseCategoricalAccuracy(name='val_accuracy')
val_accuracy = tf.keras.metrics.BinaryAccuracy(name='val_accuracy')
val_precision = tf.keras.metrics.Precision(name='val_precision')
val_recall = tf.keras.metrics.Recall(name='val_recall')
val_auc = tf.keras.metrics.AUC(name='val_auc')


def getTFRecord(x):
    shape = [4, 501, 1]
    feature = {
        'matrix': tf.io.FixedLenFeature(shape=[], dtype=tf.string),
        #         'shape': tf.io.FixedLenFeature(shape=[3], dtype=tf.int64),
        'feature': tf.io.FixedLenFeature(shape=[24], dtype=tf.float32),
        'tags': tf.io.FixedLenFeature(shape=[1], dtype=tf.int64)
    }
    parsed_example = tf.io.parse_single_example(x, feature)
    #     parsed_example['matrix'] = tf.sparse_tensor_to_dense(parsed_example['matrix'])
    parsed_example['matrix'] = tf.io.decode_raw(parsed_example['matrix'], tf.float32)
    parsed_example['matrix'] = tf.reshape(parsed_example['matrix'], shape)
    #     seq1 = tf.expand_dims(parsed_example['matrix'][:,:,0],2)
    #     seq2 = parsed_example['matrix'][:,:,1:10]
    #     seq3 = tf.expand_dims(parsed_example['matrix'][:,:,10],2) * seq1
    #     mat = tf.concat([seq1, seq2, seq3], axis = 2)
    mat = parsed_example['matrix']
    tag = parsed_example['tags']
    fea = parsed_example['feature']
    return mat, tag, fea


reader1 = tf.data.TFRecordDataset('process/matrix/HG001U0.tfrecord')
reader1 = reader1.map(getTFRecord)
train_ = reader1.batch(batch_size=bs)

reader2 = tf.data.TFRecordDataset('process/matrix/HG001U1.tfrecord')
reader2 = reader2.map(getTFRecord)
test_ = reader2.batch(batch_size=bs)

log_dir = "./logs/model_" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
train_writer = tf.summary.create_file_writer(os.path.join(log_dir, "train"))
val_writer = tf.summary.create_file_writer(os.path.join(log_dir, "val"))

best_val_recall = 0.
for epoch in range(epochs):
    train_loss.reset_states()  # clear history info
    train_accuracy.reset_states()  # clear history info
    train_precision.reset_states()
    train_recall.reset_states()
    val_loss.reset_states()  # clear history info
    val_accuracy.reset_states()  # clear history info
    val_precision.reset_states()
    val_recall.reset_states()

    # train
    train_bar = tqdm(train_)
    for mat, labels, feature in train_bar:
        # fea
        var = tf.constant(feature, dtype=tf.float32)
        # conv
        conv = mat

        train_step(conv, var, labels)
        # print train process
        train_bar.desc = "train epoch[{}/{}] loss:{:.3f}, acc:{:.3f}, precision:{:.3f}, recall:{:.3f}, auc:{:.3f}".format(
            epoch + 1, epochs, train_loss.result(),
            train_accuracy.result(), train_precision.result(), train_recall.result(), train_auc.result())

    # update learning rate
    optimizer.learning_rate = scheduler(epoch)

    # validate
    val_bar = tqdm(test_)
    for mat, labels, feature in val_bar:
        # fea
        var = tf.constant(feature, dtype=tf.float32)
        # conv
        conv = mat

        val_step(conv, var, labels)
        # print val process
        val_bar.desc = "valid epoch[{}/{}] loss:{:.3f}, acc:{:.3f}, precision:{:.3f}, recall:{:.3f}, auc:{:.3f}".format(
            epoch + 1, epochs, val_loss.result(),
            val_accuracy.result(), val_precision.result(), val_recall.result(), val_auc.result())

    # writing training loss and acc
    with train_writer.as_default():
        tf.summary.scalar("loss", train_loss.result(), epoch)
        tf.summary.scalar("accuracy", train_accuracy.result(), epoch)
        tf.summary.scalar("precision", train_precision.result(), epoch)
        tf.summary.scalar("recall", train_recall.result(), epoch)
        tf.summary.scalar("auc", train_auc.result(), epoch)

    # writing validation loss and acc
    with val_writer.as_default():
        tf.summary.scalar("loss", val_loss.result(), epoch)
        tf.summary.scalar("accuracy", val_accuracy.result(), epoch)
        tf.summary.scalar("precision", val_precision.result(), epoch)
        tf.summary.scalar("recall", val_recall.result(), epoch)
        tf.summary.scalar("auc", val_auc.result(), epoch)

    # only save best weights
    if val_recall.result() > best_val_recall:
        best_val_recall = val_recall.result()
        save_name = "./models/model/multi_NN.ckpt"
        sum_model.save_weights(save_name, save_format="tf")
