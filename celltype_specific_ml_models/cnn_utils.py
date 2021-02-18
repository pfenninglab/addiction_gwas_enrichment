from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Embedding, Conv1D, MaxPooling1D, Flatten
from tensorflow.keras.optimizers import SGD, Adam
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2
from tensorflow.keras.models import load_model
from sklearn import metrics
from Bio import SeqIO

import tensorflow as tf
import tensorflow.keras.backend as K
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import os, sys, gc, math, argparse
import pandas as pd
import numpy as np


# gc.collect()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# make sure running on GPU #
tf.config.experimental.list_physical_devices('GPU')
# tf.get_logger().setLevel('WARNING')

def onehot_seq(seq, size = 501):
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((size,4), dtype='uint8')
    # cap length with size
    for idx,letter in enumerate(seq):
        if letter not in ['N','n'] and idx < size:
            to_return[idx,letter_to_index[letter]] = 1
    return to_return


def encode_sequence(fasta_pos, fasta_neg, size, shuffleOff = True):
    x_pos = np.array([onehot_seq(seq, size) for seq in SeqIO.parse(fasta_pos, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_pos, "fasta") ], dtype="uint8") 
    x_neg = np.array([onehot_seq(seq, size) for seq in SeqIO.parse(fasta_neg, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_neg, "fasta") ], dtype="uint8")
    # concatenate positives and negatives, make sure names are unique
    ids = np.array( [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_pos, "fasta")) ] +
            [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_pos, "fasta")) ] +
            [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_neg, "fasta")) ] +
            [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_neg, "fasta")) ], dtype="object")
    print(f'There {x_pos.shape[0]} positives and {x_neg.shape[0]} negatives.')
    x = np.concatenate((x_pos, x_neg))
    y = np.concatenate((np.ones(len(x_pos)),np.zeros(len(x_neg)))).astype('uint8')
    # need to shuffle order of training set for validation splitting last
    if not shuffleOff:
        indices = np.arange(y.shape[0])
        np.random.shuffle(indices)
        x = x[indices,:]
        y = y[indices]
        ids = ids[indices]
    #
    return x, y, ids


def encode_sequence2(fasta_file, label_file, size, shuffleOff = True):
    ## read in the label, the fasta ID, and fasta sequences
    ## ids used to combine rev complements of same DNA sequences
    y = np.tile(np.loadtxt(label_file), (2, 1))
    print(f'There {sum(np.sum(y, axis = 1) > 0)} positives and {sum(np.sum(y, axis = 1) == 0)} negatives.')
    ids = np.array([str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_file, "fasta")) ] +
            [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_file, "fasta")) ], dtype="object")
    x = np.array([onehot_seq(seq, size) for seq in SeqIO.parse(fasta_file, "fasta") ] +
    [onehot_seq(seq.reverse_complement(), size) for seq in SeqIO.parse(fasta_file, "fasta") ])
    # x = np.expand_dims(x, axis=3)
    # need to shuffle order of training set for validation splitting last
    if not shuffleOff:
        indices = np.arange(y.shape[0])
        np.random.shuffle(indices)
        x = x[indices,:]
        y = y[indices]
        ids = ids[indices]
    #
    return x, y, ids



def encode_sequence3(fasta_file, size, shuffleOff = True):
    ## read in the fasta ID, and fasta sequences
    ## ids used to combine rev complements of same DNA sequences
    ids = np.array( [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_file, "fasta")) ] +
            [str(idx) + '_' + seq.id for idx, seq in enumerate(SeqIO.parse(fasta_file, "fasta")) ], dtype="object")
    print(f'There {len(ids)} sequences.')
    x = np.array([onehot_seq(seq, size) for seq in SeqIO.parse(fasta_file, "fasta") ] +
    [onehot_seq(seq.reverse_complement(), size) for seq in SeqIO.parse(fasta_file, "fasta") ])
    # x = np.expand_dims(x, axis=3)
    return x, ids



def macro_f1(y, y_hat, thresh=0.5):
    """Compute the macro F1-score on a batch of observations (average F1 across labels)
    Args:
        y (int32 Tensor): labels array of shape (BATCH_SIZE, N_LABELS)
        y_hat (float32 Tensor): probability matrix from forward propagation of shape (BATCH_SIZE, N_LABELS)
        thresh: probability value above which we predict positive
        
    Returns:
        macro_f1 (scalar Tensor): value of macro F1 for the batch
    """
    y_pred = tf.cast(tf.greater(y_hat, thresh), tf.float32)
    tp = tf.cast(tf.math.count_nonzero(y_pred * y, axis=0), tf.float32)
    fp = tf.cast(tf.math.count_nonzero(y_pred * (1 - y), axis=0), tf.float32)
    fn = tf.cast(tf.math.count_nonzero((1 - y_pred) * y, axis=0), tf.float32)
    f1 = 2*tp / (2*tp + fn + fp + 1e-16)
    macro_f1 = tf.reduce_mean(f1, axis=-1)
    return macro_f1



def macro_soft_f1(y, y_hat):
    """Compute the macro soft F1-score as a cost.
    Average (1 - soft-F1) across all labels.
    Use probability values instead of binary predictions.
    Args:
        y (int32 Tensor): targets array of shape (BATCH_SIZE, N_LABELS)
        y_hat (float32 Tensor): probability matrix of shape (BATCH_SIZE, N_LABELS)
    Returns:
        cost (scalar Tensor): value of the cost function for the batch
    """
    y = tf.cast(y, tf.float32)
    y_hat = tf.cast(y_hat, tf.float32)
    tp = tf.reduce_sum(y_hat * y, axis=0)
    fp = tf.reduce_sum(y_hat * (1 - y), axis=0)
    fn = tf.reduce_sum((1 - y_hat) * y, axis=0)
    soft_f1 = 2*tp / (2*tp + fn + fp + 1e-16)
    cost = 1 - soft_f1 # reduce 1 - soft-f1 in order to increase soft-f1
    macro_cost = tf.reduce_mean(cost, axis=-1) # average on all labels
    return macro_cost



def macro_double_soft_f1(y, y_hat):
    """Compute the macro soft F1-score as a cost (average 1 - soft-F1 across all labels).
    Use probability values instead of binary predictions.
    This version uses the computation of soft-F1 for both positive and negative class for each label.
    
    Args:
        y (int32 Tensor): targets array of shape (BATCH_SIZE, N_LABELS)
        y_hat (float32 Tensor): probability matrix from forward propagation of shape (BATCH_SIZE, N_LABELS)
        
    Returns:
        cost (scalar Tensor): value of the cost function for the batch
    """
    y = tf.cast(y, tf.float32)
    y_hat = tf.cast(y_hat, tf.float32)
    tp = tf.reduce_sum(y_hat * y, axis=0)
    fp = tf.reduce_sum(y_hat * (1 - y), axis=0)
    fn = tf.reduce_sum((1 - y_hat) * y, axis=0)
    tn = tf.reduce_sum((1 - y_hat) * (1 - y), axis=0)
    soft_f1_class1 = 2*tp / (2*tp + fn + fp + 1e-16)
    soft_f1_class0 = 2*tn / (2*tn + fn + fp + 1e-16)
    cost_class1 = 1 - soft_f1_class1 # reduce 1 - soft-f1_class1 in order to increase soft-f1 on class 1
    cost_class0 = 1 - soft_f1_class0 # reduce 1 - soft-f1_class0 in order to increase soft-f1 on class 0
    cost = 0.5 * (cost_class1 + cost_class0) # take into account both class 1 and class 0
    macro_cost = tf.reduce_mean(cost, axis=-1) # average on all labels
    return macro_cost



def get_model(input_shape, args):
    """Generate and train a binary CNN classifier on DNA sequences and binary labels.
    Use one-hot encoded DNA sequences of 501 bp long to predict 1 for positives and 0 for negatives.
    This version takes in arguments defined by the arg parsers to define defaults.
    
    Args:
        args (arg parse): arguments of the CNN
        
    Returns:
        model (Keras model): trained Keras model
    """
    ##############################
    # get weights for positive and negative sets
    K.clear_session()
    # initialize CNN model 
    model = Sequential()
    # 5 Convolutional layers, alternating w/ Drop-out
    model.add(Conv1D(filters = args.conv_filters, kernel_size = (args.conv_width), 
        activation = 'relu', kernel_regularizer = l2(l=args.l2_reg), input_shape=input_shape))
    model.add(Dropout(rate = args.dropout))
    model.add(Conv1D(filters = args.conv_filters, kernel_size = (args.conv_width), 
        activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Dropout(rate = args.dropout))
    model.add(Conv1D(filters = args.conv_filters, kernel_size = (args.conv_width), 
        activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Dropout(rate = args.dropout))
    model.add(Conv1D(filters = args.conv_filters, kernel_size = (args.conv_width), 
        activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Dropout(rate = args.dropout))
    model.add(Conv1D(filters = args.conv_filters, kernel_size = (args.conv_width), 
        activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    # 1 max pooling layer
    model.add(MaxPooling1D(pool_size=args.max_pool_size, strides=args.max_pool_stride))
    # 1 dense layer 
    model.add(Dense(units = args.dense_filters, activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Flatten())
    # dropout to output layer
    model.add(Dropout(rate = args.dropout))
    # output layer
    model.add(Dense(units = 1, activation = 'sigmoid', kernel_regularizer = l2(l=args.l2_reg)))
    # early stopping parameters & optimizer
    myoptimizer = SGD(lr=args.base_lr, momentum=args.max_m)      
    model.compile(loss = args.mylossfunc , optimizer = myoptimizer, metrics =['accuracy', macro_f1])
    return model


