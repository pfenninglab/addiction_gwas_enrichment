from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Embedding, Conv2D, MaxPooling2D, Flatten
from tensorflow.keras.optimizers import SGD, Adam
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.regularizers import l2
from tensorflow.keras.models import load_model
from clr_callback import *
from sklearn import metrics
from Bio import SeqIO

import tensorflow as tf
import tensorflow.keras.backend as K
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import math
import os


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# make sure running on GPU #
tf.config.experimental.list_physical_devices('GPU')
# tf.get_logger().setLevel('WARNING')


def onehot_seq(seq):
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((len(seq),4), dtype='int8')
    for idx,letter in enumerate(seq):
        if letter not in ['N','n']:
            to_return[idx,letter_to_index[letter]] = 1
    return to_return


def encode_sequence(fasta_pos, fasta_neg, shuffleOff = True):
    # read in FASTA files and rev. complement for positive and negatives
    x_pos = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta_pos, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_pos, "fasta") ]) 
    x_neg = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta_neg, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_neg, "fasta") ])
    # concatenate positives and negatives
    print(f'There {x_pos.shape[0]} positives and {x_neg.shape[0]} negatives.')
    x = np.expand_dims(np.concatenate((x_pos, x_neg)), axis=3) 
    y = np.concatenate((np.ones(len(x_pos)),np.zeros(len(x_neg))))
    # need to shuffle order of training set for validation splitting last
    if not shuffleOff:
        indices = np.arange(y.shape[0])
        np.random.shuffle(indices)
        x = x[indices,:]
        y = y[indices]
    #
    return x, y


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
    tf.keras.backend.clear_session()
    # initialize CNN model 
    model = Sequential()
    # 5 Convolutional layers
    model.add(Conv2D(filters = args.conv_filters, kernel_size = (args.conv_width,4), activation = 'relu', kernel_regularizer = l2(l=args.l2_reg), input_shape=input_shape))
    model.add(Conv2D(filters = args.conv_filters, kernel_size = (args.conv_width,1), activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Conv2D(filters = args.conv_filters, kernel_size = (args.conv_width,1), activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Conv2D(filters = args.conv_filters, kernel_size = (args.conv_width,1), activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Conv2D(filters = args.conv_filters, kernel_size = (args.conv_width,1), activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    # 1 max pooling layer
    model.add(MaxPooling2D(pool_size=(args.max_pool_size,1), strides=args.max_pool_stride))
    # 1 dense layer 
    model.add(Dense(units = args.dense_filters, activation = 'relu', kernel_regularizer = l2(l=args.l2_reg)))
    model.add(Flatten())
    # dropout to output layer
    model.add(Dropout(rate = args.dropout))
    # output layer
    model.add(Dense(units = 1, activation = 'sigmoid', kernel_regularizer = l2(l=args.l2_reg)))
    # early stopping parameters & optimizer
    myoptimizer = SGD()      
    model.compile(loss = args.mylossfunc , optimizer = myoptimizer, metrics =['accuracy', macro_f1])
    # train model
    # history = model.fit( xtrain, ytrain, batch_size = args.batch_size, epochs = args.epochs,
    # verbose = 1, validation_split = args.validation_split, callbacks = [ early_stopping ])
    return model



def range_test(x, y, args):
    #####################################
    # Determine the "LR range test."
    # An LR range test can be done using the triangular policy; simply set base_lr and max_lr to define 
    # the entire range you wish to test over, and set step_size to be the total number of iterations in 
    # the number of epochs you wish to test on. This linearly/exponentially increases the learning rate 
    # at each iteration over the range desired.
    #
    # scaling function to be traverse log scale between learning rates
    clr_fn = lambda x: 1/(10**(x))
    args.numCycles = .5 # default for range test, only ramp up learning rate
    args.validation_split = 0 # don't need to have validation split for range test
    iterPerEpoch = y.shape[0] * (1- args.validation_split) / args.batch_size # unit is iter
    stepsize = int(args.step_size_scale * iterPerEpoch)
    # Authors suggest setting step_size = (2-8) x (training iterations in epoch)
    epochs = math.ceil(args.numCycles*2*stepsize/iterPerEpoch)
    #
    # set cyclic learning rate w/ log-linear increasing learning rates
    clr = CyclicLR(mode=args.clr_mode, step_size=stepsize, 
        base_lr=args.base_lr, max_lr=args.max_lr,   
        scale_fn = clr_fn, scale_mode=args.scale_mode)
    clr._reset()
    model = get_model(input_shape = x.shape[1:], args = args)
    hist = model.fit(x, y, batch_size = args.batch_size, epochs = epochs,
        verbose = 2,  validation_split = args.validation_split, callbacks = [clr])
    return model, clr




def train_model_clr(x, y, args):
    # dealing w/ class imbalance
    total = y.shape[0]
    weight_for_0 = (1 / np.sum(y==0))*(total)/2.0 
    weight_for_1 = (1 / np.sum(y==1))*(total)/2.0
    class_weight = {0: weight_for_0, 1: weight_for_1}
    # An epoch is calculated by dividing the number of training images by the batchsize
    iterPerEpoch = y.shape[0] * (1- args.validation_split) / args.batch_size # unit is iter
    # number of training iterations per half cycle. 
    # Authors suggest setting step_size = (2-8) x (training iterations in epoch)
    stepsize = int(args.step_size_scale * iterPerEpoch)
    epochs = math.ceil(args.numCycles*2*stepsize/iterPerEpoch)
    #
    # set cyclic learning rate
    clr = CyclicLR(mode=args.clr_mode, step_size=stepsize, 
        base_lr=args.base_lr, max_lr=args.max_lr)
    clr._reset()
    model = get_model(input_shape = x.shape[1:], args = args)
    hist = model.fit(x, y, batch_size = args.batch_size, epochs = epochs, class_weight = class_weight,
        verbose = args.verbose,  validation_split = args.validation_split, callbacks = [clr])
    return model, clr



def predict_sequences(model_name, x, y, args):
    # Creating an empty Dataframe with column names only
    df = pd.DataFrame(columns=[
        'sample','model_class', 'model_type', 'model_species', 
        # model name
        'model', 'celltype',
        # performance metrics
        'accuracy', 'auROC','auPRC', 'f1','fhalf_score'])
    # get model values
    fold = args.label.split('_')[4].split('.')[0]
    celltype = args.label.split('_')[3]
    # predict labels
    model = load_model(model_name, compile=False)
    y_pred_score = model.predict(x)
    y_pred_class = model.predict_classes(x)
    # compute prediction statistics 
    accuracy = metrics.balanced_accuracy_score(y, y_pred_class)
    f1_score = metrics.f1_score(y, y_pred_class, average = 'weighted')
    fhalf_score = metrics.fbeta_score(y, y_pred_class, beta = 0.5, average = 'weighted')
    roc_auc = metrics.roc_auc_score(y, y_pred_score)
    precision, recall, thresholds = metrics.precision_recall_curve(y, y_pred_score)
    prc_auc = metrics.metrics.auc(recall, precision) # x, y
    print(f'Accuracy: {accuracy}. ')
    print(f'f1_score: {f1_score}.')
    print(f'roc_auc: {roc_auc}.')
    print(f'prc_auc: {prc_auc}.')
    # add this row to dataframe
    df = df.append({ 'sample' : 'BICCN_huMOp', 'model': model_name, 
        'model_class': 'CNN', 'celltype': celltype, 'model_type' : fold, 
        'model_species' : 'Human',
        'accuracy': accuracy, 'auROC': roc_auc, 'auPRC': prc_auc, 
        'f1_score': f1_score, 'fhalf_score': fhalf_score}, ignore_index=True)
    return df


def main(args):
    """Main function
    Args:
        args (argparse):
    """
    # file names
    print(f'Running cyclical learning rate for {args.label}.')
    print(f'Using {args.clr_mode} CLR mode and for {args.numCycles} cycles.')
    print(f'Learning rates range: {args.base_lr} - {args.max_lr}.')
    model_range_test = f'plots/{args.label}/{args.label}_clrRangeTest_SS{args.step_size_scale}_NC{args.numCycles}_BR{args.base_lr}_MR{args.max_lr}.pdf'
    model_name = f"cnn/singleTask/{args.label}_CM{args.clr_mode}_SS{args.step_size_scale}_NC{args.numCycles}_BR{args.base_lr}_MR{args.max_lr}.h5"
    model_train_performance = f'plots/{args.label}/{args.label}_trainingPerformance_CM{args.clr_mode}_SS{args.step_size_scale}_NC{args.numCycles}_BR{args.base_lr}_MR{args.max_lr}.pdf'
    model_test_performance = f'predictions/cnn/{args.label}/{args.label}_testPerformance_CM{args.clr_mode}_SS{args.step_size_scale}_NC{args.numCycles}_BR{args.base_lr}_MR{args.max_lr}.feather'
    # read in the training sequences
    (x, y) = encode_sequence(args.fasta_pos, args.fasta_neg, shuffleOff = args.predict_mode)
   # call main functions
    if args.range_test:
        print('In range test mode. Will only ramp up learning rate (1/2 a cycle).')
        model, clr = range_test(x, y, args)
            #
        if not os.path.exists(f'plots/{args.label}'):
            os.makedirs(f'plots/{args.label}')
        # plot accuracy by learning rate
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.0))
        fig.suptitle(f"CLR Range test for {args.label} with 'log10-triangular' Policy")
        ax1.plot(clr.history['lr'], clr.history['loss'])
        ax2.plot(clr.history['lr'], clr.history['accuracy'])
        ax1.set_xscale('log'); ax2.set_xscale('log')
        ax1.set(xlabel='Learning Rate', ylabel='Loss')
        ax2.set(xlabel='Learning Rate', ylabel='Accuracy')
        fig.savefig(model_range_test, bbox_inches='tight')
    elif not args.predict_mode:
        print('In training mode.')
        model, clr = train_model_clr(x, y, args)
        # make sure to create folder
        if not os.path.exists(f'cnn/singleTask/{args.label}'):
            os.makedirs(f'cnn/singleTask/{args.label}')
        model.save(model_name)
        #
        # plot performance by iteration
        if not os.path.exists(f'plots/{args.label}'):
            os.makedirs(f'plots/{args.label}')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.0))
        fig.suptitle(f"Training performance {args.label}, {args.clr_mode} Policy")
        ax1.plot(clr.history['iterations'], clr.history['loss'])
        ax2.plot(clr.history['iterations'], clr.history['accuracy'])
        ax1.set(xlabel='Iterations', ylabel='Loss')
        ax2.set(xlabel='Iterations', ylabel='Accuracy')
        fig.savefig(model_train_performance, bbox_inches='tight')
    else:
        print('In predict mode.')
        df = predict_sequences(model_name, x, y, args)
        # save model performances to feather object
        if not os.path.exists(f'predictions/cnn/{args.label}'):
            os.makedirs(f'predictions/cnn/{args.label}')
        df.to_feather(model_test_performance)
    #
    return



if __name__ == '__main__':  
    #
    # set cnn parameters:
    parser = argparse.ArgumentParser(description='Parse CNN parameters.')
    parser.add_argument("label", type=str, help="label of model.")
    parser.add_argument("fasta_pos", type=str, help="fasta sequence file of positives.")
    parser.add_argument("fasta_neg", type=str, help="fasta sequence file of negatives.")
    #
    # set cnn parameters:
    parser.add_argument("--conv_width", type=int, help="Convolution width.", default = 8, required=False)
    parser.add_argument("--conv_filters", type=int, help="Convolution filter width.", default = 200, required=False)
    parser.add_argument("--max_pool_size", type=int, help="Max pool size.", default = 26, required=False)
    parser.add_argument("--max_pool_stride", type=int, help="Max pool stride.", default = 26, required=False)
    parser.add_argument("--dense_filters", type=int, help="Number of dense filters.", default = 300, required=False)
    parser.add_argument("--l1_reg", type=float, help="L1 regularization rate.", default = 0, required=False)
    parser.add_argument("--l2_reg", type=float, help="L2 regularization rate.", default = 1e-5, required=False)
    parser.add_argument("--dropout", type=float, help="Percent dropout.", default = .1, required=False)
    parser.add_argument("--early_stop_num", type=int, help="Number of epochs to stop.", default = 10, required=False)
    parser.add_argument("--batch_size", type=int, help="Batch size.", default = 200, required=False)
    parser.add_argument("--verbose", type=int, help="keras verbosity", default = 2, required=False)
    parser.add_argument("--validation_split", help="Percent validation split.",
        type=float, default = 0.1, required=False)
    parser.add_argument("--mylossfunc", help="Loss function.", default = 'binary_crossentropy',
        choices=['binary_crossentropy', 'categorical_crossentropy', 'sparse_categorical_crossentropy'], required=False)
    #
    # parameters for cyclical learning Rate, see https://github.com/bckenstler/CLR for details
    parser.add_argument("--numCycles", help="Number of cyclical learning rate cycles.", 
        type=int, default = 8, required=False)
    parser.add_argument("--base_lr", help="Learning rate floor value for cyclical learning rate model.", 
        type=float, default = 2e-4, required=False)
    parser.add_argument("--max_lr", help="Learning rate ceiling value for cyclical learning rate model.", 
        type=float,default = 2e-1, required=False)
    parser.add_argument("--gamma", help="Exp_range degradation rate.", 
        type=float, default = 0.99998, required=False)
    parser.add_argument("--step_size_scale", help="Number to scale stepsize by [2-10].", 
        type=float, default = 4, required=False)
    parser.add_argument("--clr_mode", type=str, help="Learning rate cycle mode.", 
        choices=['triangular', 'triangular2', 'exp_range'], default = 'triangular', required=False)
    parser.add_argument("--scale_mode", type=str, help="Learning rate cycle mode.", 
        choices=['cycle', 'iterations'], default = 'cycle', required=False)
    parser.add_argument("--range_test",help="Compute range test between base and max learning rates.",
        action='store_true')
    parser.add_argument("--predict_mode",help="Predict on test sequences model instead of train.",
        action='store_true')
    # parse arguments
    args = parser.parse_args()
    # args = parser.parse_args(['BICCN_huMOp_SNARE-Seq2_PVALB_fold0', 
    #     'fasta/BICCN_huMOp_SNARE-Seq2_PVALB_fold0_trainPos.fa', 
    #     'fasta/BICCN_huMOp_SNARE-Seq2_PVALB_fold0_trainNeg.fa'])
    # args = parser.parse_args(['BICCN_huMOp_SNARE-Seq2_PVALB_fold0', 
    #     'fasta/BICCN_huMOp_SNARE-Seq2_PVALB_fold0_testPos.fa', 
    #     'fasta/BICCN_huMOp_SNARE-Seq2_PVALB_fold0_testNeg.fa',
    #     '--predict_mode'])
    main(args)


