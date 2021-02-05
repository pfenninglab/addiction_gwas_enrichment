from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Embedding, Conv2D, MaxPooling2D, Flatten
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

from callback_ocp import *
from cnn_utils_v2 import *

# gc.collect()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# make sure running on GPU #
tf.config.experimental.list_physical_devices('GPU')
# tf.get_logger().setLevel('WARNING')

def range_test(x, y, args):
    #####################################
    # Determine the "LR range test."
    # An LR range test can be done using the triangular policy; simply set base_lr and max_lr to define 
    # the entire range you wish to test over, and set step_size to be the total number of iterations in 
    # the number of epochs you wish to test on. This linearly/exponentially increases the learning rate 
    # at each iteration over the range desired.
    #
    # scaling function to be traverse log scale between learning rates
    args.numCycles = 1
    args.epochs = 4
    # Authors suggest setting step_size = (2-8) x (training iterations in epoch)
    iterations = list(range(0,round(y_train.shape[0]/args.batch_size*args.epochs)+1))
    step_size = len(iterations)/(args.numCycles)
    #
    # set cyclic learning rate
    clr =  CyclicLR(base_lr=args.base_lr,
                max_lr=args.max_lr,
                step_size=step_size,
                max_m=args.max_m,
                base_m=args.base_m,
                cyclical_momentum=args.cyclical_momentum)
    #
    model = get_model(input_shape = x.shape[1:], args = args)
    hist = model.fit(x, y, batch_size = args.batch_size, epochs = epochs,
        verbose = 2,  validation_split = args.validation_split, callbacks = [clr])
    return model, clr

def train_model_clr(x_train, y_train, x_valid, y_valid,  args):
    # dealing w/ class imbalance
    total = y_train.shape[0]
    weight_for_0 = (1 / np.sum(y_train==0))*(total)/2.0 
    weight_for_1 = (1 / np.sum(y_train==1))*(total)/2.0
    class_weight = {0: weight_for_0, 1: weight_for_1}
    # An epoch is calculated by dividing the number of training images by the batchsize
    iterPerEpoch = y_train.shape[0] / args.batch_size # unit is iter
    # number of training iterations per half cycle. 
    # Authors suggest setting step_size = (2-8) x (training iterations in epoch)
    iterations = list(range(0,round(y_train.shape[0]/args.batch_size*args.epochs)+1))
    step_size = len(iterations)/(args.numCycles)
    #
    # set cyclic learning rate
    clr =  CyclicLR(base_lr=args.base_lr,
                max_lr=args.max_lr,
                step_size=step_size,
                max_m=args.max_m,
                base_m=args.base_m,
                cyclical_momentum=args.cyclical_momentum)
    model = get_model(input_shape = x_train.shape[1:], args = args)
    hist = model.fit(x_train, y_train, batch_size = args.batch_size, epochs = args.epochs, 
                    verbose = args.verbose, class_weight = class_weight, 
                    validation_data=(x_valid, y_valid), callbacks = [clr])
    return model, clr, hist




def predict_sequences(model_name, x):
    # Creating an empty Dataframe with column names only
    # predict labels
    model = load_model(model_name, compile=False)
    y_pred_score = model.predict(x)
    y_pred_class = (y_pred_score > 0.5).astype("int32")
    return y_pred_class, y_pred_score



def plot_training_performance(args, hist, clr ):
    """plot function of training history and learning rate ramps
    Args:
        args (argparse):
        hist (keras training history):
        clr (argparse):
    """
    # file names
    lab = f'{args.label}_OCP_NB{args.batch_size}_NE{args.epochs}_BR{args.base_lr}_MR{args.max_lr}_BM{args.base_m}_MM{args.max_m}'
    model_train_performance = f'plots/{args.label}/{lab}_trainingPerformance.pdf'
    # plot performance by iteration
    if not os.path.exists(f'plots/{args.label}'):
        os.makedirs(f'plots/{args.label}')
    fig, axs = plt.subplots(2, 2, figsize=(8, 6))
    fig.suptitle(f"Training performance {args.label}")
    axs[0, 0].plot(clr.history['iterations'], clr.history['lr'])
    axs[0, 0].set(xlabel='Training Iterations', ylabel='Learning Rate')
    axs[0, 0].set_title("One Cycle Policy")
    #
    axs[0, 1].plot(clr.history['iterations'], clr.history['momentum'])
    axs[0, 1].set(xlabel='Training Iterations', ylabel='Momentum')
    axs[0, 1].set_title("One Cycle Policy")
    #
    val_loss = hist.history['val_loss']
    loss = hist.history['loss']
    axs[1, 0].plot(range(len(val_loss)),val_loss,'c',label='Validation loss')
    axs[1, 0].plot(range(len(loss)),loss,'m',label='Train loss')
    axs[1, 0].set(xlabel='Epoch', ylabel='Loss')
    axs[1, 0].legend()
    #
    val_f1 = hist.history['val_macro_f1']
    f1 = hist.history['macro_f1']
    axs[1, 1].plot(range(len(val_f1)),val_f1,'c',label='Validation F1 score')
    axs[1, 1].plot(range(len(f1)),f1,'m',label='Train  F1 score')
    axs[1, 1].set(xlabel='Epoch', ylabel='F1 score')
    axs[1, 1].legend()
    fig.savefig(model_train_performance, bbox_inches='tight')
    return


def main(args):
    """Main function
    Args:
        args (argparse):
    """
    # file names
    print(f'Running cyclical learning rate for {args.label}.')
    print(f'One-cycle policy training for {args.epochs} epochs.')
    print(f'Learning rates range: {args.base_lr} - {args.max_lr}.')
    print(f'Momentum rates range: {args.base_m} - {args.max_m}.')
    lab = f'{args.label}_OCP_NB{args.batch_size}_NE{args.epochs}_BR{args.base_lr}_MR{args.max_lr}_BM{args.base_m}_MM{args.max_m}_DO{args.dropout}'
    model_name = f"cnn_v2/singleTask/{args.label}/{lab}.h5"
    # call main functions
    if args.mode == 'train':
        print('In training mode.')
        (x_train, y_train) = encode_sequence(args.train_fasta_pos, args.train_fasta_neg, shuffleOff = False)
        (x_valid, y_valid) = encode_sequence(args.valid_fasta_pos, args.valid_fasta_neg, shuffleOff = False)
        #
        model, clr, hist = train_model_clr(x_train, y_train, x_valid, y_valid, args)
        plot_training_performance(args, hist, clr )
        # make sure to create folder
        if not os.path.exists(f'cnn_v2/singleTask/{args.label}'):
            os.makedirs(f'cnn_v2/singleTask/{args.label}')
        #
        model.save(model_name)
        #
    elif args.mode == 'evaluate':
        print('In evaluation mode.')
        (x_valid, y_valid) = encode_sequence(args.valid_fasta_pos, args.valid_fasta_neg, shuffleOff = True)
        df = evaluate_sequences(model_name, x_valid, y_valid, args)
        model_test_performance = f'predictions/cnn_v2/{args.label}/{lab}_testPerformance.feather'
        # save model performances to feather object
        if not os.path.exists(f'predictions/cnn_v2/{args.label}'):
            os.makedirs(f'predictions/cnn_v2/{args.label}')
        df.to_feather(model_test_performance)
        #
    elif args.mode == 'predict':
        print('In evaluation mode.')
        x = encode_sequence2(args.valid_fasta_pos)
        y_pred_class, y_pred_score = predict_sequences(model_name, x)
        model_predictions = f'predictions/cnn_v2/{args.label}/{lab}_testPredictions.txt'
        # save model performances to feather object
        if not os.path.exists(f'predictions/cnn_v2/{args.label}'):
            os.makedirs(f'predictions/cnn_v2/{args.label}')
        np.savetxt(model_predictions, np.concatenate([y_pred_class, y_pred_score], axis = 1))
    return



if __name__ == '__main__':  
    # set cnn parameters:
    parser = argparse.ArgumentParser(description='Parse CNN parameters.')
    parser.add_argument("--mode", type=str, help="Mode to perform", 
        choices=['train', 'evaluate', 'predict'], default = 'train', required=False)
    #
    parser.add_argument("label", type=str, help="label of model.")
    parser.add_argument("train_fasta_pos", type=str, help="training fasta sequence file of positives.")
    parser.add_argument("train_fasta_neg", type=str, help="training fasta sequence file of negatives.")
    parser.add_argument("valid_fasta_pos", type=str, help="validation fasta sequence file of positives.")
    parser.add_argument("valid_fasta_neg", type=str, help="validation fasta sequence file of negatives.")
    #
    # set cnn parameters:
    parser.add_argument("--conv_width", type=int, help="Convolution width.", default = 8, required=False)
    parser.add_argument("--conv_filters", type=int, help="Convolution filter width.", default = 200, required=False)
    parser.add_argument("--max_pool_size", type=int, help="Max pool size.", default = 26, required=False)
    parser.add_argument("--max_pool_stride", type=int, help="Max pool stride.", default = 26, required=False)
    parser.add_argument("--dense_filters", type=int, help="Number of dense filters.", default = 300, required=False)
    parser.add_argument("--l2_reg", type=float, help="L2 regularization rate.", default = 1e-10, required=False)
    parser.add_argument("--dropout", type=float, help="Percent dropout.", default = .2, required=False)
    parser.add_argument("--epochs", type=int, help="Number of epochs.", default = 25, required=False)
    parser.add_argument("--batch_size", type=int, help="Batch size.", default = 1000, required=False)
    parser.add_argument("--verbose", type=int, help="keras verbosity", default = 1, required=False)
    parser.add_argument("--mylossfunc", help="Loss function.", default = 'binary_crossentropy',
        choices=['binary_crossentropy', 'categorical_crossentropy', 'sparse_categorical_crossentropy'], required=False)
    #
    # parameters for cyclical learning Rate, see https://github.com/bckenstler/CLR for details
    parser.add_argument("--numCycles", help="Number of cyclical learning rate cycles.", 
        type=float, default = 2.5, required=False)
    parser.add_argument("--base_lr", help="Learning rate floor value for cyclical learning rate model.", 
        type=float, default = 2e-3, required=False)
    parser.add_argument("--max_lr", help="Learning rate ceiling value for cyclical learning rate model.", 
        type=float,default = .2, required=False)
    parser.add_argument("--base_m", help="Momentum floor value for cyclical learning rate model.", 
        type=float, default = .875, required=False)
    parser.add_argument("--max_m", help="Momentum ceiling value for cyclical learning rate model.", 
        type=float,default = .99, required=False)
    parser.add_argument("--cyclical_momentum",help="Whether to cycle the momemtum for OCP.",
        action='store_true')
    # parse arguments
    args = parser.parse_args()
    # args = parser.parse_args(['Mo2015_EXCpos_Ctx_fold1', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_trainPos.fa', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_trainNeg10x.fa', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validPos.fa', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validNeg10x.fa',
    #     '--cyclical_momentum'])
    # args = parser.parse_args(['Mo2015_EXCpos_Ctx_fold1', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_trainPos.fa', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_trainNeg10x.fa', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validPos.fa', 
    #     'FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validNeg10x.fa',
    #     '--cyclical_momentum'])
    main(args)








