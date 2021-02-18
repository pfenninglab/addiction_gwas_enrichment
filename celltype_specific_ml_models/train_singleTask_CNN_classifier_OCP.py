from tensorflow.keras.layers import Dense, Dropout, Activation, Embedding, Conv2D, MaxPooling2D, Flatten
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import SGD, Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras.models import load_model
import tensorflow.keras.backend as K
import matplotlib; matplotlib.use('agg')
import os, sys, gc, math, argparse
import matplotlib.pyplot as plt
import tensorflow as tf
import pandas as pd
import numpy as np

from sklearn import metrics
from callback_ocp import *
from cnn_utils import *
from Bio import SeqIO

# gc.collect()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# make sure running on GPU #
tf.config.experimental.list_physical_devices('GPU')
# tf.get_logger().setLevel('WARNING')
np.seterr(divide = 'ignore')


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
    clr = CyclicLR(base_lr=args.base_lr,
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




def predict_sequences(model_name, x, ids):
    # Creating an empty Dataframe with column names only
    # predict labels
    # combineStrands will average the logit of the for and rev DNA strands
    model = load_model(model_name, compile=False)
    y_pred_score = model.predict(x, verbose = args.verbose).flatten()
    y_pred_score = np.nan_to_num(y_pred_score)
    df = pd.DataFrame({'id': ids, 'y_pred_logit': np.log10(y_pred_score) - np.log10(1-y_pred_score)})
    df = df.groupby(['id']).mean()
    df['y_pred_score'] = 1 / (1 + np.exp(-df['y_pred_logit']))
    df['y_pred_class'] = df['y_pred_score'] > 0.5
    return df




def predict_sequences2(model_name, x, y, ids):
    # Creating an empty Dataframe with column names only
    # predict labels
    # combineStrands will average the logit of the for and rev DNA strands
    model = load_model(model_name, compile=False)
    y_pred_score = model.predict(x, verbose = args.verbose).flatten()
    y_pred_score = np.nan_to_num(y_pred_score)
    df = pd.DataFrame({'id': ids, 'y' : y, 'y_pred_logit': np.log10(y_pred_score) - np.log10(1-y_pred_score)})
    df = df.groupby(['id']).mean()
    df['y_pred_score'] = 1 / (1 + np.exp(-df['y_pred_logit']))
    df['y_pred_class'] = df['y_pred_score'] > 0.5
    return df



def evaluate_sequences(model_name, x, y, ids, args):
    # Creating an empty Dataframe with column names only
    df = pd.DataFrame(vars(args), index=[0])
    tmp = predict_sequences2(model_name, x, y, ids)
    # compute prediction statistics 
    accuracy = metrics.balanced_accuracy_score(tmp['y'], tmp['y_pred_class'])
    f1_score = metrics.f1_score(tmp['y'], tmp['y_pred_class'], average = 'weighted')
    fhalf_score = metrics.fbeta_score(tmp['y'], tmp['y_pred_class'], beta = 0.5, average = 'weighted')
    roc_auc = metrics.roc_auc_score(tmp['y'], tmp['y_pred_score'])
    precision, recall, thresholds = metrics.precision_recall_curve(tmp['y'], tmp['y_pred_score'])
    prc_auc = metrics.auc(recall, precision) # x, y
    print(f'Accuracy: {accuracy}. ')
    print(f'f1_score: {f1_score}.')
    print(f'roc_auc: {roc_auc}.')
    print(f'prc_auc: {prc_auc}.')
    # add this row to dataframe
    df = pd.concat([df.reset_index(drop=True),
        pd.DataFrame({'model': model_name, 'accuracy': accuracy, 'auROC': roc_auc, 'auPRC': prc_auc,
            'f1_score': f1_score, 'fhalf_score': fhalf_score}, index=[0])], axis = 1)    
    return df



def plot_training_performance(args, hist, clr ):
    """plot function of training history and learning rate ramps
    Args:
        args (argparse):
        hist (keras training history):
        clr (argparse):
    """
    # plot performance by iteration
    if not os.path.exists(f'{args.out_dir}/plots/{label}'):
        os.makedirs(f'{args.out_dir}/plots/{label}')
    fig, axs = plt.subplots(2, 2, figsize=(8, 6))
    fig.suptitle(f"Training performance {label}")
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
    print(f'Running cyclical learning rate for {label}.')
    print(f'One-cycle policy training for {args.epochs} epochs.')
    print(f'Model is: {args.model_name}')
    print(f'Batch size: {args.batch_size}.')
    print(f'Learning rates range: {args.base_lr} - {args.max_lr}.')
    print(f'Momentum rates range: {args.base_m} - {args.max_m}.')
    print(f'Dropout: {args.dropout}.')
    # call main functions
    if args.mode == 'train':
        print('In training mode.')
        if os.path.exists(args.model_name) and not args.force:
            print(f'The model exists w/o permission to overwrite. Use --force to overwrite.')
            return
        (x_train, y_train, ids_train) = encode_sequence(args.train_fasta_pos, args.train_fasta_neg, size = args.seq_length, shuffleOff = False)
        (x_valid, y_valid, ids_valid) = encode_sequence(args.valid_fasta_pos, args.valid_fasta_neg, size = args.seq_length, shuffleOff = False)
        #
        model, clr, hist = train_model_clr(x_train, y_train, x_valid, y_valid, args)
        plot_training_performance(args, hist, clr )
        # make sure to create folder
        if not os.path.exists(f'{args.out_dir}/models/{args.prefix}'):
            os.makedirs(f'{args.out_dir}/models/{args.prefix}') 
        model.save(args.model_name)
        #
    elif args.mode == 'evaluate':
        print('In evaluation mode.')
        if not os.path.exists(args.model_name):
            print('No model found with specified training parameters. Please train model first.')
            return
        (x_valid, y_valid, ids_valid) = encode_sequence(args.valid_fasta_pos, args.valid_fasta_neg, size = args.seq_length, shuffleOff = True)
        df = evaluate_sequences(args.model_name, x_valid, y_valid, ids_valid, args)
        model_test_performance = f'{args.out_dir}/predictions/cnn_v2/{args.prefix}/{label}.performance.feather'
        # save model performances to feather object
        if not os.path.exists(f'{args.out_dir}/predictions/cnn_v2/{args.prefix}'):
            os.makedirs(f'{args.out_dir}/predictions/cnn_v2/{args.prefix}')
        df.to_feather(model_test_performance)
        print(f'Model performance written to {model_test_performance}')
        #
    elif args.mode == 'predict':
        print('In prediction mode.')
        if not os.path.exists(args.model_name):
            print('No model found with specified training parameters. Please train model first.')
            return
        (x, ids) = encode_sequence3(args.predict_fasta, size = args.seq_length)
        df = predict_sequences(args.model_name, x, ids)
        model_predictions = f'{args.out_dir}/predictions/cnn_v2/{args.prefix}/{label}.predictions.txt'
        # save model performances to feather object
        if not os.path.exists(f'{args.out_dir}/predictions/cnn_v2/{args.prefix}'):
            os.makedirs(f'{args.out_dir}/predictions/cnn_v2/{args.prefix}') 
        np.savetxt(model_predictions, df, axis = 1)
        print(f'Prediction written to {model_predictions}')
    return


if __name__ == '__main__':  
    #### set cnn parameters:
    parser = argparse.ArgumentParser(description='Parse CNN parameters.')
    parser.add_argument("--mode", type=str, help="Mode to perform. Train needs all fasta. Evaluate needs validation fasta. Predict only fasta passed predict_fasta.", 
        choices=['train', 'evaluate', 'predict'], default = 'train', required=False)
    #
    parser.add_argument("--prefix", type=str, help="prefix of model.")
    parser.add_argument("--model_name", type=str, help="complete model name")
    parser.add_argument("--predict_fasta", type=str, help="fasta sequence file for predictions.")
    parser.add_argument("--train_fasta_pos", type=str, help="training fasta sequence file of positives.")
    parser.add_argument("--train_fasta_neg", type=str, help="training fasta sequence file of negatives.")
    parser.add_argument("--valid_fasta_pos", type=str, help="validation fasta sequence file of positives.")
    parser.add_argument("--valid_fasta_neg", type=str, help="validation fasta sequence file of negatives.")
    #
    #### set cnn parameters:
    parser.add_argument("--seq_length", type=int, help="DNA sequence length.", default = 501, required=False)
    parser.add_argument("--conv_width", type=int, help="Convolution width.", default = 11, required=False)
    parser.add_argument("--conv_filters", type=int, help="Convolution filter width.", default = 200, required=False)
    parser.add_argument("--max_pool_size", type=int, help="Max pool size.", default = 26, required=False)
    parser.add_argument("--max_pool_stride", type=int, help="Max pool stride.", default = 26, required=False)
    parser.add_argument("--dense_filters", type=int, help="Number of dense filters.", default = 300, required=False)
    parser.add_argument("--l2_reg", type=float, help="L2 regularization rate.", default = 1e-10, required=False)
    parser.add_argument("--dropout", type=float, help="Percent dropout.", default = .2, required=False)
    parser.add_argument("--verbose", type=int, help="keras verbosity", default = 1, required=False)
    parser.add_argument("--mylossfunc", help="Loss function.", default = 'binary_crossentropy',
        choices=['binary_crossentropy', 'categorical_crossentropy', 'sparse_categorical_crossentropy'], required=False)
    #
    #### parameters for cyclical learning Rate, see https://github.com/bckenstler/CLR for details
    parser.add_argument("--batch_size", type=int, help="Batch size.", default = 1000, required=False)
    parser.add_argument("--epochs", type=int, help="Number of epochs.", default = 23, required=False)
    parser.add_argument("--numCycles", help="Number of cyclical learning rate cycles.", 
        type=float, default = 2.35, required=False)
    parser.add_argument("--base_lr", help="Learning rate floor value for cyclical learning rate model.", 
        type=float, default = 1e-2, required=False)
    parser.add_argument("--max_lr", help="Learning rate ceiling value for cyclical learning rate model.", 
        type=float,default = 1e-1, required=False)
    parser.add_argument("--base_m", help="Momentum floor value for cyclical learning rate model.", 
        type=float, default = .85, required=False)
    parser.add_argument("--max_m", help="Momentum ceiling value for cyclical learning rate model.", 
        type=float, default = .99, required=False)
    parser.add_argument("--cyclical_momentum", help="Whether to cycle the momemtum for OCP.", default = False, 
        action='store_true')
    parser.add_argument("--force", help="Whether to overwrite previously trained model.",
        action='store_true')
    parser.add_argument("--out_dir", type=str, default = '.', help="path to ouputput directory, default is pwd")

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

    args = parser.parse_args()
    if args.model_name is None:
        label = f'{args.prefix}_OCP_NB{args.batch_size}_NE{args.epochs}_BR{args.base_lr}_MR{args.max_lr}_BM{args.base_m}_MM{args.max_m}_DO{args.dropout}'
        args.model_name = f"{args.out_dir}/models/{args.prefix}/{label}.h5"
    else:
        # parse the model name to get the OCP parameters
        label = os.path.splitext(os.path.basename(args.model_name))[0]
        args.prefix     = label.split('_OCP')[0]
        args.batch_size = int(label.split('_NB')[1].split('_')[0])
        args.epochs     = int(label.split('_NE')[1].split('_')[0])
        args.base_lr    = float(label.split('_BR')[1].split('_')[0])
        args.max_lr     = float(label.split('_MR')[1].split('_')[0])
        args.base_m     = float(label.split('_BM')[1].split('_')[0])
        args.max_m      = float(label.split('_MM')[1].split('_')[0])
        args.dropout    = float(label.split('_DO')[1].split('_')[0])

    main(args)



