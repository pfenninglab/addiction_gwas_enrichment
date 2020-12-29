import os, sys, gc, math, argparse
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from sklearn import metrics
from cnn_utils import *
from Bio import SeqIO

df = pd.read_feather('rdas/SingleTask_IDR_peakPrediction_CNN_performance.feather')

# fasta files
fastaEffect = 'FASTA/addiction_snps_allele_effect_501.fa'
fastaNonEff = 'FASTA/addiction_snps_allele_nonEffect_501.fa'

# one-hot encode the sequences and labels
x_effect, lab_effect = encode_sequence3(fastaEffect)
x_nonEff, lab_nonEff = encode_sequence3(fastaNonEff)

########################################
# create dictionary of model predictions
dict_effect = {'label' : lab_effect}; dict_nonEff = {'label' : lab_nonEff}

# predict on effect alleles
for model_name in df['model']:
	model = load_model(model_name, compile=False)
	# for effect allele
	effect_pred_score = model.predict(x_effect, verbose = 1, batch_size = 2000)
	dict_effect.update( {model_name : effect_pred_score.flatten()} )
	# for non effect allele
	nonEff_pred_score = model.predict(x_nonEff, verbose = 1, batch_size = 2000)
	dict_nonEff.update( {model_name : nonEff_pred_score.flatten()} )

# average the SNP scores between forward and reverse complement orientation DNA
df_effect = pd.DataFrame(dict_effect).groupby('label').agg([np.mean])
df_nonEff = pd.DataFrame(dict_nonEff).groupby('label').agg([np.mean])

df_effect.to_csv('rdas/SingleTask_IDR_addiction_snps_allele_effect_predictions.tsv', sep = '\t')
df_nonEff.to_csv('rdas/SingleTask_IDR_addiction_snps_allele_nonEffect_predictions.tsv', sep = '\t')
