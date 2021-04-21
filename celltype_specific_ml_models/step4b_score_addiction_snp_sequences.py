import os, sys, gc, math, argparse
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from sklearn import metrics
from cnn_utils import *
from train_singleTask_CNN_classifier_OCP import *
from Bio import SeqIO

df = pd.read_csv('tables/SingleTask_IDR_bestModels_CNN_performance.tsv', sep = '\t')

# fasta files
fastaEffect = 'FASTA/addiction_snps_allele_effect_501.fa'
fastaNonEff = 'FASTA/addiction_snps_allele_nonEffect_501.fa'

# one-hot encode the sequences and labels
x_effect, lab_effect = encode_sequence3(fastaEffect, size = 501)
x_nonEff, lab_nonEff = encode_sequence3(fastaNonEff, size = 501)

########################################
# create dictionary of model predictions
dict_effect = {'label' : lab_effect}; dict_nonEff = {'label' : lab_nonEff}

# average the SNP scores between forward and reverse complement orientation DNA
df_effect = pd.DataFrame({'label' : []})
df_nonEff = pd.DataFrame({'label' : []})

# # predict on effect alleles
# for model_name in df['model']:
# 	model = load_model(model_name, compile=False)
# 	# for effect allele
# 	effect_pred_score = model.predict(x_effect, verbose = 1, batch_size = 2000)
# 	dict_effect.update( {model_name : effect_pred_score.flatten()} )
# 	# for non effect allele
# 	nonEff_pred_score = model.predict(x_nonEff, verbose = 1, batch_size = 2000)
# 	dict_nonEff.update( {model_name : nonEff_pred_score.flatten()} )

# # average the SNP scores between forward and reverse complement orientation DNA
# df_effect = pd.DataFrame(dict_effect).groupby('label').agg([np.mean])
# df_nonEff = pd.DataFrame(dict_nonEff).groupby('label').agg([np.mean])

# predict on effect alleles
for model_name in df['model']:
	print(f'Scoring {model_name}.')
	# for effect allele
	effect_pred_score = predict_sequences(model_name, x_effect, lab_effect)
	effect_pred_score['model'] = model_name
	df_effect = df_effect.append(effect_pred_score)
	# for non effect allele
	nonEff_pred_score = predict_sequences(model_name, x_nonEff, lab_nonEff)
	nonEff_pred_score['model'] = model_name
	df_nonEff = df_nonEff.append(nonEff_pred_score)

df_effect['label'] = df_effect.index.str.split(pat = '_').str[1]
df_nonEff['label'] = df_nonEff.index.str.split(pat = '_').str[1]

df_effect.reset_index().to_csv('tables/SingleTask_IDR_bestModels_addiction_snps_allele_effect_predictions.tsv.gz', sep = '\t',compression='gzip')
df_nonEff.reset_index().to_csv('tables/SingleTask_IDR_bestModels_addiction_snps_allele_nonEffect_predictions.tsv.gz', sep = '\t',compression='gzip')
