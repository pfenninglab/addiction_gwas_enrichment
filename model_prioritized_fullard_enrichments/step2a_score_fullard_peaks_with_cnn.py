import os, sys, gc, math, argparse, glob
import pandas as pd
import numpy as np

# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../celltype_specific_ml_models')
from cnn_utils import *
from train_singleTask_CNN_classifier_OCP import *

## get all models to score sequences
models =  pd.read_csv('../celltype_specific_ml_models/tables/SingleTask_IDR_bestModels_CNN_performance.tsv', sep = '\t')

# fasta files
for fasta in sorted(glob.glob('FASTA/*501.fasta')):
	print(f'Extracting on {fasta} file.')
	fileName = f'tables/cnn_predictions/{os.path.basename(fasta).split("_501.fasta")[0]}_v3.tsv.gz'
	if not os.path.exists(fileName):
		# one-hot encode the sequences and labels
		x, lab = encode_sequence3(fasta, size = 501)
		# create dictionary of model predictions
		dt = pd.DataFrame({'label' : []})
		# for each set of fasta sequences, score them with the set of SingleTask models
		for model_name in models['model']:
			print(f'Predicting with {os.path.basename(model_name)} model.')
			# model_name2 = '../celltype_specific_ml_models/' + model_name
			# model = load_model(model_name2, compile=False)
			# y_pred_score = model.predict(x, verbose = 0, batch_size = 2000)
			# dt.update( {model_name : y_pred_score.flatten()} )
			y_pred_score = predict_sequences(model_name, x, lab)
			y_pred_score['model'] = model_name
			dt = dt.append(y_pred_score)
		# average the scores between forward and reverse complement orientation DNA
		# df = pd.DataFrame(dt).groupby('label').agg([np.mean])
		# write the filename out to file
		dt['label'] = dt.index.str.split(pat = '_').str[1]
		dt.reset_index().to_csv(fileName, sep = '\t', compression='gzip')



