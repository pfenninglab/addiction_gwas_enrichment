import os, sys, gc, math, argparse, glob
import pandas as pd
import numpy as np

# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../celltype_specific_ml_models')
from cnn_utils import *

## get all models to score sequences
models = pd.read_feather('../celltype_specific_ml_models/rdas/SingleTask_IDR_peakPrediction_CNN_performance.feather')

# fasta files
for fasta in sorted(glob.glob('FASTA/*501.fasta')):
	print(f'Extracting on {fasta} file.')
	# one-hot encode the sequences and labels
	x, lab = encode_sequence3(fasta)
	# create dictionary of model predictions
	dt = {'label' : lab}
	# for each set of fasta sequences, score them with the set of SingleTask models
	for model_name in models['model']:
		model_name2 = '../celltype_specific_ml_models/' + model_name
		print(f'Predicting with {os.path.basename(model_name)} model.')
		model = load_model(model_name2, compile=False)
		y_pred_score = model.predict(x, verbose = 0, batch_size = 2000)
		dt.update( {model_name : y_pred_score.flatten()} )
	# average the scores between forward and reverse complement orientation DNA
	df = pd.DataFrame(dt).groupby('label').agg([np.mean])
	# write the filename out to file
	fileName = f'tables/cnn_predictions/{os.path.basename(fasta).split("_501.fasta")[0]}.tsv.gz'
	df.to_csv(fileName, sep = '\t', compression='gzip')



