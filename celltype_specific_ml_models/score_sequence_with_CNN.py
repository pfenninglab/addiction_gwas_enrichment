########################
# make sure running on GPU
import pandas as pd
from Bio import SeqIO
import numpy as np
import glob
import os
import argparse

# saveFile = 'rdas/' + os.path.basename(npzFile).replace('167_1000.npz', 'cnn_scored.feather',)

def onehot_seq_ACGT(seq):
	letter_to_index = {'A':0, 'a':0,
	                 'C':1, 'c':1,
	                 'G':2, 'g':2,
	                 'T':3, 't':3}
	to_return = np.zeros((len(seq),4))
	for idx,letter in enumerate(seq):
		if letter not in ['N','n']:
			to_return[idx,letter_to_index[letter]] = 1
	return to_return

def onehot_seq_AGCT(seq):
	letter_to_index = {'A':0, 'a':0,
	                 'G':1, 'g':1,
	                 'C':2, 'c':2,
	                 'T':3, 't':3}
	to_return = np.zeros((len(seq),4))
	for idx,letter in enumerate(seq):
		if letter not in ['N','n']:
			to_return[idx,letter_to_index[letter]] = 1
	return to_return

def logit(x):
	to_return = np.ma.log10(np.divide(x,1-x, out=np.ones_like(x), where=(1-x)!=0))
	return to_return

def fasta2oneHot(fastaFile167, fastaFile1000, npzFile, force = False):
# must have both files b/c models trained w/ different sequence lengths
	if not os.path.exists() or force:
		print('Encoding {}.'.format(fastaFile))
		# one hot with ACGT
		xlab =  [seq.id for seq in SeqIO.parse(fastaFile167, "fasta") ]
		xtest167= np.array([onehot_seq_ACGT(seq) for seq in SeqIO.parse(fastaFile167, "fasta")])
		xtest167 = np.expand_dims(xtest167, axis=3)
		# one hot with AGCT
		print('Encoding {}.'.format(fastaFile1000))
		xtest1000= np.array([onehot_seq_AGCT(seq) for seq in SeqIO.parse(fastaFile1000, "fasta")])
		print('Saving one-hot encoding to {}.'.format(npzFile))
		np.savez_compressed(npzFile, xtest167 = xtest167, xtest1000 = xtest1000, xlab = xlab)
	else:
		print('NPZ file, {}, exists. Set FORCE to overwrite.'.format(fastaFile))

def scoreFastaWithCNN(npzFile, saveFile):
	import tensorflow as tf
	from tensorflow.keras.models import load_model
	tf.config.experimental.list_physical_devices('GPU')
	################################################
	# score addiction SNPs with D1/D2 & PV CNN models
	msn_models = ['../celltype_specific_ml_models/CNN/D1D2vsBulk/D1vsD2_cnn_model_lr0.0006483217546826382_freeze0_stop3.h5',
				'../celltype_specific_ml_models/CNN/D1vsBulk/D1vsBulk_cnn_model_lr0.0008183410786692336_freeze0_stop3.h5',
				'../celltype_specific_ml_models/CNN/D1vsD2/D1vsD2_cnn_model_lr0.0014952502941741953_freeze4_stop14.h5',
				'../celltype_specific_ml_models/CNN/D2vsBulk/D2vsBulk_cnn_model_lr8.909599158518026e-05_freeze0_stop4.h5']
	msn_model_type = ['D1D2vsBulk','D1vsBulk','D1vsD2','D2vsBulk']
	#
	# remember onehot for PV models is AGCT not alpha numeric
	pv_models = ['/projects/pfenninggroup/machineLearningForComputationalBiology/eramamur_stuff/ml_mouse_pv_vip_exc_atac/differential_peaks_models/pv_vs_exc_differential_model_2.hdf5',
	'/projects/pfenninggroup/machineLearningForComputationalBiology/eramamur_stuff/ml_mouse_pv_vip_exc_atac/differential_peaks_models/pv_vs_vip_differential_model_2.hdf5',
	'/projects/pfenninggroup/machineLearningForComputationalBiology/eramamur_stuff/ml_mouse_pv_vip_exc_atac/differential_peaks_models/pv_vs_pvneg_differential_model_2.hdf5']
	pv_model_type = ['PVvsEXC','PVvsVIP','PVvsPVNEG']
	#
	################################################
	# score fullard IDR peak sequences with D1/D2 & PV CNN models
	print('Processing {} ...'.format(npzFile.replace('_167_1000.npz', '',)))
	tmp = np.load(npzFile)
	xlab = tmp['xlab']
	xtest167 = tmp['xtest167']
	xtest1000 = tmp['xtest1000']
	#
	# initiate the score matrix
	scoreMat = np.zeros([len(xlab), 7])
	#
	# score with the MSN models w/ 167bp
	for msnInd in np.arange(4):
		model_file = msn_models[msnInd]
		model_type = msn_model_type[msnInd]
		print('\tScoring sequences with {} ...'.format(model_type))
		model = load_model(model_file, compile = False)
		# invert the sigmoid probabilities w/ logit function
		score = logit(model.predict(xtest167, batch_size = 2000, verbose=1))
		scoreMat[:,msnInd] = np.squeeze(score)
	#
	# score with the PV models w/ 1000bp
	for pvInd in np.arange(3):
		model_file = pv_models[pvInd]
		model_type = pv_model_type[pvInd]
		print('\tScoring sequences with {} ...'.format(model_type))
		model = load_model(model_file, compile = False)
		# invert the sigmoid probabilities w/ logit function
		score = logit(model.predict(xtest1000, batch_size = 200, verbose=1))
		scoreMat[:,pvInd+4] = np.squeeze(score)
	#
	print('Saving scored sequences to {}.'.format(saveFile))
	df = pd.DataFrame(data=scoreMat, index=xlab, columns=msn_model_type + pv_model_type)  
	df = df.reset_index()
	df.to_feather(saveFile) # save the file

def main():
	# read in the command line arguments
	usage = "usage: python score_sequence_with_CNN.py --fastaFile167 file167.fasta --fastaFile1000 file1000.fasta --saveFile file.feather"
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument("--fastaFile167", metavar= 'file167.fa',
		action="store", dest="fastaFile167", required=True, 
		help="input fasta file with length 167bp")
	parser.add_argument("--fastaFile1000", metavar= 'file1000.fa',
		action="store", dest="fastaFile1000",  required=True,
		help="input fasta file with length 1000bp")
	parser.add_argument("--npzFile", metavar= 'file.npz',
		action="store", dest="npzFile", required=True,
		help="numpy zipped file to save one-hot encoding of fasta sequences")
	parser.add_argument("--saveFile",  required=True,
		action="store", dest="saveFile", 
		help="write output to saveFile feather file")
	parser.add_argument("--force", 
		action="store_true", dest="force", 
		help="overwrite npzFile and saveFile")

	args = parser.parse_args(['--fastaFile167', 'file167.fasta', '--fastaFile1000', 'file1000.fasta', 
		'--saveFile', 'file.feather'])

	# testing arguments
	# print(args.fastaFile167)
	# print(args.fastaFile1000)
	# print(args.saveFile)	
	# print(args.npzFile)

	# make npzFile of fasta sequences
	fasta2oneHot(args.fastaFile167, args.fastaFile1000, args.npzFile)
	# score the files
	scoreFastaWithCNN(args.npzFile, args.saveFile)

if __name__ == "__main__":
    main()


