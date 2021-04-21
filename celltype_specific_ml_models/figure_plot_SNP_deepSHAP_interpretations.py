from tensorflow.keras.models import load_model
import tensorflow.keras.backend as K
import os, sys, gc, math, argparse
import tensorflow as tf
import pandas as pd
import numpy as np
import shap
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from sklearn import metrics
from callback_ocp import *
from cnn_utils import *
from deepSHAP_importance_scores import *
from modisco.visualization import viz_sequence
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# make sure running on GPU #
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
tf.config.experimental.list_physical_devices('GPU')
np.seterr(divide = 'ignore')


def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1


def normalize_scores(impscores, hyp_impscores, onehot_data):
  #normalize the hyp scores such that, at each position, hypothetical importance
  # scores that have the same sign as the original importance score all sum
  # up to the original importance score value. The rationale is that if
  # multiple different bases at a position could produce a similar score,
  # the specific identity of each individual base is less important.
  #Empirically, hypothetical scores like these appear to work better for
  # motif discovery. Using normalized importance scores derived by taking
  # the elementwise product of the normalized hypothetical scores and
  # the one-hot encoding also seems to reduce noise.
  normed_hyp_impscores = []
  normed_impscores = []
  for i in range(len(impscores)):
      imp_score_each_pos = np.sum(impscores[i],axis=-1)
      imp_score_sign_each_pos = np.sign(imp_score_each_pos)
      hyp_scores_same_sign_mask = (np.sign(hyp_impscores[i])
                                   *imp_score_sign_each_pos[:,None] > 0)
      hyp_scores_same_sign_imp_scores_sum = np.sum(
          hyp_impscores[i]*hyp_scores_same_sign_mask,axis=-1)
      norm_ratio = imp_score_each_pos/hyp_scores_same_sign_imp_scores_sum
      norm_hyp = hyp_impscores[i]*norm_ratio[:,None]
      normed_hyp_impscores.append(norm_hyp)
      normed_impscores.append(norm_hyp*onehot_data[i])
  return normed_impscores, normed_hyp_impscores


###########################################
#read in the fasta files and one-hot encode
def get_scores(fasta_file, imp_file, hyp_file):
    print("Getting sequences from:",fasta_file)
    fasta_seqs = [record.seq for record in SeqIO.parse(fasta_file, "fasta") ]
    fasta_ids = [record.id for record in SeqIO.parse(fasta_file, "fasta") ]
    #filter out any sequences that contain 'N's
    onehot_data = [np.array(one_hot_encode_along_channel_axis(x)) for x in fasta_seqs if ('N' not in x)]
    print("Num onehot sequences:",len(onehot_data))
    #read in the importance scores and hypothetical importance scores
    #filter out any sequences that contain 'N's
    hyp_impscores = [w[0] for w in zip([
        np.array( [[float(z) for z in y.split(",")]
                    for y in x.rstrip().split("\t")[2].split(";")])
        for x in open(hyp_file)
    ],fasta_seqs) if 'N' not in w[1]]
    impscores = [w[0] for w in zip([
        np.array( [[float(z) for z in y.split(",")]
                    for y in x.rstrip().split("\t")[2].split(";")])
        for x in open(imp_file)
    ],fasta_seqs) if 'N' not in w[1]]
    assert (np.max([np.max(np.abs(z*y - x)) for x,y,z in zip(impscores, onehot_data, hyp_impscores)]))==0
    # normalize important scores
    normed_impscores, normed_hyp_impscores = normalize_scores(
      impscores=impscores, hyp_impscores=hyp_impscores, onehot_data=onehot_data)
    return fasta_seqs, fasta_ids, normed_impscores, normed_hyp_impscores

###########################################
#read in the fasta files and one-hot encode

for fold in range(1,6):
    fasta_file = f"shap/top_addiction_snps_effect_allele.Pfenning_D1pos_Nac_fold{fold}.pos110.fasta"
    imp_file = f"shap/top_addiction_snps_effect_allele.Pfenning_D1pos_Nac_fold{fold}.predict.imp_SHAP_scores.txt"
    hyp_file = f"shap/top_addiction_snps_effect_allele.Pfenning_D1pos_Nac_fold{fold}.predict.hyp_SHAP_scores.txt"
    fasta_seqs, fasta_ids, normed_impscores, normed_hyp_impscores = get_scores(fasta_file, imp_file, hyp_file)
    #
    fasta_file2 = f"shap/top_addiction_snps_nonEffect_allele.Pfenning_D1pos_Nac_fold{fold}.pos110.fasta"
    imp_file2 = f"shap/top_addiction_snps_nonEffect_allele.Pfenning_D1pos_Nac_fold{fold}.predict.imp_SHAP_scores.txt"
    hyp_file2 = f"shap/top_addiction_snps_nonEffect_allele.Pfenning_D1pos_Nac_fold{fold}.predict.hyp_SHAP_scores.txt"
    fasta_seqs2, fasta_ids2, normed_impscores2, normed_hyp_impscores2 = get_scores(fasta_file2, imp_file2, hyp_file2)
    #
    index = [i for i in range(len(fasta_ids)) if fasta_ids[i] in ['rs7604640']]
    index2 = [i for i in range(len(fasta_ids2)) if fasta_ids2[i] in ['rs7604640']]
    #
    pp = PdfPages(f'plots/addiction_snp_D1pos_Nac_fold{fold}_rs7604640_delta_snp.pdf')
    for idx,idx2 in zip(index, index2):
        print("Effect Allele: ",fasta_ids[idx],"total imp",np.sum(impscores[idx]))
        viz_sequence.plot_weights(normed_impscores[idx][200:300,:], subticks_frequency=10)
        pp.savefig()
        print("Non-effect Allele: ",fasta_ids[idx2],"total imp",np.sum(impscores2[idx2]))
        viz_sequence.plot_weights(normed_impscores2[idx2][200:300,:], subticks_frequency=10)
        pp.savefig()
        delta = normed_impscores[idx][200:300,:] - normed_impscores2[idx2][200:300,:]
        viz_sequence.plot_weights(delta, subticks_frequency=10)
        pp.savefig()
    #
    pp.close()


