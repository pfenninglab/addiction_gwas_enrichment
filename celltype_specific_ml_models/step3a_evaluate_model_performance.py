from sklearn import metrics
import numpy as np
import pandas as pd
import glob, os, sys


class Namespace:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)

#####################################################
# Creating an empty Dataframe with column names only
df = pd.DataFrame()

#################################################
# read in CNN models and scored test sequences
for model_name in sorted(glob.glob('models/*/*.h5')):
	featherFile = model_name.replace('models','predictions').replace('.h5','.performance.feather')
	if os.path.exists(featherFile):
		tmp = pd.read_feather(featherFile)
		tmp['model_species'] = 'Mouse'
		tmp['sample'] = tmp['model'][0].split('/')[2]
		tmp['model_type'] = tmp['model'][0].split('/')[2].split('_fold')[0]
		df = df.append(tmp)

##############################################
# save model performances to feather object
df.reset_index().to_feather('rdas/SingleTask_IDR_peakPrediction_CNN_v2_performance.feather')
df.reset_index().to_csv('tables/SingleTask_IDR_peakPrediction_CNN_v2_performance.tsv', sep = '\t')

