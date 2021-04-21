from sklearn import metrics
import numpy as np
import pandas as pd
import glob, os, sys


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

#####################################################
# Creating an empty Dataframe with column names only
df = pd.DataFrame(columns=[
	'sample','model_class', 'model_species',  'model', 
	'accuracy', 'auROC','auPRC', 'f1_score','fhalf_score'])

#################################################
# read in CNN models and scored test sequences
for model_name in sorted(glob.glob('cnn/singleTask/*/*.h5')):
	featherFile = model_name.replace('cnn/singleTask','predictions/cnn').replace('.h5','_testSet_testPerformance.feather')
	tmp = pd.read_feather(featherFile)
	tmp['model_species'] = 'Mouse'
	tmp['sample'] = tmp['model'][0].split('/')[2]
	tmp['model_type'] = tmp['model'][0].split('/')[2].split('_fold')[0]
	df = df.append(tmp)

##############################################
# save model performances to feather object
df.to_feather('rdas/SingleTask_IDR_peakPrediction_CNN_testSet_performance.feather')
df.reset_index().to_csv('tables/SingleTask_IDR_peakPrediction_CNN_testSet_performance.tsv', sep = '\t')


