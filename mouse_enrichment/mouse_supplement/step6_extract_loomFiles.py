import loompy
import scanpy as sc
import anndata


lm = anndata.read_loom('loom/l2_neurons_cortex1.agg.loom', 
	sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', 
	obsm_names=None, var_names='Gene', varm_names=None, dtype='float32')


