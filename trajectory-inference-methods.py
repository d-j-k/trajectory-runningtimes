###################################
## Trajecctory Inference Methods ##
###################################

import time


############
## Scanpy ##
############
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy.api as sc


def runScanpyAGA(counts):
    adata = sc.AnnData(counts.transpose())
    adata.var_names = list("gene" + str(x) for x in range(counts.shape[0]))
    adata.row_names = list("cell" + str(x) for x in range(counts.shape[1]))
    start_time = time.time()
    sc.tl.tsne(adata, n_pcs=30)
    sc.tl.aga(adata,n_neighbors=4, resolution=1.2)
    end_time = time.time()
    dt = round(end_time - start_time, 2)
    return dt


def runScanpyDPT(counts):
    adata = sc.AnnData(counts.transpose())
    adata.var_names = list("gene" + str(x) for x in range(counts.shape[0]))
    adata.row_names = list("cell" + str(x) for x in range(counts.shape[1]))
    start_time = time.time()
    adata.uns['iroot'] = 0
    sc.tl.dpt(adata)
    end_time = time.time()
    dt = round(end_time - start_time, 2)
    return dt


##############
## Wishbone ##
##############

import wishbone
import pandas as pd

def runWishbone(counts):
    counts = pd.DataFrame(counts)
    scdata = wishbone.wb.SCData(counts.transpose(), data_type='sc-seq')
    start_time = time.time()
    scdata.run_pca()
    scdata.run_diffusion_map()
    wb = wishbone.wb.Wishbone(scdata)
    wb.run_wishbone(start_cell=scdata.data.index[0], num_waypoints=min(counts.shape[1], 250))
    end_time = time.time()
    dt = round(end_time - start_time, 2)
    return dt


#############
## GPfates ##
#############

import GPy
import pandas as pd
import numpy as np
from GPfates import GPfates



def runGPfates(counts):
    counts = pd.DataFrame(counts)
    logcounts = np.log10(counts + 1)

    start_time = time.time()
    bgplvm = GPy.models.BayesianGPLVM(logcounts.T.as_matrix(), input_dim=2)
    bgplvm.optimize(messages=False, max_iters=100)
    result = bgplvm.latent_space.mean.tolist()
    metadata = pd.DataFrame(result)
    metadata.columns = metadata.dtypes.index.map(str)

    m = GPfates.GPfates(metadata, logcounts)
    m.dimensionality_reduction()
    m.infer_pseudotime(s_columns=['0', '1'])
    # model 2 trajectories using BGPLVM components 0, 1
    m.model_fates(t='pseudotime', X=['0', '1'])
    end_time = time.time()
    dt = round(end_time - start_time, 2)
    return dt
