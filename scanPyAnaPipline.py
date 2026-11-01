#
# In the script, I used Scanpy (Wolf 2018) implemented with bbknn (Polanski 2019) to analyze the VEN 10X single cell sequencing data.
# But 'bbknn' looks not good enough for our 10X data to remove the batch effects. 


import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import scrublet as scr
import matplotlib.pyplot as plt

sc.settings.n_jobs = 45
sc._settings.ScanpyConfig.n_jobs = 45  # 保留10核给系统
sc._settings.ScanpyConfig.max_memory = 300  # GB，保留20G给系统

sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(6, 6))



#results_file = './write/Insula_all.h5ad'#  file that will store the analysis results
adata1  = sc.read_10x_mtx('./dataMatrix/data_Ctl/VENL20230706/', var_names='gene_symbols', cache=True)
adata2  = sc.read_10x_mtx('./dataMatrix/data_Ctl/VENR20230706/', var_names='gene_symbols', cache=True)
adata3  = sc.read_10x_mtx('./dataMatrix/data_Ctl/AIC1_2', var_names='gene_symbols', cache=True)
adata4  = sc.read_10x_mtx('./dataMatrix/data_Ctl/AIC2_1', var_names='gene_symbols', cache=True)
adata5  = sc.read_10x_mtx('./dataMatrix/data_Ctl/AIC2_2', var_names='gene_symbols', cache=True)
adata6  = sc.read_10x_mtx('./dataMatrix/data_Aut/VENR1', var_names='gene_symbols', cache=True)
adata7  = sc.read_10x_mtx('./dataMatrix/data_Aut/VENR2', var_names='gene_symbols', cache=True)
adata8  = sc.read_10x_mtx('./dataMatrix/data_Aut/Aut_3', var_names='gene_symbols', cache=True)

adata = adata1.concatenate(adata2,adata3,adata4,adata5,adata6,adata7,adata8, batch_categories=['batch-1','batch-2', 'batch-3','batch-4','batch-5','batch-6','batch-7','batch-8'])

del(adata1,adata2,adata3,adata4,adata5,adata6,adat7,adata8)

#----------------------add labels----------------------

Sex = adata.obs['batch']
batch2Sex = dict(); 
batch2Sex['batch-1']='Male'
batch2Sex['batch-2']='Male'
batch2Sex['batch-3']='Male'
batch2Sex['batch-4']='Male'
batch2Sex['batch-5']='Male'
batch2Sex['batch-6']='Male'
batch2Sex['batch-7']='Male'
batch2Sex['batch-8']='Male'
adata.obs['Sex'] = Sex.map(lambda x: batch2Sex.get(x))

Sample = adata.obs['batch']
batch2Sample=dict()
batch2Sample['batch-1']='VENL20230706'
batch2Sample['batch-2']='VENR20230706'
batch2Sample['batch-3']='AIC1_1'
batch2Sample['batch-4']='AIC2_1'
batch2Sample['batch-5']='AIC2_2'
batch2Sample['batch-6']='VENR1'
batch2Sample['batch-7']='VENR2'
batch2Sample['batch-8']='Aut_3'
adata.obs['TissueSample'] = Sample.map(lambda x: batch2Sample.get(x))

HemiBrainDict = dict={
    'batch-1':'Left',
    'batch-2':'Right',
    'batch-3':'Right',
    'batch-4':'Right',
    'batch-5':'Left',
    'batch-6':'Right',
    'batch-7':'Right',
    'batch-8':'Left',
}
adata.obs['BrainHemisphere']=adata.obs['batch'].map(HemiBrainDict).astype('category')

AutCtlDict = dict={
    'batch-1':'Ctrl',
    'batch-2':'Ctrl',
    'batch-3':'Ctrl',
    'batch-4':'Ctrl',
    'batch-5':'Ctrl',
    'batch-6':'Autism',
    'batch-7':'Autism',
    'batch-8':'Autism',
}
adata.obs['AutCtl']=adata.obs['batch'].map(AutCtlDict).astype('category')

AnimalDict = dict={
    'batch-1':'AnimalCtl_2',
    'batch-2':'AnimalCtl_2',
    'batch-3':'AnimalCtl_1',
    'batch-4':'AnimalCtl_1',
    'batch-5':'AnimalCtl_3',
    'batch-6':'AnimalAut_1',
    'batch-7':'AnimalAut_1',
    'batch-8':'AnimalAut_1',
}
adata.obs['AnimalIDs']=adata.obs['batch'].map(AnimalDict).astype('category')

SampleID_Dict = dict={
    'batch-1':'Control_3',
    'batch-2':'Control_4',
    'batch-3':'Control_1',
    'batch-4':'Control_2',
    'batch-5':'Control_5',
    'batch-6':'Model_1',
    'batch-7':'Model_1',
    'batch-8':'Model_1',
}
adata.obs['SampleIDs']=adata.obs['batch'].map(SampleID_Dict).astype('category')



adata.write('./write/scRNAseq_autVSctl_rawCounts.h5ad')



# --------------------------------calculate the doublet scores -----------------------------
def doub_score(ad_wk, m_df):
    
    df = pd.DataFrame([])
    #print(pd.unique(ad_wk.obs['sample']))
    counts_matrix = ad_wk.X
    print(counts_matrix.shape)
    
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets() 
    
    scrub.plot_histogram();
    #print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)
        
    df = pd.DataFrame({'cell_name':ad_wk.obs_names,'Doublet_scores':doublet_scores,'predicted_doublets' : predicted_doublets})
    m_df = pd.concat([m_df, df], axis=0)
    print(m_df.shape)
    
    return m_df


scores_df = pd.DataFrame([])
batch_data = adata[adata.obs['batch']=='batch-1',:]
scores_df= doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-2',:]
scores_df = doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-3',:]
scores_df= doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-4',:]
scores_df = doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-5',:]
scores_df= doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-6',:]
scores_df = doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-7',:]
scores_df= doub_score(batch_data, scores_df)
batch_data = adata[adata.obs['batch']=='batch-8',:]
scores_df = doub_score(batch_data, scores_df)

scores_df.set_index('cell_name', inplace=True)
df_admk = scores_df.loc[adata.obs_names,'Doublet_scores']
df_admk.shape
adata.obs['doublet_score']= df_admk
df_admk = scores_df.loc[adata.obs_names,'predicted_doublets']
adata.obs['predicted_doublets']= df_admk

adata.write('./write/scRNAseq_autVSctl_step1.h5ad')
# ------------------Preparing processing----------------------------------
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx
adata

sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# add mito genes...
adata.var['mt'] = adata.var_names.isin(["ND1","ND2","COX1","COX2","COX3","ATP6","ATP8","ND3","ND4","ND5","ND6","ND4L","CYTB"])
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
adata = adata[adata.obs.pct_counts_mt < 30, :]
adata = adata[~adata.obs['predicted_doublets']]

adata.layers['Counts'] = adata.X

adata.write('./write/scRNAseq_autVSctl_step2.h5ad')





sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


sc.pp.highly_variable_genes(adata, min_mean=0.015, max_mean=10, min_disp=0.05)
sc.pl.highly_variable_genes(adata)
print('# HVG: %d'%adata.var['highly_variable'].sum())

adata.raw = adata

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, use_highly_variable = True)
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, maxiter=200, min_dist=1, spread=2)
sc.pl.umap(adata, color=['batch'])

adata.write('./write/scRNAseq_autVSctl_step3.h5ad')

# ------------Integration datasets batches with harmony-----------------
import harmonypy as hm
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, ['batch'])
adata.obsm['X_harmonypca'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_harmonypca')

adata.obsm['X_umapraw'] = adata.obsm['X_umap']
sc.tl.umap(adata,random_state=1)
adata.obsm['X_umapharmony'] = adata.obsm['X_umap']
sc.pl.umap(adata, color=['batch'], ncols=2, frameon=False)
cmap=scr.darken_cmap(plt.cm.Reds, 0.9)
sc.pl.umap(adata, color=['doublet_score'], legend_loc='on data', title='', frameon=False,use_raw=True,color_map=cmap,size=50)

sc.pl.umap(adata,color=['SLC17A7','GAD1','AQP4','C1QA','C1QL1','MBP','MEIS2','FLT1','BGN','CFAP126','DCN','LYVE1'],use_raw=True,legend_loc='right margin', color_map=cmap, frameon=False,size=5)

sc.tl.tsne(adata, n_pcs=10, use_rep='X_harmonypca')
sc.pl.tsne(adata,color='batch',use_raw=True,legend_loc='right margin', title='', frameon=False,size=5)


# -------------clustering after batch integration--------------------
sc.tl.leiden(adata, resolution=0.6)
#sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata
adata.write('./write/scRNAseq_autVSctl_rawluster.h5ad')







