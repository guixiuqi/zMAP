# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:52:42 2019

@author: guixiuqi
E-mail:  guixiuqi@gmail.com
"""

import scipy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import argparse
import os
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.patheffects import withStroke
import matplotlib.gridspec as gridspec
import pandas as pd
import random
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke
import matplotlib.gridspec as gridspec

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines




def pca_analysis(pearsonr_cor):
    from sklearn.decomposition import PCA
    X_new =  pearsonr_cor.T
    pca = PCA(n_components=5)
    X_r = pca.fit(X_new).transform(X_new)
    print('explained variance ratio (first three components): %s'
      % str(pca.explained_variance_ratio_))
    pc_df = pd.DataFrame(X_r,index = X_new.index,columns = ["PC1","PC2","PC3","PC4","PC5"])
    return pca,pc_df

def get_condition_mark_dict(batch_df,col_=""):

    num = len(set(batch_df[col_].values))
    condition_marker_dict = {}
    marker_list = ['o','^','P','s','p','*','H','D','d','h','X','>']

    i=0
    for index in list(set(batch_df[col_].values)):

        condition_marker_dict[index] = marker_list[i]
        i+=1
    return condition_marker_dict


def get_batch_color_dict(batch_df,col_=""):
    
    num = len(set(batch_df[col_].values))
    batch_color_dict = {}
    color_list = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])for i in range(num)]
                  
    i=0
    for index in list(set(batch_df[col_].values)):
        
        batch_color_dict[index] = color_list[i]
        i+=1
    return batch_color_dict




def get_color_list(columns,sample_info_df,batch_color_dict,condition_color_dict):
    condition_color_list = []
    batch_color_list = []
    for i in columns:
        
        condition_color_list.append(condition_color_dict[sample_info_df.loc[i,'Sample_condition']])
            
        batch_color_list.append(batch_color_dict[sample_info_df.loc[i,'MS_run']])
        
    color_df = pd.DataFrame([batch_color_list,condition_color_list],index=["MS_run","Sample_condition"],columns=columns)
    return color_df.T


def correlation(df,method="pearsonr"):
    correlation_array=[]
    for i in df.columns:
        cc_list=[]
        for j in df.columns:
            if i==j:
                cc_list.append(1)
            else:
                two_columns_df = df.loc[:,[i,j]]
                #detected_index =two_columns_df.index[(1-(np.isnan(two_columns_df[i]) | np.isnan(two_columns_df[j])).values).astype(np.bool)]
                new_two_columns_df = two_columns_df.dropna(how="any")
                
                
                corr = pearsonr(new_two_columns_df.loc[:,i],new_two_columns_df.loc[:,j])[0]
                cc_list.append(corr)
        correlation_array.append(cc_list)
    correlation_array = pd.DataFrame(correlation_array,index=df.columns,columns=df.columns)
    correlation_array.to_csv(outdir+"/pearsonr_correlation_coefficient_of_z_statistic.txt",sep="\t")
    return correlation_array


def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass


def cluster_(correlation_df,color_df):
    patch_list = []

    for key in condition_color_dict.keys():
        patch = mpatches.Patch(color=condition_color_dict[key],label=key)
        patch_list.append(patch)

    plt.figure(figsize=(3,3))
    g = sns.clustermap(correlation_df,cmap = "bwr",linewidths = 0.5,
                   cbar_kws = {"shrink": .5,'label':"PCC"},col_colors = color_df,
                   xticklabels= False,yticklabels=True)
    g.fig.suptitle("Hierarchical clustering based on the PCC between each pair of samples",size=20)
    g.ax_row_dendrogram.set_visible(False)
    legends = patch_list
    l2=g.ax_heatmap.legend(loc='center left',bbox_to_anchor=(0.7,1.25),handles=legends,frameon=True)
    l2.set_title(title='Sample condition',prop={'size':6})

    
    plt.savefig(outdir + "/pearsonr_correlation_coefficient_of_z_statistic.pdf",dpi = 300,bbox_inches="tight")
    plt.savefig(outdir + "/pearsonr_correlation_coefficient_of_z_statistic.png",dpi = 300,bbox_inches="tight")
    plt.close("all")
def pc_scatter(pc_df,pca,batch_color_dict,condition_marker_dict,sample_info_df,outdir):

    batch_color_list = []
    condition_shape_list =[]
    for index in pc_df.index:

        condition = sample_info_df.loc[index,"Sample_condition"]
        batch = sample_info_df.loc[index,'MS_run']
        condition_shape_list.append(condition_marker_dict[condition])
        batch_color_list.append(batch_color_dict[batch])

    pc_number = len(pc_df.columns)
    for x in range(pc_number-1):
        for y in range(x+1,pc_number):

            x_pc = pc_df.columns[x]
            y_pc = pc_df.columns[y]

            x_pos = pc_df[x_pc]
            y_pos = pc_df[y_pc]

            fig = plt.figure(figsize=(5,5),dpi=300)
            ax = plt.gca()

            for i in range(len(pc_df)):
                ax.scatter(x_pos[i], y_pos[i], color=batch_color_list[i], marker=condition_shape_list[i])

            color_handles = [mpatches.Patch(color=list(batch_color_dict.values())[i]) for i in range(len(batch_color_dict))]
            shape_handles = [mlines.Line2D([0], [0], marker=list(condition_marker_dict.values())[i], color='w', markerfacecolor='black', markersize=10) for i in range(len(condition_marker_dict))]


            color_legend = plt.legend(color_handles, list(batch_color_dict.keys()), loc='upper left', title='MS run',bbox_to_anchor=(1.05, 0.5),ncol=2)
            shape_legend = plt.legend(shape_handles, list(condition_marker_dict.keys()), loc='upper left', title='Condition',bbox_to_anchor=(1.05, 1))


            ax.add_artist(color_legend)
            ax.add_artist(shape_legend)



            ax.set_xlabel(x_pc+"(%s%%)"%("%.1f"%(pca.explained_variance_ratio_[x]*100)),size=6)
            ax.set_ylabel(y_pc+"(%s%%)"%("%.1f"%(pca.explained_variance_ratio_[y]*100)),size=6)

            plt.savefig(outdir+"/%s_%s_scatterplot.pdf"%(x_pc,y_pc),dpi=300,bbox_inches="tight")

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--z_statistic_matrix",help="z-transformed log2-intensity",required=True)
    parser.add_argument("--sample_info",help="sample_information",required=True)
    parser.add_argument("--outdir",help="outdir for result",required=True)
    
    argv = vars(parser.parse_args())
    z_statistic_matrix = argv['z_statistic_matrix']
    sample_info = argv['sample_info']
    outdir = argv['outdir']

    mkdir(outdir)
    z_df=pd.read_csv(z_statistic_matrix,sep="\t",index_col=0)
    sample_info_df = pd.read_csv(sample_info,sep="\t",index_col=0)

    overlap_sample = list(set(z_df.columns)&set(sample_info_df.index))
    sample_info_df = sample_info_df.loc[overlap_sample]    

    z_df = z_df[overlap_sample]
    #Sample-based hierarchical cluster
    condition_color_dict = get_batch_color_dict(sample_info_df,col_="Sample_condition")
    batch_color_dict = get_batch_color_dict(sample_info_df,col_="MS_run")    
    correlation_df = correlation(z_df,method="pearsonr")
    color_df = get_color_list(correlation_df.index,sample_info_df,batch_color_dict,condition_color_dict)    
    cluster_(correlation_df,color_df)

    #principal component analysis
    z_df=z_df.dropna(how="any")
    condition_marker_dict = get_condition_mark_dict(sample_info_df,col_="Sample_condition")
    z_df=z_df.dropna(how="any")
    pca,pc_df = pca_analysis(z_df)
    pc_df.to_csv(outdir+"/pca_df.txt",sep="\t")
    pc_scatter(pc_df,pca,batch_color_dict,condition_marker_dict,sample_info_df,outdir)
