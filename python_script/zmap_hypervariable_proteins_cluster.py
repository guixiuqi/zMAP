# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 10:12:48 2020

@author: guixiuqi
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

#import Cluster_Ensembles as CE
from sklearn.cluster import KMeans
from collections import defaultdict
import os
import seaborn as sns
import argparse
from sklearn.decomposition import PCA 
from matplotlib.lines import Line2D
import scipy

import random
import matplotlib.pyplot as plt
from collections import Counter

def read_file(filein):
    filein_df = pd.read_csv(filein,sep = "\t",index_col = 0)
    return filein_df

def get_hypervariable_df(filein_df,pvalue_cutoff = 0.05):
    
    hypervariable_df = filein_df.loc[filein_df["BH-corrected combined pvalue"] < pvalue_cutoff]
    
    return hypervariable_df


def get_hyper_z_trans_df_mean(hypervariable_df,rearrange_z_trans_df,sample_info_df):
    
    hyper_z_df = rearrange_z_trans_df.loc[hypervariable_df.index]
    mean_hypervariable_df = pd.DataFrame(index = hyper_z_df.index)

    for condition in set(list(sample_info_df['Sample_condition'].values)):
        sample_list = list(sample_info_df.loc[sample_info_df['Sample_condition']==condition].index)
        mean_hypervariable_df["average %s"%condition] = rearrange_z_trans_df.loc[:,sample_list].mean(axis=1)
    
    return mean_hypervariable_df
    
    
def gene_cluster_by_mean_value(mean_hypervariable_df,rearrange_z_trans_df,filein_df,outdir,vmin_value = -6,vmax_value = 6,
                               cluster_number = 8,cutoff = 50,pvalue_cutoff = 0.01):
    
    ## cluster
    map_color = "bwr"
    g1=sns.clustermap(mean_hypervariable_df,cmap=map_color,metric='correlation',col_cluster = False,
                 figsize=(3,12),vmin = vmin_value,vmax = vmax_value,xticklabels = mean_hypervariable_df.columns)
    
    proteins_indexs=scipy.cluster.hierarchy.fcluster(g1.dendrogram_row.calculated_linkage,
                                                     t = cluster_number,criterion='maxclust')
    
    results=Counter(proteins_indexs)
    sub_protein_clusters={}
    
    for protein ,cluster in zip(mean_hypervariable_df.index,proteins_indexs):
        if cluster in sub_protein_clusters.keys():
            sub_protein_clusters[cluster].append(protein)
        else:
            sub_protein_clusters[cluster] = [protein]
            
    cluster_list = []
    combine_cluster = []
    for cluster in sub_protein_clusters.keys():
        if len(sub_protein_clusters[cluster]) > cutoff:
            cluster_list.append(cluster)
        else:
            combine_cluster.append(cluster)
    
    cluster_list.append(combine_cluster)
    gene_list = []
    for cluster in cluster_list:
        if isinstance(cluster,list):
            for small_cluster in cluster:
                gene_list += sub_protein_clusters[small_cluster]
        else:
            gene_list += sub_protein_clusters[cluster]
            
    ##cluster_heatmap
    plt.figure(figsize=(5,15))
    plt.box(False)
    c_map = "coolwarm"        
    vmin = -6
    vmax=6
    heat = plt.imshow(rearrange_z_trans_df.loc[gene_list],aspect="auto",cmap=c_map,vmin=vmin, vmax=vmax,interpolation='nearest')
    #sns.clustermap(rearrange_z_trans_df.loc[gene_list],col_cluster=False,row_cluster=False,vmin=vmin,vmax=vmax)
#    plt.xticks(np.arange(len(rearrange_z_trans_df.columns)), labels=rearrange_z_trans_df.columns)
#    plt.yticks(np.arange(len(gene_list)), labels=gene_list)
    a = -0.5    
    extend = 10
    x_value = len(rearrange_z_trans_df.columns) - 0.2    
    text_value = len(rearrange_z_trans_df.columns)
    
    i = 1
    for key in cluster_list:
        if isinstance(key,list):
            b=0
            for j in key:
                b+=len(sub_protein_clusters[j])
            plt.plot([x_value,x_value],[a+extend,a+b-extend],linewidth=5)
            plt.text(text_value,a+b/2,"cluster_0")
            a=a+b          
        else:
            b=len(sub_protein_clusters[key])
            plt.plot([x_value,x_value],[a+extend,a+len(sub_protein_clusters[key])-extend],
                     linewidth=5)
            plt.text(text_value,a+b/2,"cluster_%s"%(i))
            
            i +=1
            a = a+len(sub_protein_clusters[key])
    plt.xticks(range(0,len(rearrange_z_trans_df.columns)),rearrange_z_trans_df.columns,rotation=90)
    plt.yticks([])
    cbar = plt.colorbar(heat, shrink=0.5,orientation="horizontal")
    cbar.set_ticks([-6,-3,0,3,6])
    cbar.set_ticklabels([-6,-3,0,3,6])    
    plt.savefig(outdir + "/hypervariable_proteins_cluster_heatmap_bh_cutoff_%s.pdf"%pvalue_cutoff,dpi=300,bbox_inches="tight")
    plt.savefig(outdir + "/hypervariable_proteins_cluster_heatmap_bh_cutoff_%s.png"%pvalue_cutoff,dpi=300,bbox_inches="tight")
    
   # mkdir(outdir+"/clusters")
    new_dir = outdir+"/clusters"
    mkdir(new_dir)
    ## cluster results
    i=1
    for cluster in cluster_list:
        
        if isinstance(cluster, list):
            sub_gene_list = []
            
            for small_cluster in cluster:
                sub_gene_list += sub_protein_clusters[small_cluster]
                
            cluster_df = filein_df.loc[sub_gene_list]
            cluster_df.to_csv(new_dir+ "/" +"cluster_0.txt",sep="\t")
            
        else:
            cluster_df = filein_df.loc[sub_protein_clusters[cluster]]
            cluster_df.to_csv(new_dir + "/" + "cluster_%s.txt"%(i),sep="\t")
            i+=1
    return sub_protein_clusters,cluster_list

import enrichment_protein

def pathway_background(filein_df):
    gmt_dir = os.getenv('gmtPATH')
    kegg_pathway_file = "%s/enrichr.KEGG_2016.gmt"%gmt_dir
    go_bp_pathway_file = "%s/human_GO_Biological_Process_2015.gmt"%gmt_dir
    go_cc_pathway_file = "%s/human_GO_Cellular_Component_2015.gmt"%gmt_dir
    go_mf_pathway_file = "%s/human_GO_Molecular_Function_2015.gmt" %gmt_dir
    kegg_pathway_dict = enrichment_protein.get_kegg_pathway_information(kegg_pathway_file,filein_df.index)
    go_pathway_bp_dict = enrichment_protein.get_go_pathway_information(go_bp_pathway_file,filein_df.index)
    go_pathway_cc_dict = enrichment_protein.get_go_pathway_information(go_cc_pathway_file,filein_df.index)
    go_pathway_mf_dict = enrichment_protein.get_go_pathway_information(go_mf_pathway_file,filein_df.index)
    return [kegg_pathway_dict,go_pathway_bp_dict,go_pathway_cc_dict,go_pathway_mf_dict]


def pathway_analysis(sub_protein_clusters,cluster_list,outdir,pathway_dict_list,filein_df):
    
    enrich_dir = outdir + "/enrichment_results/"
    mkdir(enrich_dir)
    str_list = ["kegg","go_bp","go_cc","go_mf"]
    str_list1 = ["KEGG","GO_Biological_Process","GO_Cellular_Component","GO_Molecular_Function"]
    
    for cluster in cluster_list:
        if isinstance(cluster,list):
            pass
        else:
            sub_gene_list = sub_protein_clusters[cluster]
            i = 0
            for pathway_dict in pathway_dict_list:
                enrich_result = enrichment_protein.enrich(pathway_dict,sub_gene_list,filein_df.index)
                enrich_result_sort = enrich_result.sort_values(by=["padjust_by_BH"])
        
                enrich_file = enrich_dir + "cluster%s_%s_enrich_result.pdf"%(cluster,str_list[i])
                enrich_file1 = enrich_dir + "cluster%s_%s_enrich_result.png"%(cluster,str_list[i])
        
                enrichment_protein.enrich_plot(enrich_result_sort,title= str(cluster) + "_" + str_list1[i] ,pvalue=0.05,fig_file = enrich_file)
                enrichment_protein.enrich_plot(enrich_result_sort,title= str(cluster) + "_" + str_list1[i],pvalue=0.05,fig_file = enrich_file1)
        
                enrich_result.to_csv(enrich_dir+"cluster%s_%s.txt"%(cluster,str_list[i]),sep = "\t")
                i+=1



def top_gene(rearrange_z_trans_df,hypervariable_df,mean_hypervariable_df,outdir,top_gene = 100,cluster_num = 8):
    top_hypervariable_dir = outdir +"/top_differentially_expressed_proteins_%s"%(top_gene)
    mkdir(top_hypervariable_dir)
    mean_hypervariable_df_copy = mean_hypervariable_df.copy()
    mean_hypervariable_df_copy["standard deviation"] = mean_hypervariable_df_copy.std(axis=1)
    hypervariable_df_sort =  mean_hypervariable_df_copy.sort_values(by="standard deviation")
    top_gene_list = mean_hypervariable_df_copy.iloc[0:top_gene].index
    top_gene_df = rearrange_z_trans_df.loc[top_gene_list]
    top_mean_df = mean_hypervariable_df.loc[top_gene_list]
    
    estimator = KMeans(n_clusters=cluster_num,random_state =222)
    estimator.fit(top_mean_df)
    label_pred = estimator.labels_
    #color=["lightcoral","darkorange","yellowgreen","lime","lightseagreen","grey","lightskyblue","blue","black"]
    #m=331
    #plt.figure(figsize=(4,15))
    protein = []
    for i in range(0,cluster_num):
        cluster_df = top_mean_df[label_pred == i]
        if len(cluster_df)>2:
            sort_rank_dict = defaultdict(list)
            for gene in cluster_df.index:
                expression_value = list(cluster_df.loc[gene].values)
                sorted_expression_value = sorted(expression_value)
                str_rank=""
                for j in sorted_expression_value[-2:]:
                    str_rank += str(expression_value.index(j))
                sort_rank_dict[str_rank].append(gene)
            for key in sort_rank_dict.keys():
                protein +=sort_rank_dict[key]
        else:                    
            protein += list(top_mean_df[label_pred == i].index)
    
    new_top_df = top_gene_df.loc[protein]
    new_top_df.to_csv(top_hypervariable_dir + "/top_%s_differentially_expressed_proteins.txt"%top_gene,sep="\t")
    plt.figure(figsize=(4,15))      
    sns.heatmap(new_top_df,cmap="coolwarm",vmin=-10,vmax=10,square=True,cbar_kws = {"shrink":0.2},xticklabels=True,yticklabels=True)
    plt.savefig(top_hypervariable_dir+"/top_%s_differentially_expressed_proteins_heatmap.pdf"%top_gene,dpi=300,bbox_inches="tight")    
    plt.savefig(top_hypervariable_dir+"/top_%s_differentially_expressed_proteins_heatmap.png"%top_gene,dpi=300,bbox_inches="tight")
    
   
def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass    
    
        
    
    
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--pvalue_results",help="pvalue_results",required=True)
    parser.add_argument("--z_statistic_matrix",help="z_statistic_matrix",required=True)
    parser.add_argument("--sample_info",help="sample_information",required=True)
    parser.add_argument("--outdir",help="outdir",required=True)
    parser.add_argument("--fdr",type=float, default=0.05,help="BH-corrected pvalue cutoff used to define hypervariable proteins,default is 0.05")
    parser.add_argument("--cluster_number_for_hypervariable",help="cluster number of hypervariable proteins defined by user",required=True)
    parser.add_argument("--minclustersize",help="minimum protein number for each cluster",required=True)
    parser.add_argument("--top",help="number of top hypervariable proteins",required=True)
    parser.add_argument("--cluster_number_for_top_proteins",help="cluster number for top hypervariable proteins",required=True)
    argv = vars(parser.parse_args())
    pvalue_results = argv['pvalue_results']
    z_statistic_matrix = argv['z_statistic_matrix']
    sample_info = argv['sample_info']
    outdir = argv['outdir']
    mkdir(outdir)
    bh_pvalue_cutoff = float(argv['fdr'])
    cluster_num = int(argv['cluster_number_for_hypervariable'])
    minclustersize = int(argv['minclustersize'])
    top = int(argv['top'])
    cluster_number2 = int(argv["cluster_number_for_top_proteins"])    
    pvalue_results_df = read_file(pvalue_results) 

    hypervariable_df = get_hypervariable_df(pvalue_results_df,pvalue_cutoff=bh_pvalue_cutoff)

    z_trans_df = pd.read_csv(z_statistic_matrix,sep="\t",index_col=0)
    sample_info_df = pd.read_csv(sample_info,sep="\t",index_col=0)

    sample_info_df_sort = sample_info_df.sort_values(by="Sample_condition")

    rearrange_z_trans_df = z_trans_df[list(sample_info_df_sort.index)]
    
    mean_hypervariable_df  = get_hyper_z_trans_df_mean(hypervariable_df,rearrange_z_trans_df,sample_info_df)
    
    print("%s hypervariable proteins were identified(FDR<0.05).Hypervariable proteins will be clustered into %s clusters."%(str(len(mean_hypervariable_df)),str(cluster_num)))    
    sub_protein_clusters,cluster_list = gene_cluster_by_mean_value(mean_hypervariable_df,rearrange_z_trans_df,pvalue_results_df,outdir,vmin_value = -6,vmax_value = 6,
                               cluster_number = cluster_num,cutoff = minclustersize,pvalue_cutoff=0.05)
    
    pathway_dict_list = pathway_background(pvalue_results_df)
    
    print("Start pathway enrichment analysis for each cluster.")
    pathway_analysis(sub_protein_clusters,cluster_list,outdir,pathway_dict_list,pvalue_results_df)
    
    print("Top %s hypervariable proteins will be clustered into %s clusters."%(str(top),str(cluster_number2)))
    top_gene(rearrange_z_trans_df,hypervariable_df,mean_hypervariable_df,outdir,top_gene = top,cluster_num = cluster_number2)
    print("All finished.")
    
    
    
    
    
    
    
    
    
    
    
    












