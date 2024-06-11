# -*- coding: utf-8 -*-
#"""
#Created on Wed Jul 24 09:26:15 2019

#@author: guixiuqi
#E-mail:  guixiuqi@gmail.com
#"""

#import gseapy
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact,hypergeom,binom
import numpy as np
import seaborn as sns

#go_bp_pathway_file = open("H:\project\enrich_dataset\enrichr.GO_Biological_Process_2018.gmt","r")
def get_kegg_pathway_information(gmt_file,all_detected_gene):
    pathway_dict={}
    with open(gmt_file) as f:
        for line in f:
            line_elements=line.strip().split("\t")
            #print (line_elements)

            kegg_term = line_elements[0]
            for gene in line_elements[1:]:
                if gene == "":
                    pass
                else:
                    if kegg_term in pathway_dict.keys():
                        if gene in all_detected_gene:
                            pathway_dict[kegg_term].append(gene)
                        else:
                            pass
                    else:
                        if gene in all_detected_gene:
                            pathway_dict[kegg_term]=[gene]
                        else:
                            pass
    return pathway_dict

def get_go_pathway_information(gmt_file,all_detected_gene):
    pathway_dict={}
    with open(gmt_file) as f:
        for line in f:
            line_elements =  line.strip().split("\t")
            #print (line_elements)
            go_term = line_elements[0]
            for gene in line_elements[1:]:
                if gene =="":
                    pass
                else:
                    if go_term in pathway_dict.keys():
                        if gene in all_detected_gene:
                            pathway_dict[go_term].append(gene)
                        else:
                            pass
                    else:
                        if gene in all_detected_gene:
                            pathway_dict[go_term]=[gene]
                        else:
                            pass
    return pathway_dict                            
#######################################################################
#                  target_path             not in target path      
#interest              a                       b            a+b(number of interest genes)
#
#not interest          c                       d             c+d (number 0f [all_gene - interest gene])
#         number of target pathway gene       number of [allgene-target_gene_number]
#######################################################################
def get_pathway_database_gene_number(pathway_dict):
    pathway_gene_list=[]
    for key in pathway_dict.keys():
        for gene in pathway_dict[key]:
            if gene in pathway_gene_list:
                pass
            else:
                pathway_gene_list.append(gene)
    return pathway_gene_list,len(pathway_gene_list)

        
def enrich(pathway_dict,gene_list,all_detected_gene):
    pathway_database_gene_list,pathway_database_gene_number = get_pathway_database_gene_number(pathway_dict)
#    print (pathway_database_gene_number)
    enrich_result=[]
    for key in pathway_dict.keys():
        target_pathway_gene_list = pathway_dict[key]  ###all genes in path
        a = len(set(gene_list)&set(target_pathway_gene_list)) ###insterest gene numbers in target path  
        b = len(gene_list) - a   ### insterest gene numbers not in target path
        c = len(target_pathway_gene_list) - a
        d = len(all_detected_gene) - (a+b) - c
        
        if a > 0:
            pvalue = fisher_exact([[a,b],[c,d]],alternative="greater")[1]
            genes=",".join(list(set(gene_list)&set(target_pathway_gene_list)))
            enrich_result.append([key,str(a)+"/"+str(len(target_pathway_gene_list)),genes,pvalue])
        else:
            pass
    enrich_result_df = pd.DataFrame(enrich_result,columns=["pathway","overlap","gene list","pvalue"])
    enrich_result_df.drop(columns=["pathway"])
    sorted_df = enrich_result_df.sort_values(by=["pvalue"])
    padjust_value = [len(enrich_result_df)*i if len(enrich_result_df)*i<=1 else 1  for i in sorted_df["pvalue"]]
    sorted_df["padjust_by_bonferroni"]=padjust_value
    sorted_df["padjust_by_BH"] =[i if i<1 else 1 for i in sorted_df["pvalue"]/sorted_df["pvalue"].rank(axis=0)*len(sorted_df)]
   # sorted_df = sorted_df.sort_values(by="padjust_by_BH")
    sorted_df.index = np.arange(len(sorted_df))
    return sorted_df

def enrich_plot(kegg_enrich_result_adjust_df,title="",pvalue=0.05,fig_file = "",top = 20):
    
    significant_df = kegg_enrich_result_adjust_df.iloc[0:20,:].loc[kegg_enrich_result_adjust_df["padjust_by_BH"]<=pvalue]
    if len(significant_df) >0:
        plt.figure(figsize=(5,len(significant_df)/2 ))
        #plt.scatter(-np.log(significant_df["padjust_by_BH"].values),[i for i in range(len(significant_df))],
        #            s = significant_df.apply(lambda x:int(x["overlap"].split("/")[0]),axis = 1).values*3,
        #            c = -np.log(significant_df["padjust_by_BH"].values),cmap="Reds")
        plt.barh([i for i in range(len(significant_df))],-np.log10(significant_df["padjust_by_BH"].values))
        
        plt.xlabel('$-log_{10}{ (padjust\ by\ BH)}$',size=20)
        plt.yticks([i for i in range(len(significant_df))],significant_df["pathway"] + " ["+significant_df["overlap"]+"]",size=15)
        plt.ylim(len(significant_df)-0.5,-0.5)
        #cbar=plt.colorbar()
        #cbar.ax.tick_params(labelsize=10)
        #cbar.solids.set_edgecolor("face")
        #cbar.set_label('$-log_{10}{(padjust\ by\ BH)}$',size=20)
        plt.title(title,size=20)
        plt.savefig(fig_file,dpi=300,bbox_inches="tight")
        plt.close("all")
    else:
        pass
 
def mean_heatmap(mean_df,gene_list,title = "",fig_file = "" ):
    plt.figure(figsize=(6,10))
    plt.title(title,size=15)
    new = mean_df.loc[gene_list]
    sns.heatmap(new,cmap="bwr",vmin=-8,vmax=8,xticklabels = ["A0","A6","F0","F6"], cbar_kws = {"shrink": 0.6})
    plt.savefig(fig_file,dpi=300,bbox_inches="tight")
    plt.close("all")
    

def get_mean_df(table):
    mean_column_list = [] 
    for i in table.columns:
        if "mean" in i:
            mean_column_list.append(i)
    return table[mean_column_list]
    
    

def get_z_transformed_columns(table):
    log_two_intensity = []
    z_transformed = []
    z_score=[]
    for i in table.columns:
        if "Z-score" in i :
            z_score.append(i)
        elif "Z-tranformed" in i  and "mean" not in i:
            z_transformed.append(i)
        elif "log2-intensity" in i and "mean" not in i :
            log_two_intensity.append(i)
        else:
            pass
    return z_transformed
    

def rearrangement(z_transformed):
    A0=[]
    F0=[]
    A6=[]
    F6=[]
    for i in z_transformed:
        if "A0" in i:
            A0.append(i)
        elif "A6" in i:
            A6.append(i)
        elif "F0" in i:
            F0.append(i)
        else:
            F6.append(i)
    return A0,A6,F0,F6,A0+A6+F0+F6

def get_gene_detectd_more_than_2_each_condition(table):
    keep_gene_list = []
    remove_gene_list = []
    z_transformed = get_z_transformed_columns(table)
    A0,A6,F0,F6,rearrangement_z_transformed = rearrangement(z_transformed)
    for gene in table.index:
        num = 0
        for condition in [A0,A6,F0,F6]:
            if list(table.loc[gene,condition].values).count("Not detected") > 3:
                num+=1
        if num > 0 :
            remove_gene_list.append(gene)
        else:
            keep_gene_list.append(gene)
    return keep_gene_list

"""
kegg_pathway_file = r"H:\project\enrich_dataset\enrichr.KEGG_2016.gmt"
go_bp_pathway_file = r"H:/project/enrich_dataset/human_GO_Biological_Process_2015.gmt"
go_cc_pathway_file = r"H:/project/enrich_dataset/human_GO_Cellular_Component_2015.gmt"
go_mf_pathway_file = r"H:/project/enrich_dataset/human_GO_Molecular_Function_2015.gmt" 
            
file=r'H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\statistics_all_comparisons_saved.xls'
data=pd.ExcelFile(file)
table=data.parse('Sheet1')
table.index=table.loc[:,'Gene Name']

gene_detectd_more_than_2_each_condition_list = get_gene_detectd_more_than_2_each_condition(table)

print("gene number detectd more than 2 each condition list :" +str(len(gene_detectd_more_than_2_each_condition_list)))
kegg_pathway_dict = get_kegg_pathway_information(kegg_pathway_file,gene_detectd_more_than_2_each_condition_list)
go_pathway_bp_dict = get_go_pathway_information(go_bp_pathway_file,gene_detectd_more_than_2_each_condition_list)
go_pathway_cc_dict = get_go_pathway_information(go_cc_pathway_file,gene_detectd_more_than_2_each_condition_list)
go_pathway_mf_dict = get_go_pathway_information(go_mf_pathway_file,gene_detectd_more_than_2_each_condition_list)


mean_df = get_mean_df(table)

def all_mean_heatmap(mean_df,gene_list,title = ""):
   # plt.figure(figsize=(6,10))
   # plt.title(title,size=15)
    new = mean_df.loc[gene_list]
    sns.heatmap(new,cmap="bwr",vmin=-10,vmax=10,xticklabels = ["A0","A6","F0","F6"], cbar_kws = {"shrink": 0.6})
   # plt.savefig(fig_file,dpi=300,bbox_inches="tight")
  

cluster_list = ["cluster"+str(i) for i in range(1,13)] 

num=1
plt.figure(figsize=(20,15))
for cluster in cluster_list:
    cluster_file = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\%s_all_hypervariable_gene.xls"%(cluster)
    data = pd.read_csv(cluster_file,sep = ",")
    cluster_gene = list(data["Gene Name"])
    cluster = cluster_file.split("\\")[-1].split("_")[0]    
    heatmap_file = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\mean.heatmap.pdf"
    plt.subplot(3,4,num)
    plt.title(cluster)
    all_mean_heatmap(mean_df,cluster_gene,title = cluster)
    num += 1
plt.savefig(heatmap_file,dpi=300)
    



import glob
cluster_files = glob.glob(r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\*_all_hypervariable_gene.xls")
### fisher exact test
for cluster_file in cluster_files:

    data = pd.read_csv(cluster_file,sep = ",")
    cluster_gene = list(data["Gene Name"])
    cluster = cluster_file.split("\\")[-1].split("_")[0]
    heatmap_file = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\%s.heatmap.pdf"%(cluster)
    heatmap_file1 = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\%s.heatmap.png"%(cluster)
    mean_heatmap(mean_df,cluster_gene,title = cluster,fig_file=heatmap_file)
    mean_heatmap(mean_df,cluster_gene,title = cluster,fig_file=heatmap_file1)
    str_list = ["kegg","go_bp","go_cc","go_mf"]
    str_list1 = ["KEGG","GO_Biological_Process","GO_Cellular_Component","GO_Molecular_Function"]
    i = 0
    for pathway_dict in [kegg_pathway_dict,go_pathway_bp_dict,go_pathway_cc_dict,go_pathway_mf_dict]:
        
        enrich_result = enrich(pathway_dict,cluster_gene,gene_detectd_more_than_2_each_condition_list)
        enrich_result_sort = enrich_result.sort_values(by=["padjust_by_BH"])
        #if len(enrich_result)>0:
        enrich_file = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\enrich_result\%s_%s_enrich_result.pdf"%(cluster,str_list[i])
        enrich_file1 = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\enrich_result\%s_%s_enrich_result.png"%(cluster,str_list[i])
        
        
        
        
        enrich_plot(enrich_result_sort,title=cluster + "_" + str_list1[i] ,pvalue=0.05,fig_file = enrich_file)
        enrich_plot(enrich_result_sort,title=cluster + "_" + str_list1[i],pvalue=0.05,fig_file = enrich_file1)
        
        
        
        enrich_result.to_csv(r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\enrich_result\%s_%s_enrich_result.txt"%(cluster,str_list[i]),sep = "\t")
        i+=1


#####gseapy enrich result
for cluster_file in cluster_files:
    
    data = pd.read_csv(cluster_file,sep=",")
    cluster_gene = list(data["Gene Name"])
    cluster = cluster_file.split("\\")[-1].split("_")[0]
    str_list = ["kegg","go_bp","go_cc","go_mf"]
    i = 0
    for pathway_dict in [kegg_pathway_dict,go_pathway_bp_dict,go_pathway_cc_dict,go_pathway_mf_dict]:
        enr = gseapy.enrichr(gene_list = cluster_gene,gene_sets = pathway_dict,format = "png",top_term = 20,
                             outdir = r"H:\project\protein\4Xiuqi\Coding_proteomics_zMAP_intensity\guixiuqi_results\gene_cluster\cluster_gene_all_hypervariable_include_not_detected\gseapy_enrich_result\%s\%s"%(cluster,str_list[i]))
        i += 1        

"""


