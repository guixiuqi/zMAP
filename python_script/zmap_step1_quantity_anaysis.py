# -*- coding: utf-8 -*-
#"""
#Created on Fri Oct 11 18:08:21 2019
#
#@author: guixiuqi
#E-mail:  guixiuqi@gmail.com
#"""
#
import scipy
import numpy as np
import matplotlib
matplotlib.use('Agg')
#import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#matplotlib.use('Agg')
#import matplotlib
import os
import math
import scipy.stats as stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import glob
import argparse
import pandas as pd
import natural_cubic_model
#matplotlib.use('Agg')
####################log2 normalization##################

def remove_bad_gene(raw_intensity_df,cutoff=10):
    good_index = []
    for index in raw_intensity_df.index:
        index_df = raw_intensity_df.loc[index] > cutoff
        num = sum(index_df.values)
        if num > 0:
            good_index.append(index)
        else:
            pass
    raw_intensity_df = raw_intensity_df.loc[good_index]
    return raw_intensity_df

def normalization(raw_intensity_df,output_dir,l=1.5,batch=""):
    #raw_intensity_df = pd.read_csv(raw_intensity_table,sep="\t",index_col = "Gene Name")
    sum_trimmed_intensity_list = []
    for column in raw_intensity_df.columns:
        q3 = np.percentile(raw_intensity_df[column],75)
        q1 = np.percentile(raw_intensity_df[column],25)
        top_threshold_intensity = q3 + l*(q3 - q1)
        trimmed_intensity_df = raw_intensity_df[column].loc[raw_intensity_df[column]<=top_threshold_intensity]
        sum_trimmed_intensity = trimmed_intensity_df.sum(axis=0)
        sum_trimmed_intensity_list.append(sum_trimmed_intensity)
    mean_of_sum_tirmmed_intensity = np.mean(sum_trimmed_intensity_list)
    size_factor = np.array(sum_trimmed_intensity_list)/mean_of_sum_tirmmed_intensity
    normalized_intensity_df = raw_intensity_df.div(size_factor)
    log2_normalized_intensity_df = np.log2(normalized_intensity_df + 1)
    log2_raw_intensity_df = np.log2(raw_intensity_df + 1)
    
    plt.figure(figsize=(8,8))
    
    ax = plt.subplot(211)
    ax1 = plt.subplot(212)

    boxprops = dict(linestyle = '-', linewidth = 2, color = 'darkorange')
    medianprops = dict(linestyle = '-', linewidth = 1.5, color = 'firebrick')
    flierprops = dict(marker = 'o', markerfacecolor = 'grey',
                  linestyle = 'none')
    whiskerprops = dict(linestyle = '-', linewidth = 2, color = 'grey')

    capprops = dict(linestyle='-', linewidth=2, color='grey')
    
    ax.boxplot(log2_raw_intensity_df.values,labels= log2_raw_intensity_df.columns,
                vert = False,widths=0.7,
            boxprops = boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            capprops=capprops,flierprops=flierprops)
    ax.set_xlabel("Raw log$_2$-intensity",size=15)
    ax.set_ylabel("Sample",size=15)
	
    
    ax1.boxplot(log2_normalized_intensity_df.values,labels= log2_normalized_intensity_df.columns,
                vert = False,widths=0.7,
            boxprops = boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            capprops=capprops,flierprops=flierprops)
    ax1.set_xlabel("Normalized log$_2$-intensity",size=15)
    ax1.set_ylabel("Sample",size=15)
   # plt.title(batch,size=20)
    plt.suptitle(batch,size=20)
    plt.savefig(output_dir+r"/"+"%s_protein_intensity_boxplot.png"%(batch),dpi=300,bbox_inches="tight")
    plt.savefig(output_dir+r"/"+"%s_protein_intensity_boxplot.pdf"%(batch),dpi=300,bbox_inches="tight")
    plt.close('all')
    return log2_normalized_intensity_df

######################## mkdir ################################

def mkdir(path_dir):
    isexists = os.path.exists(path_dir)
    if not isexists:
        os.makedirs(path_dir)
    else:
        pass
		
########### calculate mean log2-intensity among samples and observed variance#####
def calculate_a_variance_value(log2_normalized_intensity_df):
    #log2_normalized_intensity_df = pd.read_csv(new_table_txt,sep="\t", index_col = 'Gene Name')
    sample_num = len(log2_normalized_intensity_df.columns)
    A_value_list = log2_normalized_intensity_df.mean(axis=1).values
    variance_list = log2_normalized_intensity_df.var(axis=1).values
    log2_normalized_intensity_df["A-value"] = A_value_list
    log2_normalized_intensity_df["Observed Variance"] = variance_list
    return log2_normalized_intensity_df,sample_num

################## sort log2-intensity table by A-value###########################
def get_all_sliding_region(table):
    protein_num = len(table)
    sorted_table = table.sort_values(by = "A-value")
 
    return sorted_table


################################## Michael and Schucany ##########################
def MS_plotting_position(start_pos, end_pos, total_pos, aa):
    ### start_Pos=1 end_pos = window size+1 total_pos = window size aa =0.375
    from functools import reduce
    constant = (total_pos - aa + 1) / (total_pos - 2*aa + 1)
    
    def get_inner(x):
        r_digital = (x - aa) / (x - aa + 1)
        return r_digital

    pp_list = []
    for j in range(start_pos, end_pos):
        p_i_list = map(get_inner, range(j, total_pos+1))
        p_i = constant * reduce(lambda m,n: m*n, p_i_list)
        pp_list.append(p_i)
    return pp_list

#####################linear regression for sliding window#########################
def window_vaplot(sorted_table,num_of_sample,window_size = 400,sliding_size = 100,result_dir="",top_gene_used_to_fit_linear_model = 0.3,batch=""):
    
    sliding_region_left = 0
    sliding_region_right = len(sorted_table)
    
  #  vaplot_dir = result_dir + "/VAplot"
  #  linear_regression = result_dir + "/linear_model"
    mkdir(result_dir)
    #mkdir(vaplot_dir)
    #mkdir(linear_regression)
    window_mean_A_list = []
    window_slope_list = []
    r2_list = []
    left = sliding_region_left
    while left < sliding_region_right:
        if left + window_size > sliding_region_right:
            non_target_window_df = pd.concat([sorted_table.iloc[0:left,:],sorted_table.iloc[sliding_region_right:,:]],axis=0)
        else:
            non_target_window_df = pd.concat([sorted_table.iloc[0:left,:],sorted_table.iloc[left+window_size:,:]],axis=0)
        
        target_window_df = sorted_table.iloc[left:min(left+window_size,sliding_region_right),:]
        
        target_window_df_variance = target_window_df["Observed Variance"]
        
        gene_number_in_window = len(target_window_df_variance)
        
        if gene_number_in_window > 100:
            
            top_percent_num = int(gene_number_in_window * top_gene_used_to_fit_linear_model)
        
            sorted_target_window_df_variance = target_window_df_variance.sort_values(ascending=True)
        
            gene_used_to_fit_curve = sorted_target_window_df_variance.index[0:top_percent_num]
        
            other_genes = sorted_target_window_df_variance.index[top_percent_num:]
        
            target_window_bottom_gene_df = target_window_df.loc[gene_used_to_fit_curve]
        
            target_window_up_gene_df = target_window_df.loc[other_genes]
        
            window_mean_A = np.mean(target_window_df["A-value"])
        
            window_mean_A_list.append(window_mean_A)
        
            degree_of_freedom = num_of_sample - 1
            
            pp_list = MS_plotting_position(1,gene_number_in_window+1,gene_number_in_window,0.375)

            theoretical_quantiles = [scipy.stats.chi2.ppf(j,degree_of_freedom)/degree_of_freedom for j in pp_list] 
        
            p = np.polyfit(theoretical_quantiles[0:top_percent_num],sorted_target_window_df_variance[0:top_percent_num],1)
            regr = LinearRegression()
            model_x = np.asarray(theoretical_quantiles[0:top_percent_num]).reshape(len(theoretical_quantiles[0:top_percent_num]),1)
            model_y = np.asarray(sorted_target_window_df_variance[0:top_percent_num]).reshape(len(sorted_target_window_df_variance[0:top_percent_num]),1)
            regr.fit(model_x,model_y)
            r2 = regr.score(model_x, model_y)
            r2_list.append(r2)
            window_slope_list.append(p[0])
            if p[1] > 0:
                p_function = "y= %s x + %s, r-square = %s" %(str(p[0]), str(p[1]), str(r2))
            elif p[1] < 0:
                p_function = "y= %s x - %s, r-square = %s" %(str(p[0]), str(-p[1]), str(r2))
            else:
                p_function = "y= %s x, r-square = %s" %(str(p[0]), str(r2))
        left += sliding_size
    mean_a_scope_df = pd.DataFrame([window_mean_A_list,window_slope_list],index = ["mean_A_value","estimated_variance"]).T
    mean_a_scope_df.to_csv(result_dir + "/"+"%s_mean_A_estimated_variance.txt"%(batch),sep ="\t")
    plt.figure(figsize= (5,6))
    plt.boxplot(r2_list)
    plt.ylabel("$R^2$")
    plt.title("Boxplot for %s $R^2$ \nof linear regression"%(batch),size=15)
    
    plt.savefig(output_dir + "/"+"%s_boxplot_for_r2_of_linear_regression.png"%(batch),dpi=300,bbox_inches="tight")
    
    plt.savefig(output_dir + "/"+"%s_boxplot_for_r2_of_linear_regression.pdf"%(batch),dpi=300,bbox_inches="tight")
    plt.close('all')
    
    return window_mean_A_list,window_slope_list


#############calculate R^2###################
def getIndexes(y_predict, y_data):
    
    n = y_data.size
    
    SSE = ((y_data-y_predict)**2).sum()
    
    MSE = SSE/n
    
    RMSE = np.sqrt(MSE)
    
    u = y_data.mean()
    SST = ((y_data -u)**2).sum()
    SSR = SST -SSE
    R_square = SSR/SST
    
    return SSE,MSE,RMSE,R_square


###########expotential_function##########
def func_nl_lsq(x,a,b,c):
    return np.exp(a + (b*x))+ c
def func_nl_lsp_no_constant_term(x,a,b):
    return np.exp(a + (b*x))


################################fit non_linear model###########################################
def fit_exp_model(table,mean_a_list,slope_list,used_for_model_fitting=[0.05,0.95],result_dir="",batch=""):
    
    #nonlinear_dir = result_dir+"/nonlinear_model"
    #mkdir(nonlinear_dir)
    
    a_value = np.asarray(table["A-value"].sort_values(ascending=True))
    
    start = int(len(mean_a_list)*used_for_model_fitting[0])
    end = int(len(mean_a_list)*used_for_model_fitting[1])
    x_model = np.asarray(mean_a_list[start:end])
    y_model = np.asarray(slope_list[start:end])
    popt, pcov = curve_fit(func_nl_lsq, x_model, y_model,bounds=([-1,-5,-1], [2, 0, 0.1]))
    if popt[2] > 0:
        SSE,MSE,RMSE,R_square = getIndexes(func_nl_lsq(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(mean_a_list[0:start]+mean_a_list[end:],slope_list[0:start]+slope_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value, func_nl_lsq(a_value, *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)+%s$"%(str("%.3f"%(popt[0])),
                                                                   str("%.3f"%(-popt[1])),str("%.3f"%(popt[2]))))
    else:
        popt, pcov = curve_fit(func_nl_lsp_no_constant_term, x_model, y_model,bounds=([-1,-5], [2, 0]))
        SSE,MSE,RMSE,R_square = getIndexes(func_nl_lsp_no_constant_term(x_model,*popt),y_model)
        plt.figure(figsize=(5,5))
        plt.scatter(mean_a_list[0:start]+mean_a_list[end:],slope_list[0:start]+slope_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
        plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )
        plt.plot(a_value, func_nl_lsp_no_constant_term(a_value, *popt), '--',linewidth = 3,color = "darkorange",
             label = "Fitted exponential funcation:\n$\sigma^2=exp(%s-%s*S)$"%(str("%.3f"%(popt[0])),str("%.3f"%(-popt[1]))))
        
    plt.xlabel("$Average\ log_2-intensity\ S$",size=15)
    plt.ylabel("$Estimated\ variance\ \sigma^2$",size=15)
    plt.title("%s $R^2$=%s"%(batch,str("%.3f"%(R_square))),fontsize = 15)
    plt.legend()
    plt.savefig(result_dir + r"/"+"%s_nonlinear_exp_model.pdf"%(batch),dpi=300,bbox_inches='tight')
    plt.savefig(result_dir + r"/"+"%s_nonlinear_exp_model.png"%(batch),dpi=300,bbox_inches='tight')
    plt.close("all")
    return popt,R_square

#########################################fit natural cubic spline model##########################
def fit_natural_cubic_spline_model(table,mean_a_list,slope_list,used_for_model_fitting=[0.05,0.95],result_dir="",batch=""):
   # nonlinear_dir = result_dir+"/nonlinear_model"
   # mkdir(nonlinear_dir)
    a_value = np.asarray(table["A-value"].sort_values(ascending=True))
    start = int(len(mean_a_list)*used_for_model_fitting[0])
    end = int(len(mean_a_list)*used_for_model_fitting[1])
    x_model = np.asarray(mean_a_list[start:end])
    y_model = np.asarray(slope_list[start:end])

    nonlinear_model = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=3)

    R_square = getIndexes(nonlinear_model.predict(x_model),y_model)[-1]

    nonlinear_model1 = natural_cubic_model.get_natural_cubic_spline_model(x_model,y_model,minval=min(x_model),maxval=max(x_model),n_knots=4)

    R_square1 = getIndexes(nonlinear_model1.predict(x_model),y_model)[-1]

    if R_square1 > R_square:
        nonlinear_model=nonlinear_model1
        R_square = R_square1
#    print (R_square,R_square1)

    plt.figure(figsize=(5,5))
    plt.scatter(mean_a_list[0:start]+mean_a_list[end:],slope_list[0:start]+slope_list[end:],
                    facecolors="none",s = 50,edgecolors = "red",label = "Not used for model fitting" )
    plt.scatter(x_model, y_model,s=50,facecolors='none', edgecolors="green",label = "Used for model fitting" )

    plt.plot(a_value,nonlinear_model.predict(a_value),c="darkorange",linewidth=3,label = "Fitted natural cubic spline")
        
    plt.xlabel("Average log$_2$-intensity",size=15)
    plt.ylabel("$Estimated\ variance\ \sigma^2$",size=15)
    plt.title("%s $R^2$=%s"%(batch,str("%.3f"%(R_square))),fontsize = 15)
    plt.legend()
    plt.savefig(result_dir + r"/"+"%s_nonlinear_natural_cubic_spline_model.pdf"%(batch),dpi=300,bbox_inches='tight')
    plt.savefig(result_dir + r"/"+"%s_nonlinear_natural_cubic_spline_model.png"%(batch),dpi=300,bbox_inches='tight')
    plt.close("all")
    return nonlinear_model,R_square


    
#########################################calculate chi^2 pvalue based on expotential_function model###################################
def calculate_exp_chi2(sorted_table,popt,result_dir="",num_of_sample=4,batch=""):
    
    sample_columns = sorted_table.columns[0:num_of_sample]
    for sample_column in sample_columns:
        sorted_table["M-value of " + sample_column] = sorted_table[sample_column]-sorted_table["A-value"]
        sorted_table["z-score of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/sorted_table["Observed Variance"]**0.5
    if len(popt)==3:
        exp_function = func_nl_lsq
    else:
        exp_function = func_nl_lsp_no_constant_term
    sorted_table["Estimated Variance"] = [exp_function(mean_value,*popt) for mean_value in sorted_table["A-value"].values]
    for sample_column in sample_columns:
        sorted_table["z-transformed of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/sorted_table["Estimated Variance"]**0.5
    
    sorted_table["Chi^2-statistic of %s"%batch] = sorted_table["Observed Variance"]*(num_of_sample-1)/sorted_table["Estimated Variance"]
    sorted_table["pvalue of %s"%batch] = [scipy.stats.chi2.sf(x,num_of_sample-1) for x in sorted_table["Chi^2-statistic of %s"%batch].values]
    sorted_by_pvalue_table = sorted_table.sort_values(by = "pvalue of %s"%batch)
    sorted_by_pvalue_table["BH-corrected P-value of %s"%batch] =[i if i<1 else 1 for i in sorted_by_pvalue_table["pvalue of %s"%batch]/sorted_by_pvalue_table["pvalue of %s"%batch].rank(axis=0)*len(sorted_by_pvalue_table)]
    sorted_by_pvalue_table.to_csv(result_dir+"/"+"%s_zmap_result_output.txt"%(batch),sep="\t")
    return sorted_by_pvalue_table

##########################################calculate chi^2 pvalue based on natural cubic spline model##################################
def calculate_natural_chi2(sorted_table,natural_cubic_spline_model,result_dir="",num_of_sample=4,batch=""):

    sample_columns = sorted_table.columns[0:num_of_sample]

    for sample_column in sample_columns:
        sorted_table["M-value of " + sample_column] = sorted_table[sample_column]-sorted_table["A-value"]
        sorted_table["z-score of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/(sorted_table["Observed Variance"]**0.5)

    sorted_table["Estimated Variance"%batch] = natural_cubic_spline_model.predict(sorted_table["A-value"].values)

    for sample_column in sample_columns:
        sorted_table["z-transformed of " + sample_column] = (sorted_table[sample_column]-sorted_table["A-value"])/(sorted_table["Estimated Variance"]**0.5)

    sorted_table["Chi^2-statistic of %"%batch] = sorted_table["Observed Variance"]*(num_of_sample-1)/sorted_table["Estimated Variance"]
    sorted_table["pvalue of %s"%batch] = [max(scipy.stats.chi2.sf(x,num_of_sample-1),1e-200) for x in sorted_table["Chi^2-statistic of %s"%batch].values]
    sorted_by_pvalue_table = sorted_table.sort_values(by = "pvalue of %s"%batch)
    sorted_by_pvalue_table["BH-corrected P-value %s"%batch] =[i if i<1 else 1 for i in sorted_by_pvalue_table["pvalue of %s"%batch]/sorted_by_pvalue_table["pvalue of %s"%batch].rank(axis=0)*len(sorted_by_pvalue_table)]
    sorted_by_pvalue_table.to_csv(result_dir+"/"+"%s_zmap_result_output.txt"%(batch),sep="\t")
    return sorted_by_pvalue_table
	

	

###################################scatterplot of theoretical quantiles###############################
def theoretical_quantiles_scatterplot(zmap_result_df,result_dir = "",batch=""):
    z_score_columns=[]
    z_transformed_columns=[]
    for column in zmap_result_df.columns:
        if "z-score" in column:
            z_score_columns.append(column)
        elif "z-trans" in column:
            z_transformed_columns.append(column)
    
    plt.figure(figsize=(8,10))
    plt.subplot(211) 
    plt.title(batch,size=20)
    for z_score_column in z_score_columns:
        z_score = zmap_result_df[z_score_column].sort_values(ascending=True)
        #z_score_quantiles_percent = [stats.percentileofscore(z_score,x) for x in z_score]
        z_score_quantiles_percent = MS_plotting_position(1,len(z_score)+1,len(z_score),0.375)
        theoretical_quantiles = [scipy.stats.norm.ppf(j) for j in z_score_quantiles_percent]
        plt.scatter(theoretical_quantiles,z_score,label=z_score_column)
    plt.xlabel("Theoretical quantiles of N(0,1)",size=15)
    plt.ylabel("Z-score of $log_2$-intensity",size=15)
    plt.plot(theoretical_quantiles,theoretical_quantiles,'k--',label="y = x")
    plt.title(batch,size=20)
    plt.legend()
    plt.subplot(212)
    for z_transformed_column in z_transformed_columns:
        z_transformed = zmap_result_df[z_transformed_column].sort_values(ascending=True)
        #z_transformed_quantiles_percent = [stats.percentileofscore(z_transformed,x) for x in z_transformed]
        
        z_transformed_quantiles_percent= MS_plotting_position(1,len(z_transformed)+1,len(z_transformed),0.375)
        
        theoretical_quantiles = [scipy.stats.norm.ppf(j) for j in z_transformed_quantiles_percent]
        plt.scatter(theoretical_quantiles,z_transformed,label=z_transformed_column)
    plt.xlabel("Theoretical quantiles of N(0,1)",size=15)
    plt.ylabel("Z-transformed of $log_2$-intensity",size=15)
    plt.plot(theoretical_quantiles,theoretical_quantiles,'k--',label="y = x")
    plt.legend()
    plt.savefig(result_dir+"/"+"%s_theoretical_quantiles.pdf"%(batch),dpi=300)
    plt.savefig(result_dir+"/"+"%s_theoretical_quantiles.png"%(batch),dpi=300)
    plt.close('all')

def combine_result(all_rep_zmap_result_xls_list,sample_info_df,ms_run_list,out_dir):
    results_df_list=[]
    batch_id_list = []
    i=0
    z_transformed_df_list=[]
    chi_pvalue_df_list = []
    for result in all_rep_zmap_result_xls_list:
        
        result_df = pd.read_csv(result,sep="\t",index_col = 0)
        batch_id = ms_run_list[i]
        i+=1
        batch_id_list.append(batch_id)
        z_trans_list = []
        for column in result_df.columns:
            if "z-transformed" in column:
                z_trans_list.append(column)
            else:
                pass
        z_df = result_df[z_trans_list]
        z_transformed_df_list.append(z_df)
        chi_pvalue_df = result_df.iloc[:,-3:]
        chi_pvalue_df_list.append(chi_pvalue_df)
	
    all_z_df = pd.concat(z_transformed_df_list,axis=1,sort=True)
    all_z_df.columns = [i.strip("z-transformed of ") for i in all_z_df.columns]
    all_z_df.to_csv(out_dir+"/"+"z_statistic_table.txt",sep="\t")
    all_chi_pvalue_df = pd.concat(chi_pvalue_df_list,axis=1,sort=True)
#    all_chi_pvalue_df.to_csv(out_dir+"/"+"chi_square_pvalue.xls",sep="\t")
    
    chi_list = []
    for column in all_chi_pvalue_df.columns:
        if "Chi" in column:
            chi_list.append(column)
    chi_df = all_chi_pvalue_df[chi_list]
        

     
    all_chi_pvalue_df["Average Chi^2-statistic"] = all_chi_pvalue_df.loc[:,chi_list].mean(axis=1)

    N = len(set(sample_info_df['Sample_condition'].values))
#    print (N)

    all_chi_pvalue_df["Number of detected"] = chi_df.count(axis=1)

    all_chi_pvalue_df["combined_pvalue"]= [scipy.stats.gamma.sf(all_chi_pvalue_df["Average Chi^2-statistic"].loc[index],
                       all_chi_pvalue_df["Number of detected"].loc[index]*(N-1)/2,
                       scale = 2/all_chi_pvalue_df["Number of detected"].loc[index]) for index in all_chi_pvalue_df.index]
    all_sort_df = all_chi_pvalue_df.sort_values(by = "combined_pvalue")
    all_sort_df["BH-corrected combined pvalue"] =[i if i<1 else 1 for i in \
               all_sort_df["combined_pvalue"]/all_sort_df["combined_pvalue"].rank(axis=0)*len(all_sort_df)]


    all_sort_df.to_csv(out_dir+"/"+"zmap_chi_square_pvalue.txt",sep="\t")

def pvalue_scatterplot(table,result_dir="",batch=""):
	
	#table = pd.read_csv(table,sep="\t",index_col=0)
	sort_table = table.sort_values(by = "Observed Variance",ascending=True)
	
	plt.figure(figsize=(6,5))
	plt.scatter(sort_table["A-value"],sort_table["Observed Variance"],s=20,c = -np.log10(sort_table["pvalue of %s"%batch].values),
                edgecolors="black",cmap="Reds",vmin=0, vmax=7,linewidth=0.2,alpha=0.8)
	plt.xlabel("Mean of $log_2-intensities$",size = 12)
	plt.ylabel("Variance of $log_2-intensities$" , size = 12)
	cbar = plt.colorbar()    
	cbar.ax.tick_params(labelsize=10)
	cbar.solids.set_edgecolor("face")    
	cbar.set_label('$-log_{10}{(P-value)}$',size=10)
	plt.title(batch,size=20)
    #all_hemoglobin_protein = ["HBE1","HBM","HBB","HBD","HBQ1","HBA1","AHSP","HBG1","HBG2","HBZ"]
    #for gene in all_hemoglobin_protein:
    #    plt.annotate(s=gene,xy=(table.loc[gene]["A-value"],table.loc[gene]["Observed Variance"]),
    #                 xytext = (table.loc[gene]["A-value"]+1,table.loc[gene]["Observed Variance"]+0.5),
    #                 arrowprops={'arrowstyle':'fancy'})
	plt.savefig(result_dir+"/"+"%s_pvalue.pdf"%(batch),dpi=300,bbox_inches="tight")
	plt.savefig(result_dir+"/"+"%s_pvalue.png"%(batch),dpi=300,bbox_inches="tight")
	plt.close('all')

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="")
    parser.add_argument("--protein_intensity_file",help="protein intensity file",required=True)
    parser.add_argument("--sample_info",help="sample_information",required=True)
    parser.add_argument("--outdir",help="outdir for result",required=True)
    parser.add_argument("--window_size",help="window size for linear regression,default is 400",required=True)
    parser.add_argument("--step_size",help="step for sliding window,default is 100",required=True)
    parser.add_argument("--percent",help="left percent proteins used to linear fitting, 20<=percent<=50",required=True)
    parser.add_argument("--method",help="method used to fit nonlinear model,the method can be exponential_function,natural_cubic_spline",required=True)
    argv = vars(parser.parse_args())
    protein_intensity_file = argv['protein_intensity_file']
    sample_info = argv['sample_info']
    outdir = argv['outdir']
    method = argv['method']
    window_size=int(argv['window_size'])
    step_size=int(argv['step_size'])
    percent = int(argv['percent'])/100
    
    mkdir(outdir)
    raw_intensity_data = pd.read_csv(protein_intensity_file,sep="\t",index_col=0)
#    sheet_names= raw_intensity_data.sheet_names
    sample_info_df = pd.read_csv(sample_info,sep="\t",index_col=0)
    ms_run_list = list(set(list(sample_info_df['MS_run'].values)))
    
    all_rep_zmap_result_xls_list =[]
    
    for ms_run in ms_run_list:
        batch = ms_run
        output_dir = outdir +"/"+batch
        mkdir(output_dir)
        all_rep_zmap_result_xls_list.append(output_dir+"/"+"%s_zmap_result_output.txt"%(batch))
		# """
        sample_ = list(sample_info_df.loc[sample_info_df['MS_run']==batch].index)
        raw_intensity_df = raw_intensity_data[sample_]
        raw_intensity_df = raw_intensity_df.dropna(how="any")
        raw_intensity_df = remove_bad_gene(raw_intensity_df,cutoff=10)
        normalized_df = normalization(raw_intensity_df,output_dir,l=1.5,batch=batch)

        table_df,sample_num = calculate_a_variance_value(normalized_df)

        sorted_table = get_all_sliding_region(table_df)

	
        mean_a_list,slope_list = window_vaplot(sorted_table,sample_num,window_size = window_size,
                                           sliding_size = step_size,result_dir = output_dir,
                                           top_gene_used_to_fit_linear_model = percent,batch=batch)

        if method == "exponential_function":
            popt, exp_R_square = fit_exp_model(sorted_table,mean_a_list,slope_list,
                         used_for_model_fitting=[0.05,0.95],result_dir=output_dir,batch=batch)
            zmap_result_df = calculate_exp_chi2(sorted_table,popt,result_dir = output_dir,num_of_sample=sample_num,batch=batch)
        else:# method == "natural_cubic_spline":
            natural_cubic_spline_model,ncs_R_square = fit_natural_cubic_spline_model(sorted_table,mean_a_list,slope_list,used_for_model_fitting=[0.05,0.95],result_dir=output_dir,batch=batch)
            zmap_result_df = calculate_natural_chi2(sorted_table,natural_cubic_spline_model,result_dir = output_dir,num_of_sample = sample_num,batch=batch)

        pvalue_scatterplot(zmap_result_df,result_dir=output_dir,batch=batch)
        theoretical_quantiles_scatterplot(zmap_result_df,result_dir = output_dir,batch=batch)
		
     
    combine_result(all_rep_zmap_result_xls_list,sample_info_df,ms_run_list,out_dir=outdir)
    
		
        

    
    
    
    
    
    
