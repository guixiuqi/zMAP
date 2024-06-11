#export scriptdir

scriptPATH=/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/python_script
export scriptPATH
gmtPATH=/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data/gmt
inputdataPATH=/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data

#zMAP
#python $scriptPATH/zmap_step1_quantity_anaysis.py --protein_intensity_file $inputdataPATH/test_data/raw_protein_intensity_in_gene_level_for_web.txt --sample_info $inputdataPATH/test_data/zmap_sample_info_for_web.txt --outdir $inputdataPATH/test_data/zMAP_results --window_size 400 --step_size 100 --percent 30 --method exponential_function

######## downstream analysis ###########

#Sample quanlity control
#python $scriptPATH/sample_quality_control.py --z_statistic_matrix $inputdataPATH/test_data/zMAP_results/z_statistic_table.txt --sample_info $inputdataPATH/test_data/zmap_sample_info_for_web.txt --outdir $inputdataPATH/test_data/sample_quality_control



#HVPs identification and clustering
#python $scriptPATH/zmap_hypervariable_proteins_cluster.py --pvalue_results $inputdataPATH/test_data/zMAP_results/zmap_chi_square_pvalue.txt --z_statistic_matrix $inputdataPATH/test_data/zMAP_results/z_statistic_table.txt --sample_info $inputdataPATH/test_data/zmap_sample_info_for_web.txt --outdir $inputdataPATH/test_data/zmap_hypervariable_proteins_cluster --cluster_number_for_hypervariable 15 --minclustersize 20 --top 100 --cluster_number_for_top_proteins 8 --fdr 0.05


#Gene set variation analysis
python $scriptPATH/gsva.py --z_statistic_matrix $inputdataPATH/test_data/zMAP_results/z_statistic_table.txt --sample_info $inputdataPATH/test_data/gsva_sample_info.txt --outdir $inputdataPATH//test_data/gsva --top_n 50 --fdr 0.05
