# zMAP toolset
Isobaric labeling-based mass spectrometry (ILMS) has been widely used to quantify, on a proteome-wide scale, the relative protein abundance in different biological conditions. 
However, large-scale ILMS data sets typically involve multiple runs of mass spectrometry (MS), bringing great computational difficulty to the integration of ILMS samples. 
We present zMAP toolset that makes ILMS intensities comparable across MS runs by modeling the associated mean-variance dependence and accordingly applying a variance stabilizing z-transformation.

This application provides two computation modules (zMAP and reverse-zMAP) for transforming ILMS intensities into z-statistics and thus integrating samples across MS runs,
as well as a series of interfaces to common downstream analyses on the integrated z-statistic matrix.
![zMAP toolset](https://github.com/guixiuqi/zMAP/blob/main/imgs/zMAP_overview.png "zMAP toolset")

zMAP is provided at https://github.com/guixiuqi/zMAP/

reverse-zMAP is provided at https://github.com/guixiuqi/reverse-zMAP

Contact：guixiuqi@gmail.com
# zMAP

The design of zMAP aims to simultaneously compare protein profiles of multiple samples and integrate samples from different MS runs for identifying hypervariable 
proteins across samples. zMAP models the mean-variance dependence associated with ILMS intensities and accordingly applies a variance-stabilizing z-transformation
(z-statistic), which dramatically increases the comparability between samples generated by different MS runs. The z-statistic can be widely applied in downstream 
analyses, including PCA, clustering analysis, GSVA, and so on.
zMAP requires that all the involved MS runs are associated with the same biological design, such that the average intensities across all samples in each run are biologically identical.

# Workflow of zMAP

![Workflow of zMAP](https://github.com/guixiuqi/zMAP/blob/main/imgs/zMAP_workflow.png "zMAP Workflow")

# Try it

A Web-based application of zMAP is provided at http://bioinfo.cemps.ac.cn/shaolab/zMAP. 


# Installation
## From source code
The scripts of zMAP  require no installation and can be used in-place. Just install the dependencies (see below)

```python
git clone https://github.com/guixiuqi/zMAP.git
cd ./zmap
scriptPATH=./python_script
export scriptPATH
gmtPATH=./data/gmt
inputdataPATH=./data/test_data
```
## Python dependencies
```bat
conda install matplotlib  #version 3.6.2
conda install conda-forge::pandas #version 1.5.1
conda install anaconda::scipy #version 1.10.1
conda install anaconda::seaborn #version 0.12.2
conda install -c conda-forge scikit-learn #version 1.1.3
```

## R dependencies
```bat
library(GSVA)  #version 1.28.0
library(limma) #version 3.36.5
library(GSEABase,quietly=TRUE) #version 1.42.0

```


# Usage

## zMAP
zMAP models the mean-variance dependence associated with ILMS intensities and accordingly applies a variance-stabilizing z-transformation (z-statistic), which dramatically increases the comparability between samples generated by different MS runs. The z-statistic can be widely applied in downstream analyses.

Two input files provided by user:

1. Protein intensity file
   
    A tab-delimited file containing raw gene-level protein intensity with samples in columns, and gene symbols in rows.
Note:(1). The protein intensity matrix does not require normalization.(2). Sample names can only consist of letters, numbers, and underscores.

2. Sample information file

    Sample information file is a three-column, tab-delimited file with the first line identifying the columns. The column names are ```Sample_id```, ```Sample_condition```, and ```MS_run```.

Use the command shown as below:

```bat
python $scriptPATH/zmap.py --protein_intensity_file $inputdataPATH/raw_protein_intensity_in_gene_level_for_web.txt --sample_info $inputdataPATH/zmap_sample_info_for_web.txt --outdir $inputdataPATH/zMAP_results --window_size 400 --step_size 100 --percent 30 --method exponential_function
```
## Downstream analyses based on z-statistic

### Sample Quality Control

Performing principal component and hierarchical clustering analysis on the z-statistic matrix of samples provides a concise visualization of the overall impact of sample conditions and MS runs.
![QC of zMAP](https://github.com/guixiuqi/zMAP/blob/main/imgs/zmap_QC.png "zMAPQC")

Two input files provided by user:

1. Output file z_statistic_table.txt from zMAP


    Note:(1). The protein intensity matrix does not require normalization.(2). Sample names can only consist of letters, numbers, and underscores.

2. Sample information file

    Sample information file is a three-column, tab-delimited file with the first line identifying the columns. The column names are ```Sample_id```, ```Sample_condition```, and ```MS_run```.

Use the command shown as below:
```bat
python $scriptPATH/sample_quality_control.py --z_statistic_matrix $inputdataPATH/zMAP_results/z_statistic_table.txt --sample_info $inputdataPATH/zmap_sample_info_for_web.txt --outdir $inputdataPATH/sample_quality_control
```
### Hypervariable proteins(HVPs) identification and clustering

The chi-square statistics derived by zMAP along with the associated numbers of degrees of freedom were summed across MS runs for each protein, giving rise to a p-value that assessed the overall expression variability of the protein. This functionality is designed toOut select hypervariable proteins (HVPs) from the output files of zMAP and cluster them to obtain multiple expression signatures. Then, pathway enrichment analysis was conducted on these signatures to reveal biological insights.
![HVPs](https://github.com/guixiuqi/zMAP/blob/main/imgs/zmap_hvps_cluster.png "zMAP_Hvps")

Three input files provided by user:

1. Output file zmap_chi_square_pvalue.txt from zMAP 


2. Output file z_statistic_table.txt from zMAP


    Note: (1). The protein intensity matrix does not require normalization. (2). Sample names can only consist of letters, numbers, and underscores.

3. Sample information file

    Sample information file is a three-column, tab-delimited file with the first line identifying the columns. The column names are ```Sample_id```, ```Sample_condition```, and ```MS_run```.
   
Parameters:

1. Number of HVP clusters (```--cluster_number_for_hypervariable```)

   HVPs are hierarchically clustered into multiple clusters revealing diverse expression signatures across different sample conditions.

2. Minimum proteins within each cluster (```--minclustersize```)

   If the number of proteins within certain clusters falls below the specified minimum, these clusters will be merged into a single cluster labeled as cluster_0. As a result, the final number of clusters may be fewer than what the user initially specified.

3. Number of top-ranked DEPs (```--top```)

   Top-ranked proteins are selected for heatmap visualization.

4. Number of DEP clusters (```--cluster_number_for_top_proteins```)

   Top-ranked DEPs are clustered into multiple clusters based on K-means.

5. FDR (```--fdr```)

   BH-adjusted pvalue cutoff for significance.

Use the command shown as below:
```bat
python $scriptPATH/zmap_hypervariable_proteins_cluster.py --pvalue_results $inputdataPATH/zMAP_results/zmap_chi_square_pvalue.txt --z_statistic_matrix $inputdataPATH/zMAP_results/z_statistic_table.txt --sample_info $inputdataPATH/zmap_sample_info_for_web.txt --outdir $inputdataPATH/zmap_hypervariable_proteins_cluster --cluster_number_for_hypervariable 15 --minclustersize 20 --top 100 --cluster_number_for_top_proteins 8 --fdr 0.05
```
### Gene set variation analysis
GSVA calculates gene set enrichment scores (GSVA scores) for each sample using the z-statistic matrix.
Differential expression analysis is then conducted on these GSVA scores using limma, aiming to identify differentially regulated pathways across sample groups. Finally, the differential pathway activities across sample groups are visualized using a heatmap.
![GSVA of zMAP](https://github.com/guixiuqi/zMAP/blob/main/imgs/zmap_gsva.png "zMAP_GSVA")

Two input files provided by user:

1. Output file z_statistic_table.txt from zMAP


    Note:(1). The protein intensity matrix does not require normalization.(2). Sample names can only consist of letters, numbers, and underscores.

2. Sample information file

    Sample information file is a three-column, tab-delimited file with the first line identifying the columns. The column names are ```Sample_id```, ```Sample_condition```, and ```MS_run```.

Parameters:

1. Top N(number) Pathways (```--top_n```)
 
   Extract a table of the top-ranked differentially regulated pathways across sample groups.

2. FDR (```--fdr```)

   BH-adjusted pvalue cutoff for significance.



Use the command shown as below:
```bat
python $scriptPATH/gsva.py --z_statistic_matrix $inputdataPATH/zMAP_results/z_statistic_table.txt --sample_info $inputdataPATH/gsva_sample_info.txt --outdir $inputdataPATH/gsva --top_n 50 --fdr 0.05
```


# Citation

To cite the zMAP toolset in publications, please cite

Gui, X., Huang, J., Ruan, L. et al. zMAP toolset: model-based analysis of large-scale proteomic data via a variance stabilizing z-transformation. Genome Biol 25, 267 (2024). https://doi.org/10.1186/s13059-024-03382-9

# License
zMAP is licensed under the GNU General Public License v3.0 (GPLv3) and is available for free use in academic and research settings, provided that all terms of the GPLv3 are adhered to, including the requirement that any modifications or derivative works must also be distributed under the same license.


