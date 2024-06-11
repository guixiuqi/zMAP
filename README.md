# zMAP

The design of zMAP aims to simultaneously compare protein profiles of multiple samples and integrate samples from different MS runs for identifying hypervariable 
proteins across samples. zMAP models the mean-variance dependence associated with ILMS intensities and accordingly applies a variance-stabilizing z-transformation
(z-statistic), which dramatically increases the comparability between samples generated by different MS runs. The z-statistic can be widely applied in downstream 
analyses, including PCA, clustering analysis, GSVA, and so on.
zMAP requires that all the involved MS runs are associated with the same biological design, such that the average intensities across all samples in each run are biologically identical.

# Workflow of zMAP

![Workflow of zMAP](https://github.com/guixiuqi/zMAP/blob/main/imgs/zMAP_workflow.png "zMAP Workflow")

# Try it

A Web-based application of zMAP is provided at http://bioinfo.sibs.ac.cn/shaolab/zMAP. 


# Installation
## From source code
The scripts of zMAP  require no installation and can be used in-place. Just install the dependencies (see below)

```python
git clone https://github.com/guixiuqi/zMAP.git
cd ./zmap
scriptPATH=./python_script
export scriptPATH
gmtPATH=./data/gmt
inputdataPATH=./data
```
## Dependencies
```bat
conda install matplotlib
conda install conda-forge::pandas
conda install anaconda::scipy
conda install anaconda::seaborn
conda install -c conda-forge scikit-learn
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
python $scriptPATH/zmap_step1_quantity_anaysis.py --protein_intensity_file $inputdataPATH/test_data/raw_protein_intensity_in_gene
_level_for_web.txt --sample_info $inputdataPATH/test_data/zmap_sample_info_for_web.txt --outdir $inputdataPATH/test_data/zMAP_res
ults --window_size 400 --step_size 100 --percent 30 --method exponential_function
```
## Downstream analyses based on z-statistic

### Sample Quality Control

Performing hierarchical clustering and principal component analysis on the z-statistic matrix of samples provides a concise visualization of the overall impact of sample conditions and MS runs.

Two input files provided by user:

1. Output file z_statistic_table.txt from zMAP


    Note:(1). The protein intensity matrix does not require normalization.(2). Sample names can only consist of letters, numbers, and underscores.

2. Sample information file

    Sample information file is a three-column, tab-delimited file with the first line identifying the columns. The column names are ```Sample_id```, ```Sample_condition```, and ```MS_run```.

Use the command shown as below:
```bat
python $scriptPATH/zmap_step1_quantity_anaysis.py --protein_intensity_file $inputdataPATH/test_data/raw_protein_intensity_in_gene
_level_for_web.txt --sample_info $inputdataPATH/test_data/zmap_sample_info_for_web.txt --outdir $inputdataPATH/test_data/zMAP_res
ults --window_size 400 --step_size 100 --percent 30 --method exponential_function
```




```r
# Get a motif pwms
example_motifs <- getJasparMotifs(species = "Homo sapiens",
                                  collection = "CORE")
# Make a set of input regions
Input <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
                                ranges = IRanges::IRanges(start = c(76585873,42772928,100183786),
                                                          width = 500))
# Make a set of control regions
Control <- GenomicRanges::GRanges(seqnames = c("chr1","chr3","chr5"),
                                  ranges = IRanges::IRanges(start = c(453123,6524593,100184233),
                                                            width = 500))
# Scan motif for example motifs
motif_ix_input <- motifScan(example_motifs, Input, genome = "BSgenome.Hsapiens.UCSC.hg19")
motif_ix_control <- motifScan(example_motifs, Control, genome = "BSgenome.Hsapiens.UCSC.hg19")

# Find Enrichment motif of input by control
Enrichment_result <- motifEnrichment(motif_ix_input, motif_ix_control)
```
You could type `?motifEnrichment` for more detail.

# Citation

To cite the `zMAP` package in publications, please use


