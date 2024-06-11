library(GSVA)
library(limma)
library(GSEABase,quietly=TRUE)
group_info <- read.csv("/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data/test_data/gsva_sample_info.txt",sep="\t",row.names = 1,stringsAsFactors = FALSE)
data <- read.csv("/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data/test_data/zMAP_results/z_statistic_table.txt",sep="\t",row.names = 1)
data <- na.omit(as.matrix(data))
keggset <- getGmt("/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data/gmt/human_GO_Biological_Process_2015.gmt")
keggs_kcdf_none <- gsva(data,gset.idx.list=keggset,kcdf="Gaussian", parallel.sz=1,min.sz=2)
write.csv(keggs_kcdf_none,file="/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data//test_data/gsva/go_gsva_kcdf_Gaussian.csv")
design <- factor(c(group_info$Sample_condition))
design1 <- model.matrix(~0+design)
colnames(design1) <- levels(design)
fit <- lmFit(keggs_kcdf_none, design1)
contrast.matrix <- makeContrasts(ConditionB-ConditionD,ConditionB-ConditionC,ConditionB-ConditionA,ConditionD-ConditionC,ConditionD-ConditionA,ConditionC-ConditionA, levels = design1)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, number= 50)
write.csv(fit2,file="/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data//test_data/gsva/go_gsva_differential_pathway_activity.csv")
write.csv(top,file="/mnt/storage/user/guixiuqi/rsgeno/protein/zmap_github/data//test_data/gsva/go_gsva_differential_pathway_activity_top_50.csv")