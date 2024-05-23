# DESeq Analysis for 358 Cells
# Define the project name, organize folder structure. Script in project root.
project <- "RNASeqAnalysis"
# Get the Rstudio Directory
rstudio_dir <- rstudioapi::getActiveProject()
# Define the project directory inside Rstudio directory
project_dir <- file.path(rstudio_dir, project)
if (!file.exists(project_dir)) {dir.create(project_dir, recursive = TRUE)}
# Define input and output directories
input_dir <- file.path(project_dir, "input")
output_dir <- file.path(project_dir, "output")
# Check if input and output directories exist, if not, create them
if (!file.exists(input_dir)) {dir.create(input_dir, recursive = TRUE)}
if (!file.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
# Set the working directory to the project directory
setwd(project_dir)
getwd()
#Print Confirmation of Correct Folder Structure
if(basename(getwd()) == project) "Folder setup correctly" else "Fix folder structure"
#-----------
# Project specific output for V0.4 UnTrimmedNewIndexHierarchicalCategory
output_dir <- "/Users/i/Dropbox/Clinic3.0/RStudio/RNASeqAnalysis/output/v0.41_PerCellLine"
#-----------

#-----------
#Define functions
saveFigure <- function(figure, fileName, h = 7, w = 7, dpi = 300) {
  currentDate <- format(Sys.Date(), "%Y%m%d") #current date in YYYYMMDD format
  # Define the directory for saving figures
  figuresDir <- file.path(output_dir, "figures")
  if (!dir.exists(figuresDir)) { dir.create(figuresDir, recursive = TRUE) }
  fullFilePath <- file.path(figuresDir, paste0(currentDate, "_", fileName, ".pdf"))
  # Save the figure
  pdf(file = fullFilePath, height = h, width = w, pointsize = dpi / 72)
  print(figure)
  dev.off()
}
############################
#BiocManager::install("genefilter")
#library(apeglm)
sapply(c("tidyverse", 
         "devtools", 
         "annotate", 
         "org.Hs.eg.db", 
         "DESeq2",
         "apeglm"), 
       require, character.only = TRUE)
#Get Inputs
colData_358_file<-paste0(input_dir,"/samplesheet_358.csv")
colData_358<-read.csv(file = colData_358_file, 
                      header=TRUE, 
                      stringsAsFactors = FALSE, 
                      check.names = FALSE,
                      row.names = "sample")
colData_358 <- colData_358[,c("cellLine", "drug")]
colData_358$cellLine <- factor(colData_358$cellLine)
colData_358$drug <- factor(colData_358$drug)

fc_file<-paste0(input_dir,"/featureCounts_0")
fc<-read.delim(fc_file,row.names=1,check.names = FALSE)
#Keep rows where more than one genes are expressed
genes_to_keep <- rowSums(fc[,7:ncol(fc)]) > 1
fc_3<-fc[genes_to_keep,]
fc_3<-fc_3[,6:ncol(fc_3)] #Remove the initial columns
fc_3$EnsembleID<-gsub("\\..*$","",rownames(fc_3)) #New column with removed suffix of Ensemble ID with the dot
fc_3<-fc_3[order(fc_3[,1],decreasing = TRUE),] # Order the rows based on gene length
fc_3<-fc_3[!duplicated(fc_3$EnsembleID),] #Remove duplicated rows based on EnsembleID keeping the first occurrence

rownames(fc_3) <- fc_3$EnsembleID #Rownames as the dot removed ensembleID
fc_3 <- subset(fc_3, select= -Length)
fc_3 <- subset(fc_3, select = -EnsembleID)

fc_358 <- fc_3[, grepl("^358", names(fc_3))]


desired_order_358 <- rownames(colData_358)
fc_358 <- fc_358[, desired_order_358]
all(colnames(fc_358) %in% rownames(colData_358))
all(colnames(fc_358) == rownames(colData_358))

dds_358 <- DESeqDataSetFromMatrix(countData = fc_358,
                                  colData = colData_358,
                                  design = ~ drug)


keep_358 <-  rowSums(counts(dds_358))>= 10 #Keep genes which have at least 10 detected across samples
dds_358 <- dds_358[keep_358,]
#dds_358$drug
dds_358$drug <- relevel(dds_358$drug, ref="Control")

dds358 <- DESeq(dds_358)
res358 <- results(dds358)
resultsNames(dds358)

plotMA(res358, ylim=c(-2, 2))


#Get Colors########
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#FEAEAF", "#FF7B2A"))(99)


sapply(c("genefilter", "clusterProfiler", "EnhancedVolcano", "GSVA", "DOSE", "RColorBrewer", "pheatmap", "xCell"), require, character.only = TRUE)
rld358<-vst(dds358) #estimate dispersion trend and apply a variance stabilizing transformationrld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformation

## creating distance matrix
sampleDists_subset358 <- as.matrix(dist(t(assay(rld358))))
hm358<-pheatmap::pheatmap(as.matrix(sampleDists_subset358),annotation_col = colData_358, col=my_colors,annotation_legend=TRUE)

## creating PCA plot
pca<-plotPCA(rld358, intgroup="drug")
# pca<-pca + geom_text(aes(label=name),vjust=2, size = 3)
# saveFigure(figure=pca,fileName="PCAPlot",h=10,w=20)
#If needs to label samples using ggrepel
library(ggrepel)
zz358 <- plotPCA(rld358, intgroup="drug", returnData=TRUE)
pca_repel358 <- pca+geom_label_repel(data=zz358, aes(label=name))+
  ggtitle(label="PCA plot of 358 Cells")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

saveFigure(figure=pca_repel358, fileName = "PCAPlot_Repel", h=15, w=20)
colData(dds358)


##### Plotting DEG
ComparisonColumn <- "drug"
factor1 <- "Control"
factor2 <- "128-10"
resultsNames(dds358)
e <- as.character(c(ComparisonColumn, factor1, factor2))
res358 <- lfcShrink(dds358, contrast = e, type = "normal")

resdata_subset358 <- merge(as.data.frame(res358), as.data.frame(assay(rld358)), by="row.names", sort=FALSE)
write.csv(resdata_subset358, paste0(output_dir,"/","DifferentialExpressionAnalysis358", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset358)[1] <- "Gene"
resdata_subset358 <- resdata_subset358[order(resdata_subset358$padj),]
resdata_subset358 <- resdata_subset358[!is.na(resdata_subset358$padj),]

#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db


resdata_subset358$GeneSymbol<- mapIds(sp, keys=resdata_subset358$Gene, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset358$EntrezID<- mapIds(sp, keys=resdata_subset358$Gene, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
write.csv(resdata_subset358, paste0(output_dir,"/","DifferentialExpressionAnalysis_cleaned358", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

resdatagenes_subset358 <- resdata_subset358[complete.cases(resdata_subset358$GeneSymbol),]
resdatagenes_subset358 <- resdatagenes_subset358[!duplicated(resdatagenes_subset358$GeneSymbol), ]

resdatagenes_subset358 <- resdatagenes_subset358[order(resdatagenes_subset358$log2FoldChange),]
Top50genesdown_subset358<-head(resdatagenes_subset358, 50)
Top50genesup_subset358<-tail(resdatagenes_subset358, 50)
Top_subset358<-rbind(Top50genesdown_subset358, Top50genesup_subset358)


R358 <- as.character(rownames(colData_358))
R358<-c(R358, "GeneSymbol")
Top_subset358<-Top_subset358%>% dplyr::select(all_of(R358))
TopD_subset358<-data.frame(Top_subset358, check.names = FALSE)
#name the rows as genesymbols
rownames(TopD_subset358)<-Top_subset358$GeneSymbol

## top 50 fold changes
hmt50358<-pheatmap::pheatmap(TopD_subset358[1:(length(TopD_subset358)-1)],scale="row", annotation_col=colData_358,
                          annotation_legend =TRUE, color= colorRampPalette(c("blue","white","red"))(99), fontsize_row = 8, main="Top50 FoldChange in 358 Cell Line")
saveFigure(figure=hmt50358,fileName="Top50FoldChange_heatmap",h=12,w=12)


## variable genes
Ds358<-resdatagenes_subset358%>% dplyr::select(all_of(R358), "GeneSymbol")
Ds358<-Ds358[!duplicated(Ds358$GeneSymbol), ]
rownames(Ds358)<-Ds358$GeneSymbol
Ds358<-Ds358%>% dplyr::select(-"GeneSymbol")
#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds358)),100)
mat358 <- Ds358[topVarGenes, ]
mat358 <- mat358 - rowMeans(mat358)

#plot the variable genes in heatmap
hmVariable358<-pheatmap::pheatmap(mat358,
                annotation_col=colData_358,
                color= colorRampPalette(c("blue","white","red"))(99), 
                annotation_colors = list(drug=c("128-10"="#F0978D", 
                                                "128-13"="#63D7DE",
                                                "130"="#007DEE",
                                                "Control"="#999000")),
                annotation_legend =TRUE, 
                scale="row", 
                fontsize_row = 8, 
                show_rownames=T, 
                main="Top100 Variable Genes in 358 Cell Line")
saveFigure(figure=hmVariable,fileName="Top100VariableGenes",h=12,w=12)

# Volcano Plot
vp358<-EnhancedVolcano(resdata_subset358,
                    lab = resdata_subset358$GeneSymbol,x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 10e-4,
                    FCcutoff = 1)
