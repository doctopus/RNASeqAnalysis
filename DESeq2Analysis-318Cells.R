#DESeq2 Analsysis In Individual Cell Lines
#Based on DESeq2 Vignette

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
BiocManager::install("apeglm")
library(apeglm)
sapply(c("tidyverse", 
         "devtools", 
         "annotate", 
         "org.Hs.eg.db", 
         "DESeq2",
         "apeglm"), 
       require, character.only = TRUE)
#Get Inputs
colData_318_file<-paste0(input_dir,"/samplesheet_318.csv")
colData_318<-read.csv(file = colData_318_file, 
                  header=TRUE, 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE,
                  row.names = "sample")
colData_318 <- colData_318[,c("cellLine", "drug")]

colData_318_file<-paste0(input_dir,"/samplesheet_318.csv")
colData_318<-read.csv(file = colData_318_file, 
                  header=TRUE, 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE,
                  row.names = "sample")
colData_318 <- colData_318[,c("cellLine", "drug")]
#colData <- colData %>% dplyr::select(-fastq_1, -fastq_2, -strandedness, -replicate, -group)
#colData <- colData %>% column_to_rownames(var="sample")

colData_318$cellLine <- factor(colData_318$cellLine)
colData_318$drug <- factor(colData_318$drug)

colData_318$cellLine <- factor(colData_318$cellLine)
colData_318$drug <- factor(colData_318$drug)
#head(colData)

fc_file<-paste0(input_dir,"/featureCounts_0")
fc<-read.delim(fc_file,row.names=1,check.names = FALSE)
#Keep rows where more than one genes are expressed
genes_to_keep <- rowSums(fc[,7:ncol(fc)]) > 1
fc_3<-fc[genes_to_keep,]
fc_3<-fc_3[,6:ncol(fc_3)] #Remove the initial columns
fc_3$EnsembleID<-gsub("\\..*$","",rownames(fc_3)) #New column with removed suffix of Ensemble ID with the dot
fc_3<-fc_3[order(fc_3[,2],decreasing = TRUE),] # Order the rows
fc_3<-fc_3[!duplicated(fc_3$EnsembleID),] #Remove duplicated rows based on EnsembleID keeping the first occurrence

rownames(fc_3) <- fc_3$EnsembleID #Rownames as the dot removed ensembleID
fc_318 <- subset(fc_3, select= -Length)
fc_318 <- subset(fc_318, select = -EnsembleID)

fc_318 <- fc_318[, grepl("^318", names(fc_318))]
#Get Bioconductor Annotation Database
#sp <- org.Hs.eg.db

#fc_3_new <- subset(fc_3, select = -Length) #Remove the length column
#fc_3_new$GeneSymbol <- mapIds(sp, keys=fc_3_new$EnsembleID, column = c("SYMBOL"), keytype = c("ENSEMBL"), multiVals = "first")

#fc_3_new <- fc_3_new %>% dplyr::select(-tail(names(.), 1)) #Remove the EnsembleID column in the end

desired_order_318 <- rownames(colData_318)
fc_318 <- fc_318[, desired_order_318]

all(colnames(fc_318) %in% rownames(colData_318))
all(colnames(fc_318) == rownames(colData_318))

dds_318 <- DESeqDataSetFromMatrix(countData = fc_318,
                                  colData = colData_318,
                                  design = ~ drug)

keep_318 <-  rowSums(counts(dds_318))>= 10
dds_318 <- dds_318[keep,]
#dds_318$drug
dds_318$drug <- relevel(dds_318$drug, ref="Control")

dds318 <- DESeq(dds_318)
res318 <- results(dds318)
resultsNames(dds318)

plotMA(res318, ylim=c(-2, 2))

#Get Colors########
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#FEAEAF", "#FF7B2A"))(99)


sapply(c("genefilter", "clusterProfiler", "EnhancedVolcano", "GSVA", "DOSE", "RColorBrewer", "pheatmap", "xCell"), require, character.only = TRUE)
rld318<-vst(dds318) #estimate dispersion trend and apply a variance stabilizing transformationrld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformation

## creating distance matrix
sampleDists_subset318 <- as.matrix(dist(t(assay(rld318))))
hm318<-pheatmap::pheatmap(as.matrix(sampleDists_subset318),annotation_col = colData_318, col=my_colors,annotation_legend=TRUE)

## creating PCA plot
pca<-plotPCA(rld318, intgroup="drug")
# pca<-pca + geom_text(aes(label=name),vjust=2, size = 3)
# saveFigure(figure=pca,fileName="PCAPlot",h=10,w=20)
#If needs to label samples using ggrepel
library(ggrepel)
zz318 <- plotPCA(rld318, intgroup="drug", returnData=TRUE)
pca_repel318 <- pca+geom_label_repel(data=zz318, aes(label=name))+
  ggtitle(label="PCA plot of 318 Cells")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

saveFigure(figure=pca_repel318, fileName = "PCAPlot_Repel", h=15, w=20)
colData(dds318)


##### Plotting DEG
ComparisonColumn <- "drug"
factor1 <- "Control"
factor2 <- "128-10"
resultsNames(dds318)
e <- as.character(c(ComparisonColumn, factor1, factor2))
res318 <- lfcShrink(dds318, contrast = e, type = "normal")

resdata_subset318 <- merge(as.data.frame(res318), as.data.frame(assay(rld318)), by="row.names", sort=FALSE)
write.csv(resdata_subset318, paste0(output_dir,"/","DifferentialExpressionAnalysis318", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset318)[1] <- "Gene"
resdata_subset318 <- resdata_subset318[order(resdata_subset318$padj),]
resdata_subset318 <- resdata_subset318[!is.na(resdata_subset318$padj),]

#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db


resdata_subset318$GeneSymbol<- mapIds(sp, keys=resdata_subset318$Gene, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset318$EntrezID<- mapIds(sp, keys=resdata_subset318$Gene, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
write.csv(resdata_subset318, paste0(output_dir,"/","DifferentialExpressionAnalysis_cleaned318", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

resdatagenes_subset318 <- resdata_subset318[complete.cases(resdata_subset318$GeneSymbol),]
resdatagenes_subset318 <- resdatagenes_subset318[!duplicated(resdatagenes_subset318$GeneSymbol), ]

resdatagenes_subset318 <- resdatagenes_subset318[order(resdatagenes_subset318$log2FoldChange),]
Top50genesdown_subset318<-head(resdatagenes_subset318, 50)
Top50genesup_subset318<-tail(resdatagenes_subset318, 50)
Top_subset318<-rbind(Top50genesdown_subset318, Top50genesup_subset318)


R318 <- as.character(rownames(colData_318))
R318<-c(R318, "GeneSymbol")
Top_subset318<-Top_subset318%>% dplyr::select(all_of(R318))
TopD_subset318<-data.frame(Top_subset318, check.names = FALSE)
#name the rows as genesymbols
rownames(TopD_subset318)<-Top_subset318$GeneSymbol

#Keep only the compared columns
# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, sep = "|"), names(TopD_subset), value = TRUE)

# Subset the data frame to keep only those columns
Top50_f1vsf2 <- TopD_subset[, columns_to_keep]


## top 50 fold changes
hmt50318<-pheatmap::pheatmap(Top50_f1vsf2, scale="row", annotation_col=colData_318,
                             annotation_legend =TRUE, color= colorRampPalette(c("blue","white","red"))(99), fontsize_row = 8, main="Top50 FoldChange in 318 Cell Line")
saveFigure(figure=hmt50318,fileName="Top50FoldChange_heatmap",h=12,w=12)


## variable genes
Ds318<-resdatagenes_subset318%>% dplyr::select(all_of(R318), "GeneSymbol")
Ds318<-Ds318[!duplicated(Ds318$GeneSymbol), ]
rownames(Ds318)<-Ds318$GeneSymbol
Ds318<-Ds318%>% dplyr::select(-"GeneSymbol")
#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds318)),100)
mat318 <- Ds318[topVarGenes, ]
mat318 <- mat318 - rowMeans(mat318)

#plot the variable genes in heatmap
hmVariable318<-pheatmap::pheatmap(mat318,
                                  annotation_col=colData_318,
                                  color= colorRampPalette(c("blue","white","red"))(99), 
                                  annotation_colors = list(drug=c("128-10"="#F0978D", 
                                                                  "128-13"="#63D7DE",
                                                                  "130"="#007DEE",
                                                                  "Control"="#999000")),
                                  annotation_legend =TRUE, 
                                  scale="row", 
                                  fontsize_row = 8, 
                                  show_rownames=T, 
                                  main="Top100 Variable Genes in 318 Cell Line")
saveFigure(figure=hmVariable,fileName="Top100VariableGenes",h=12,w=12)

# Volcano Plot
vp318<-EnhancedVolcano(resdata_subset318,
                       lab = resdata_subset318$GeneSymbol,x = 'log2FoldChange',
                       y = 'pvalue',
                       pCutoff = 10e-4,
                       FCcutoff = 1)
