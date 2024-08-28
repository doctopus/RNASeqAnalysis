# DESeq2 Analysis of 358 Cell Lines
# Define the project name, organize folder structure. Script in project root.
project <- "RNASeq-Kinase"
# Get the Rstudio Directory
rstudio_dir <- rstudioapi::getActiveProject()
rstudio_base <- basename(rstudio_dir) # Get the base name of the Rstudio Directory
if (rstudio_base != project) { # If the base name != project, create a folder named as the project inside rstudio_dir
  project_dir <- file.path(rstudio_dir, project)
    if (!file.exists(project_dir)) {dir.create(project_dir, recursive = TRUE)}
  } else {
  # If the base name is the same as the project, use the rstudio_dir as the project_dir
  project_dir <- rstudio_dir
}
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
# Project specific output directory
output_dir <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/RNASeqAnalysis/output/v0.41_PerCellLine"
#-----------
#Define functions
savePNG <- function(figure, fileName, w = 1200, h = 900) {
  currentDate <- format(Sys.Date(), "%Y%m%d") # current date in YYYYMMDD format
  # Define the directory for saving figures
  figuresDir <- file.path(output_dir, "figures")
  if (!dir.exists(figuresDir)) { dir.create(figuresDir, recursive = TRUE) }
  fullFilePath <- file.path(figuresDir, paste0(currentDate, "_", fileName, ".png"))
  # Save the figure as PNG with dimensions in pixels
  png(file = fullFilePath, width = w, height = h, units = "px")
  # Render the plot
  print(figure)
  dev.off()
}
####Analysis Specific Code
#### Install packages if needed----
list.of.packages.cran <- c("tidyverse", "devtools", "annotate", "ggrepel", "pheatmap", "RColorBrewer", "EnhancedVolcano")
new.packages.cran <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages.cran)>0) install.packages(new.packages.cran)
# Install not-yet-installed Bioconductor packages
list.of.packages.bioc <- c("org.Hs.eg.db", "DESeq2", "apeglm","genefilter", "clusterProfiler", "GSVA", "DOSE", "xCell")
new.packages.bioc <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)>0)if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(new.packages.bioc, update = FALSE)
# Load packages
sapply(c(list.of.packages.cran, list.of.packages.bioc), require, character.only=TRUE)

#Get Inputs
colData_358_file<-paste0(input_dir,"/samplesheet_358.csv")
colData_358<-read.csv(file = colData_358_file, 
                      header=TRUE, 
                      stringsAsFactors = FALSE, 
                      check.names = FALSE,
                      row.names = "sample")
colData_358 <- colData_358[,c("cellLine", "drug")]
colData_358[, c("cellLine", "drug")] <- lapply(colData_358[, c("cellLine", "drug")], factor)
head(colData_358)

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

#### Calculations for DESeq2 ----
ddsObject_358 <- DESeqDataSetFromMatrix(countData = fc_358,
                                  colData = colData_358,
                                  design = ~ drug)

keep_358 <-  rowSums(counts(ddsObject_358))>= 10 #Keep genes which have at least 10 detected across samples
ddsObject_358 <- ddsObject_358[keep_358,]
#ddsObject_358$drug
ddsObject_358$drug <- relevel(ddsObject_358$drug, ref="Control")

dds_358 <- DESeq(ddsObject_358)
res_358 <- results(dds_358)
resultsNames(dds_358)

plotMA(res_358, ylim=c(-2, 2))


#Get Colors########
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#FF7B2A"))(99)


# Extract Transformed values
rld_358<-vst(dds_358) #estimate dispersion trend and apply a variance stabilizing transformation
# Here vst algorithm applied instead of rlog

## creating distance matrix
sampleDists_subset_358 <- as.matrix(dist(t(assay(rld_358))))
hm_358<-pheatmap::pheatmap(as.matrix(sampleDists_subset_358),
                           annotation_col = colData_358,
                           annotation_legend=TRUE,
                           color= colorRampPalette(c("blue","white","red"))(99), 
                           annotation_colors = list(drug=c("128-10"="#9FD900"
                                                           ,"128-13"="#FAA800"
                                                           ,"130"="#FC5AA2"
                                                           ,"Control"="#696969"),
                                                    cellLine=c(#"318"="#6500B1",
                                                      "358"="#006E18"
                                                    )),
                           #fontsize_row = 8, 
                           #cluster_cols = F,
                           #cluster_rows = T,
                           #scale = "row",
                           main="Distance Matrix 358 Cell Line"
                          )
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(min(sampleDists_subset_358), median(sampleDists_subset_358), max(sampleDists_subset_358)), c("blue", "white", "red"))
annotation_col <- data.frame(drug = colData_358$drug, cellLine = colData_358$cellLine)
rownames(annotation_col) <- rownames(colData_358)

# Define annotation colors
annotation_colors <- list(
  drug = c("128-10"="#9FD900", "128-13"="#FAA800", "130"="#FC5AA2", "Control"="#696969"),
  cellLine = c("358"="#006E18")
)

# Create HeatmapAnnotation object
ha <- HeatmapAnnotation(df = annotation_col, col = annotation_colors)

# Create and plot the heatmap
Heatmap(sampleDists_subset_358,
        name = "Distance",
        #top_annotation = ha,
        column_names_rot = 45,  # Rotate the column names by 45 degrees
        col = col_fun,
        show_row_names = T,
        show_column_names = T,
        heatmap_legend_param = list(title = "Distance"),
        column_title = "Distance Matrix 358 Cell Line")
# savePNG(hm_358, "DistanceMatrix-358")
# savePDF(hm_358, "DM")
## creating PCA plot
pca_358<-plotPCA(rld_358, intgroup="drug")
# pca<-pca + geom_text(aes(label=name),vjust=2, size = 3)
# saveFigure(figure=pca,fileName="PCAPlot",h=10,w=20)
#If needs to label samples using ggrepel
zz_358 <- plotPCA(rld_358, intgroup="drug", returnData=TRUE)
pca_repel_358 <- pca_358 + geom_label_repel(data=zz_358, aes(label=name))+
  ggtitle(label="PCA plot of 358 Cells")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

#saveFigure(figure=pca_repel358, fileName = "PCAPlot_Repel", h=15, w=20)
colData(dds_358)


#### Plotting DEG  (res has contrast data)
ComparisonColumn <- "drug"
factor1 <- "Control"
factor2 <- "130"
resultsNames(dds_358)
e <- as.character(c(ComparisonColumn, factor1, factor2))
res_358 <- lfcShrink(dds_358, contrast = e, type = "normal")

resdata_subset_358 <- merge(as.data.frame(res_358), as.data.frame(assay(rld_358)), by="row.names", sort=FALSE)
#write.csv(resdata_subset_358, paste0(output_dir,"/","DifferentialExpressionAnalysis_358", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset_358)[1] <- "Gene"
resdata_subset_358 <- resdata_subset_358[order(resdata_subset_358$padj),]
resdata_subset_358 <- resdata_subset_358[!is.na(resdata_subset_358$padj),]

#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db


resdata_subset_358$GeneSymbol<- mapIds(sp, keys=resdata_subset_358$Gene, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset_358$EntrezID<- mapIds(sp, keys=resdata_subset_358$Gene, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
# write.csv(resdata_subset_358, paste0(output_dir,"/","DifferentialExpressionAnalysis_cleaned_358", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

resdatagenes_subset_358 <- resdata_subset_358[complete.cases(resdata_subset_358$GeneSymbol),]
resdatagenes_subset_358 <- resdatagenes_subset_358[!duplicated(resdatagenes_subset_358$GeneSymbol), ]

resdatagenes_subset_358 <- resdatagenes_subset_358[order(resdatagenes_subset_358$log2FoldChange),]
Top50genesdown_subset_358<-head(resdatagenes_subset_358, 50)
Top50genesup_subset_358<-tail(resdatagenes_subset_358, 50)
Top_subset_358<-rbind(Top50genesdown_subset_358, Top50genesup_subset_358)


R_358 <- as.character(rownames(colData_358))
R_358<-c(R_358, "GeneSymbol")
Top_subset_358<-Top_subset_358%>% dplyr::select(all_of(R_358))
TopD_subset_358<-data.frame(Top_subset_358, check.names = FALSE)
#name the rows as genesymbols
rownames(TopD_subset_358)<-Top_subset_358$GeneSymbol

#Keep only the compared columns
# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, sep = "|"), names(TopD_subset_358), value = TRUE)

# Subset the data frame to keep only those columns
Top50_f1vsf2_358 <- TopD_subset_358[, columns_to_keep]


## top 50 fold changes
hmt50_358<-pheatmap::pheatmap(Top50_f1vsf2_358, scale="row", 
                          annotation_col=colData_358,
                          annotation_legend =TRUE, 
                          color= colorRampPalette(c("#0101FF","white","#FF0101"))(99),
                          annotation_colors = list(drug=c(#"128-10"="#9FD900",
                                                          #"128-13"="#FAA800",
                                                          "130"="#FC5AA2",
                                                          "Control"="#696969"),
                                                  cellLine=c(#"318"="#6500B1",
                                                            "358"="#006E18"
                                                  )),
                          fontsize_row = 8, 
                          cluster_cols = F,
                          cluster_rows = T,
                          main="Top50 FoldChange [Control vs 130] in 358 Cell Line")
# saveFigure(figure=hmt50358,fileName="Top50FoldChange_heatmap",h=12,w=12)


## variable genes
Ds_358<-resdatagenes_subset_358%>% dplyr::select(all_of(R_358), "GeneSymbol")
Ds_358<-Ds_358[!duplicated(Ds_358$GeneSymbol), ]
rownames(Ds_358)<-Ds_358$GeneSymbol
Ds_358<-Ds_358%>% dplyr::select(-"GeneSymbol")
#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds_358)),100)
mat_358 <- Ds_358[topVarGenes, ]
mat_358 <- mat_358 - rowMeans(mat_358)

#plot the variable genes in heatmap
hmVariable358<-pheatmap::pheatmap(mat_358,
                annotation_col=colData_358,
                color= colorRampPalette(c("#0101FF","white","#FF0101"))(99),
                annotation_colors = list(drug=c("128-10"="#9FD900"
                                                ,"128-13"="#FAA800"
                                                ,"130"="#FC5AA2"
                                                ,"Control"="#696969"), #48B2F9Blue
                                         cellLine=c(#"318"="#6500B1",
                                                    "358"="#006E18"
                                         )),
                annotation_legend =TRUE, 
                scale="row", 
                fontsize_row = 8,
                cluster_cols = F,
                cluster_rows = T,
                show_rownames=T, 
                main="Top100 Variable Genes in 358 Cell Line")
saveFigure(figure=hmVariable,fileName="Top100VariableGenes",h=12,w=12)

# Volcano Plot
vp_358<-EnhancedVolcano(resdata_subset_358,
                    lab = resdata_subset_358$GeneSymbol,x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 10e-4,
                    FCcutoff = 1)

#Pathway Analysis With Hallmark Geneset
human_hall_file<-paste0(input_dir,"/GSEA/h.all.v2023.2.Hs.symbols.gmt")
hall <- read.gmt(human_hall_file)

resdatagenes_gsea <- resdatagenes_subset_358
resdatagenes_gsea <- resdatagenes_gsea[!is.na(resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea <- resdatagenes_gsea[order(-resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea$FC_pval <- (resdatagenes_gsea$log2FoldChange)
logFC.l2n <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$FC_pval
names(logFC.l2n) <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$GeneSymbol

gsea.hall.l2n <- GSEA(logFC.l2n, TERM2GENE=hall, verbose=FALSE, pvalueCutoff=1)
gsea.hall.l2n.df <- as.data.frame(gsea.hall.l2n@result)
# write.csv(gsea.hall.l2n.df, paste0(outsdir,"/","GSEA_output", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

dp<-dotplot(gsea.hall.l2n, x="NES", showCategory=50, orderBy= "NES",font.size = 7) 
saveFigure(figure=dp,fileName="HallmarkPathwayAnalysis")