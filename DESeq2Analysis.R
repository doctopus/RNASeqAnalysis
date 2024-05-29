#DESeq2 Vignette

#Define Functions########################
# function to create project folder if not same as R Project folder and io folders in it 
#...project folder could be the R project or within it with the script & io within the project
setupProject <- function(project) {
  rstudio_dir <- rstudioapi::getActiveProject() # Get the Rstudio Directory
  rstudio_base <- basename(rstudio_dir) # Get the base name of the Rstudio Directory
  if (rstudio_base != project) { # If the base name != project, create a folder named as the project inside rstudio_dir
    project_dir <- file.path(rstudio_dir, project)
    if (!file.exists(project_dir)) {
      dir.create(project_dir, recursive = TRUE)
    }
  } else {
    # If the base name is the same as the project, use the rstudio_dir as the project_dir
    project_dir <- rstudio_dir
  }
  # Define input and output directories inside the project directory
  input_dir <- file.path(project_dir, "input")
  output_dir <- file.path(project_dir, "output")
  # Check if input and output directories exist, if not, create them
  if (!file.exists(input_dir)) {
    dir.create(input_dir, recursive = TRUE)
  }
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Set the working directory to the project directory
  setwd(project_dir)
  # Print Confirmation of Correct Folder Structure
  if (basename(getwd()) == project) {
    print("Folder setup correctly")
  } else {
    print("Fix folder structure")
  }
  # Export input_dir and output_dir to the global environment
  assign("input_dir", input_dir, envir = .GlobalEnv)
  assign("output_dir", output_dir, envir = .GlobalEnv)
}
savePDF <- function(figure, fileName, w = 7, h = 10) {
  currentDate <- format(Sys.Date(), "%Y%m%d") #current date in YYYYMMDD format
  # Define the directory for saving figures
  figuresDir <- file.path(output_dir, "figures")
  if (!dir.exists(figuresDir)) { dir.create(figuresDir, recursive = TRUE) }
  fullFilePath <- file.path(figuresDir, paste0(currentDate, "_", fileName, ".pdf"))
  # Save the figure
  pdf(file = fullFilePath, width = w, height = h, pointsize = 300 / 72)
  print(figure)
  dev.off()
}
savePNG <- function(figure, fileName, w = 900, h = 1300) {
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

####Analysis Specific Code ----
#Initiate project
setupProject("RNASeqAnalysis")
getwd()
# Project specific override: output folder V0.4 UnTrimmedNewIndexHierarchicalCategory
output_dir <- "/Users/i/Dropbox/Clinic3.0/Developer/RStudio/RNASeqAnalysis/output/v0.4_UnTrimmedNewIndexHierarchicalCategory"

#### Install packages if needed----
list.of.packages.cran <- c("annotate", "devtools", "EnhancedVolcano", "ggpubr", "ggrepel", "matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "viridis")
new.packages.cran <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages.cran)>0) install.packages(new.packages.cran)
# Install not-yet-installed Bioconductor packages
list.of.packages.bioc <- c("apeglm", "clusterProfiler", "DESeq2", "DOSE", "genefilter", "GSVA", "org.Hs.eg.db", "xCell")
new.packages.bioc <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)>0)if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(new.packages.bioc, update = FALSE)
# Load packages
sapply(c(list.of.packages.cran, list.of.packages.bioc), require, character.only=TRUE)

#rm(list=ls()) #remove all existing lists

#### Source & Process Input files ----
colData_file<-paste0(input_dir,"/samplesheet.csv")
fc_file<-paste0(input_dir,"/featureCounts_0")

##Process input files
colData<-read.csv(file = colData_file, 
                  header=TRUE, 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE,
                  row.names = "sample")
colData <- colData[,c("cellLine", "drug")]
colData[, c("cellLine", "drug")] <- lapply(colData[, c("cellLine", "drug")], factor)
head(colData)

fc <- read.delim(fc_file, row.names = NULL, check.names = FALSE)
genes_to_keep <- rowSums(fc[, 8:ncol(fc)]) > 1        #Keep Genes which are expressed 
desired_order <- rownames(colData)                    #Order of columns per colData

fc <- fc %>%
  filter(genes_to_keep) %>%                           # Keep genes with expression in >1 sample
  mutate(EnsembleID = gsub("\\..*$", "", Geneid)) %>% # Add column of removed dot of EnsembleID
  dplyr::select(7:ncol(.)) %>%                        # Remove the initial columns except length
  arrange(desc(Length)) %>%                           # Order the rows based on gene length
  filter(!duplicated(EnsembleID)) %>%                 # Remove duplicated EnsembleID keep 1st
  dplyr::select(-Length) %>%                          # Remove the Length column
  column_to_rownames("EnsembleID")                    # Set EnsembleID as row names

# Order the columns according to desired_order
fc <- fc[, desired_order]

head(fc)

#Check Duplicates in fc
rows_in_fc <- length(rownames(fc)) #sum(duplicated(rownames(fc)))#number of duplicate EnsembleIDs
#Number of complete cases based on the EnsembleID (rownames) (not empty, not NA)
completeEnsembleIDs <- sum(rownames(fc) != "" & !is.na(rownames(fc)))
if (rows_in_fc == completeEnsembleIDs) "Data OK! No duplicates" else paste(sum(duplicated(rownames(fc))),"Duplicates present")

#Make fc as the countData
countData <- fc

# fc<-read.delim(fc_file,
#                row.names=NULL,
#                check.names = FALSE)
# genes_to_keep <- rowSums(fc[,8:ncol(fc)]) > 1 #Keep rows where >1 genes expressed
# fc<-fc[genes_to_keep,]
# fc$EnsembleID<-gsub("\\..*$","",fc$Geneid) #+column of removed dot of EnsembleID
# fc<-fc[,7:ncol(fc)] #Remove the initial columns
# fc <- fc %>% arrange(desc(Length)) # Order the rows based on gene length
# fc<-fc[!duplicated(fc$EnsembleID),] #Remove duplicated rows based on EnsembleID keeping the first occurrence
# fc <- fc %>% dplyr::select(-Length) %>% column_to_rownames("EnsembleID")
# 
# desired_order <- rownames(colData)
# fc <- fc[, desired_order]



if (all(colnames(countData) %in% rownames(colData)) && 
    all(colnames(countData) == rownames(colData))) "Data ready for DESeq2" else "Data not ready"

#### Calculations for DESeq2 ----
ddsObject <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ cellLine + drug)

#Keeping rows that have at least 10 reads for a minimum number of samples
#Minimal number of samples is to specify the smallest group size, eg here 12 of each cellLine
# So 12 would be a good minimal number of samples.
# So row sum of minimal number of samples should have more than 10 reads at least
smallestGroupSize <- 12
counts(ddsObject)
#keep <-  rowSums(counts(dds2))>= 10
keep <-  rowSums(counts(ddsObject)>= 10) >= smallestGroupSize

ddsObject_filtered <- ddsObject[keep,]

ddsObject_filtered$drug

#Since we are primarily comparing between different drugs, so our primary level
#of comparison is drugs, and here reference level is Control
ddsObject_filtered$drug <- relevel(ddsObject_filtered$drug, ref="Control")

ddsObject_filtered$drug <- droplevels(ddsObject_filtered$drug) #remove the levels (of drug) 
# ...which do not have samples in the current data set. Here nothing removed
dds <- DESeq(ddsObject_filtered)

## Total number of raw counts per sample
colSums(counts(dds))
## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

#The results function without any argument will make log fold change
#...and p values for the last variable. In this case drugs
res <- results(dds)
head(res)
#We can optionally specify the coefficient or contrast to compare against
#========Need to do if needed for QC
#Shrinkage of Log Fold Change using apeglm algorithm
# We provide the name or number of the coefficient we want to shrink,
#..where the number refers to the order of the coefficient in the following command

resultsNames(dds) #coef in the lfcShrink should be from the output of this command.

resLFC <- lfcShrink(dds3, coef = "drug_128.10_vs_Control", type = "apeglm")

#Order the result table by smallest p value:
resOrdered <- res[order(res$pvalue),]
summary(res)
#summary(resOrdered) #Is same, it is only reordering

#How many adjusted p-values were less than 0.1?
sum(res$padj <0.1, na.rm = TRUE)


#To change the adjusted p-value cut off change alpha
result0.01 <- results(dds3, alpha=0.01)
summary(result0.01)
#How many adjusted p-values were less than 0.01?
sum(result0.01$padj <0.01, na.rm = TRUE)

###See the comparison group
resultsNames(dds)

#Specify comparison levels in results to get diferentially expressed results
#specificContrast <- results(dds3, contrast = c("drug", "Control", "128.10"))
#specificContrast
#Filter to find significant changes
#sigs <- na.omit (specificContrast)
#sigs <- sigs[sigs$padj <0.05,]
#============
#Counts of any single Gene (here, gene with min padj value)gene=which.min(res$padj) OR "ENSG00000197355"
plotCount <- plotCounts(dds, gene=which.min(res$padj), intgroup = "drug", returnData = TRUE)
ggplot(plotCount, aes(x=drug, y=count, color =drug))+
  geom_point(position = position_jitter(w=0.1, h=0))+
  geom_text_repel(aes(label = rownames(plotCount)))+
  theme_bw()+
  ggtitle("UAP1L1 Gene Expression")+
  theme(plot.title = element_text(hjust = 0.5))

#MA Plot
plotMA(res, ylim=c(-2, 2))
#plotMA(resOrdered)

#It is more useful to visualize the shrunken log2 fold changes
#..which remove the noise associated with log2 fold changes from low 
#..count genes without requiring arbitrary filtering
plotMA(resLFC)

#Estimate dispersion trend and apply variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#PCA Plot using VSD Transformation
plotPCAvst <- plotPCA(vsd, intgroup = "drug", returnData=FALSE)
plotPCAvstData <- plotPCA(vsd, intgroup="drug", returnData=TRUE)
pca_repel <- plotPCAvst+geom_label_repel(data=plotPCAvstData, aes(label=name))+
  ggtitle(label="PCA plot of All Cells")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

savePNG(figure=pca_repel, fileName = "PCAPlot_Repel", h=1200, w=900)
savePDF(figure=pca_repel, fileName = "PCAPlot_Repel", h=12, w=9)

#PCA Plot using RLOG Transformation
plotPCArld(rld, intgroup="drug") #Almost a similar plot

#Plot dispersion estimate
plotDispEsts(dds)

#Plot Heatmmap of Hierarchical Clustering
rld_mat <- assay(rld) #Extract the trasnformed matrix from the object
rld_cor <- cor(rld_mat) #Compute pairwise correlation values
head(rld_cor)
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation=colData, color = heat.colors, border_color = NA,
         fontsize = 10, fontsize_row = 10, height = 20)


#Get Colors########
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#FF7B2A"))(99)


# Extract Transformed values
rld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformationrld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformation

## creating distance matrix
sampleDists_subset <- as.matrix(dist(t(assay(rld))))
hm<-pheatmap::pheatmap(as.matrix(sampleDists_subset),
                       annotation_col = colData, 
                       col=my_colors,
                       annotation_legend=TRUE)

## creating PCA plot
pca<-plotPCA(rld, intgroup="drug")
# pca<-pca + geom_text(aes(label=name),vjust=2, size = 3)
# saveFigure(figure=pca,fileName="PCAPlot",h=10,w=20)
#If needs to label samples using ggrepel
library(ggrepel)
zz <- plotPCA(rld, intgroup="drug", returnData=TRUE)
pca_repel <- pca+geom_label_repel(data=zz, aes(label=name))+
  ggtitle(label="PCA plot of 318 Cells")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

#saveFigure(figure=pca_repel318, fileName = "PCAPlot_Repel", h=15, w=20)
colData(dds)


#### Plotting DEG
ComparisonColumn <- "drug"
factor1 <- "Control"
factor2 <- "130"
resultsNames(dds)
e <- as.character(c(ComparisonColumn, factor1, factor2))
res <- lfcShrink(dds, contrast = e, type = "normal")

resdata_subset <- merge(as.data.frame(res), as.data.frame(assay(rld)), by="row.names", sort=FALSE)
write.csv(resdata_subset, paste0(output_dir,"/","DifferentialExpressionAnalysis", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset)[1] <- "Gene" #Renaming EnsembleID column as Gene
# resdata_subset <- resdata_subset[order(resdata_subset$padj),]
# resdata_subset <- resdata_subset[!is.na(resdata_subset$padj),]
resdata_subset <- resdata_subset %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db


resdata_subset$GeneSymbol<- mapIds(sp, keys=resdata_subset$Gene, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset$EntrezID<- mapIds(sp, keys=resdata_subset$Gene, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
write.csv(resdata_subset, paste0(output_dir,"/","DifferentialExpressionAnalysis_cleaned", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

resdatagenes_subset <- resdata_subset[complete.cases(resdata_subset$GeneSymbol),]
resdatagenes_subset <- resdatagenes_subset[!duplicated(resdatagenes_subset$GeneSymbol), ]

resdatagenes_subset <- resdatagenes_subset[order(resdatagenes_subset$log2FoldChange),]

#↑Resdatagenes last modification

Top50genesdown_subset<-head(resdatagenes_subset, 50)
Top50genesup_subset<-tail(resdatagenes_subset, 50)
Top_subset<-rbind(Top50genesdown_subset, Top50genesup_subset)


R <- as.character(rownames(colData))
R<-c(R, "GeneSymbol")
Top_subset<-Top_subset %>% dplyr::select(all_of(R))
TopD_subset<-data.frame(Top_subset, check.names = FALSE)
#name the rows as genesymbols
rownames(TopD_subset)<-Top_subset$GeneSymbol

#Keep only the compared columns
# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, sep = "|"), names(TopD_subset), value = TRUE)

# Subset the data frame to keep only those columns
Top50_f1vsf2 <- TopD_subset[, columns_to_keep]


## top 50 fold changes
hmt50f1vsf2<-pheatmap::pheatmap(Top50_f1vsf2, scale="row", 
                              annotation_col=colData,
                              annotation_legend =TRUE, 
                              color= colorRampPalette(c("blue","white","red"))(99), 
                              annotation_colors = list(drug=c("128-10"="#48B2F9",
                                                             "128-13"="#FAA800",
                                                                "130"="#FC5AA2",
                                                             "Control"="#9FD900"),
                                                      cellLine=c("318"="#6500B1",
                                                                 "358"="#006E18")),
                              fontsize_row = 8, 
                              cluster_cols = F,
                              cluster_rows = T,
                              main="Top 50 FoldChange - Control vs 128-13")
saveFigure(figure=hmt50,fileName="Top50FoldChange_heatmap_Control_128-13",h=12,w=12)


## variable genes
Ds<-resdatagenes_subset%>% dplyr::select(all_of(R), "GeneSymbol")
Ds<-Ds[!duplicated(Ds$GeneSymbol), ]
rownames(Ds)<-Ds$GeneSymbol
Ds<-Ds%>% dplyr::select(-"GeneSymbol")
#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds)),100)
mat <- Ds[topVarGenes, ]
mat <- mat - rowMeans(mat)

#plot the variable genes in heatmap
hmVariable<-pheatmap::pheatmap(mat,
                                  annotation_col=colData,
                                  color= colorRampPalette(c("#1A1AFF","white","#FF1A1A"))(99), 
                                  annotation_colors = list(drug=c("128-10"="#48B2F9", 
                                                                  "128-13"="#FAC000",
                                                                      "130"="#FC5AA2",
                                                                  "Control"="#9FD900"),
                                                           cellLine=c("318"="#6500B1",
                                                                      "358"="#006E18")),
                                  annotation_legend =TRUE, 
                                  scale="row", 
                                  fontsize_row = 8, 
                                  show_rownames=T, 
                                  cluster_rows = T,
                                  cluster_cols = T,
                                  main="Top100 Variable Genes - All Cell Line")
saveFigure(figure=hmVariable,fileName="Top100VariableGenes All Cell Lines",h=12,w=12)

# Volcano Plot
vp<-EnhancedVolcano(resdata_subset,
                       lab = resdata_subset$GeneSymbol,x = 'log2FoldChange',
                       y = 'pvalue',
                       pCutoff = 10e-4,
                       FCcutoff = 1)

#Pathway Analysis With Hallmark Geneset
human_hall_file<-paste0(input_dir,"/GSEA/h.all.v2023.2.Hs.symbols.gmt")
hall <- read.gmt(human_hall_file)

resdatagenes_gsea <- resdatagenes_subset
resdatagenes_gsea <- resdatagenes_gsea[!is.na(resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea <- resdatagenes_gsea[order(-resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea$FC_pval <- (resdatagenes_gsea$log2FoldChange)
logFC.l2n <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$FC_pval
names(logFC.l2n) <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$GeneSymbol

gsea.hall.l2n <- GSEA(logFC.l2n, TERM2GENE=hall, verbose=FALSE, pvalueCutoff=1)

gsea.hall.l2n@result$Description <- gsub('HALLMARRK_', '', gsea.hall.l2n@result$Description)
gsea.hall.l2n@result$Description <- gsub('_', ' ', gsea.hall.l2n@result$Description)

# #If Needed to save
# gsea.hall.l2n.df <- as.data.frame(gsea.hall.l2n@result)
# write.csv(gsea.hall.l2n.df, paste0(outsdir,"/","GSEA_output", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

dotplot(gsea.hall.l2n, 
            x="NES", 
            showCategory=50, 
            orderBy= "NES",
            #title = "Cohort1 vs Cohort2",
            font.size = 7) +
  scale_color_gradient2(limits = c(0, 0.05),
                        breaks = c(0, 0.05),
                        labels = c("0.0", "0.05"),
                        low = "red",
                        mid = "red",
                        high = "red",
                        na.value = "blue") +
  ggtitle(paste(factor1, "vs", factor2)) +
  theme(plot.title = element_text(hjust = 0.5))

# saveFigure(figure=dp,fileName="HallmarkPathwayAnalysis")