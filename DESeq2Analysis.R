#DESeq2 Vignette

#Define Functions########################
# function to create project folder if not same as R Project folder and io folders in it 
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

# fx to save figures as 300dpi pdfs
savePDF <- function(figure, fileName, h = 7, w = 7, dpi = 300) {
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

savePNG <- function(figure, fileName, h = 1200, w = 900, dpi = 300) {
  currentDate <- format(Sys.Date(), "%Y%m%d") # current date in YYYYMMDD format
  # Define the directory for saving figures
  figuresDir <- file.path(output_dir, "figures")
  if (!dir.exists(figuresDir)) { dir.create(figuresDir, recursive = TRUE) }
  fullFilePath <- file.path(figuresDir, paste0(currentDate, "_", fileName, ".png"))
  # Save the figure as PNG with dimensions in pixels
  png(file = fullFilePath, height = h, width = w, units = "px", res = dpi)
  print(figure)
  dev.off()
}

# Example usage
# Assuming 'figure' is your plot object and 'output_dir' is defined
# savePNG(figure, "my_plot", h = 700, w = 700, dpi = 300)


############################
#Initiate project
setupProject("RNASeqAnalysis")
getwd()
# Project specific output for V0.4 UnTrimmedNewIndexHierarchicalCategory
output_dir <- "/Users/i/Dropbox/Clinic3.0/RStudio/RNASeqAnalysis/output/v0.4_UnTrimmedNewIndexHierarchicalCategory"

############################ Install packages
# BiocManager::install("apeglm")
# library(apeglm)
sapply(c("tidyverse", 
         "devtools", 
         "annotate", 
         "org.Hs.eg.db", 
         "DESeq2",
         "apeglm"), 
       require, character.only = TRUE)
############Get Input files
colData_file<-paste0(input_dir,"/samplesheet.csv")
fc_file<-paste0(input_dir,"/featureCounts_0")
########Process input files
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
  dplyr::select(7:ncol(.)) %>%                        # Remove the initial columns
  arrange(desc(Length)) %>%                           # Order the rows based on gene length
  filter(!duplicated(EnsembleID)) %>%                 # Remove duplicated EnsembleID keep 1st
  dplyr::select(-Length) %>%                          # Remove the Length column
  column_to_rownames("EnsembleID")                    # Set EnsembleID as row names

# Order the columns according to desired_order
fc <- fc[, desired_order]
head(fc)
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

if (all(colnames(fc) %in% rownames(colData)) && 
    all(colnames(fc) == rownames(colData))) "Data ready for DESeq2" else "Data not ready"

# all(colnames(fc) %in% rownames(colData))
# all(colnames(fc) == rownames(colData))

#Get Bioconductor Annotation Database
#sp <- org.Hs.eg.db
#Calculations for DESeq2 starts, need to repeat if groups change
ddsObject <- DESeqDataSetFromMatrix(countData = fc,
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

#Shrinkage of Log Fold Change using apeglm algorithm
# We provide the name or number of the coefficient we want to shrink,
#..where the number refers to the order of the coefficient in the following command
resultsNames(dds3)

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
resultsNames(dds3)

#Specify comparison levels in results to get diferentially expressed results
#specificContrast <- results(dds3, contrast = c("drug", "Control", "128.10"))
#specificContrast
#Filter to find significant changes
#sigs <- na.omit (specificContrast)
#sigs <- sigs[sigs$padj <0.05,]

#MA Plot
plotMA(res, ylim=c(-2, 2))
plotMA(resOrdered)

#It is more useful to visualize the shrunken log2 fold changes
#..which remove the noise associated with log2 fold changes from low 
#..count genes without requiring arbitrary filtering
plotMA(resLFC)
#Estimate dispersion trend and apply variance stabilizing transformation
#vsdata <- vst(dds3, blind = FALSE)

plotPCA(vsdata, intgroup = "group")

#Plot dispersion estimate
plotDispEsts(dds3)






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

## top 50 fold changes
hmt50318<-pheatmap::pheatmap(TopD_subset318[1:(length(TopD_subset318)-1)],scale="row", annotation_col=colData_318,
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
