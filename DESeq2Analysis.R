# Downstream Analysis of RNASeq

#### Define Functions ########################
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

#### Analysis Specific Code ----
#Initiate project
setupProject("RNASeqAnalysis")
getwd()
# Project specific override: output folder 1.1_HierarchicalCategory_Top100DEGs
output_dir <- paste0(output_dir, "/1.1_HierarchicalCategory_Top100DEGs")
#### Install & Load Packages ----
list.of.packages.cran <- c("annotate", "circlize", "devtools", "EnhancedVolcano", "ggpubr", "ggrepel", "matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "viridis")
new.packages.cran <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages.cran)>0) install.packages(new.packages.cran)
# Install not-yet-installed Bioconductor packages
list.of.packages.bioc <- c("apeglm", "clusterProfiler", "ComplexHeatmap", "DESeq2", "DOSE", "genefilter", "GSVA", "org.Hs.eg.db", "xCell")
new.packages.bioc <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)>0)if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(new.packages.bioc, update = FALSE)
# Load packages
sapply(c(list.of.packages.cran, list.of.packages.bioc), require, character.only=TRUE)

#### Source & Process Input files ----
colData_file<-paste0(input_dir,"/samplesheet.csv")
fc_file<-paste0(input_dir,"/featureCounts_0")

### Process input files
colData<-read.csv(file = colData_file, 
                  header=TRUE, 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE,
                  row.names = "sample")
colData <- colData[,c("cellLine", "drug")]
colData[, c("cellLine", "drug")] <- lapply(colData[, c("cellLine", "drug")], factor)
head(colData)

fc <- read.delim(fc_file, row.names = NULL, check.names = FALSE)
genes_to_keep <- rowSums(fc[, 8:ncol(fc)]) > 1        #Keep Genes which are expressed in >1 sample
desired_order <- rownames(colData)                    #Order of columns per colData

fc1 <- fc %>%
  filter(genes_to_keep) %>%                           # Keep genes with expression in >1 sample
  mutate(EnsembleID = gsub("\\..*$", "", Geneid)) %>% # Add column of removed dot of EnsembleID
  dplyr::select(7:ncol(.)) %>%                        # Remove the initial columns except length
  arrange(desc(Length)) %>%                           # Order the rows based on gene length
  filter(!duplicated(EnsembleID)) %>%                 # Remove duplicated EnsembleID keep 1st
  dplyr::select(-Length) %>%                          # Remove the Length column
  column_to_rownames("EnsembleID")                    # Set EnsembleID as row names

# Order the columns according to desired_order
fc2 <- fc1[, desired_order]

head(fc2)

#Check Duplicates in fc
rows_in_fc <- length(rownames(fc2)) #sum(duplicated(rownames(fc)))#number of duplicate EnsembleIDs
#Number of complete cases based on the EnsembleID (rownames) (not empty, not NA)
completeEnsembleIDs <- sum(rownames(fc2) != "" & !is.na(rownames(fc2)))
if (rows_in_fc == completeEnsembleIDs) {
  paste("Data OK! No Empty/NA. Total Duplicates:", sum(duplicated(rownames(fc2))))
} else {
  paste(sum(duplicated(rownames(fc2))), "Duplicates present")
}


#Make fc as the countData
countData <- fc2

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
##Prepare DESeq Object
ddsObject <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ cellLine + drug)

#Keeping rows that have at least 10 reads for a minimum number of samples
#Minimal number of samples is the smallest group size, eg here 12 of each cellLine
# So row sum of 12 samples should have more than 10 reads at least
smallestGroupSize <- 12
#counts(ddsObject)
#keep <-  rowSums(counts(dds2))>= 10
keep <-  rowSums(counts(ddsObject)>= 10) >= smallestGroupSize

ddsObject_filtered <- ddsObject[keep,]

ddsObject_filtered$drug

#Since we are primarily comparing between different drugs, so our primary level
#of comparison is drugs, and here reference level is Control
ddsObject_filtered$drug <- relevel(ddsObject_filtered$drug, ref="Control")

ddsObject_filtered$drug <- droplevels(ddsObject_filtered$drug) #remove the levels (of drug) 
# ...which do not have samples in the current data set. Here nothing removed

##Run DESeq2 analysis
dds <- DESeq(ddsObject_filtered)

##DDS Quality Control (Optional)----
## Total number of raw counts per sample
colSums(counts(dds))
## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

#The results function without any argument will make log fold change & p values
#...for the last variable in design formula, between defined reference (Control) and last factor in it.
res <- results(dds)
head(res)

# The result function will by default pull the last variable in design (drug)
#...unless "name" or "contrast" argument is provided as follows.

#To define contrast pairs in the results argument as follows
resControlx128.10 <- results(dds, contrast = c("drug", "128-10", "Control"))
head(resControlx128.10)

#Alternatively pick a cobination from 
resultsNames(dds)

resControlx128.13 <- results(dds, name="drug_128.13_vs_Control")
head(resControlx128.13)

# Visualize the result to detect distribution anomaly
plotMA(res, ylim = c(-3, 3))
hist(res$pvalue, breaks=20, col="grey")
hist(res$padj, breaks=20, col="grey")

# Since we want rank and visualization, shrink data by either apeglm, ashr or "normal" algorithm
# We provide the name or number of the coefficient we want to shrink,
#..where the number refers to the order of the coefficient in resultsNames(dds)
# coef in the lfcShrink should be from the output of resultsNames(dds)

# See the comparison groups
resultsNames(dds)

resLFCapeglm <- lfcShrink(dds, coef = "drug_130_vs_Control", type = "apeglm")
resLFCnormal <- lfcShrink(dds, coef = "drug_130_vs_Control", type = "normal")

#Plot dispersion estimate
plotDispEsts(dds)

#MA Plot
plotMA(res)
#plotMA(resOrdered)

#It is more useful to visualize the shrunken log2 fold changes
#..which remove the noise associated with log2 fold changes from low 
#..count genes without requiring arbitrary filtering
plotMA(resLFCapeglm)
plotMA(resLFCnormal)

#Order the result table by smallest p value:
resOrdered <- res[order(res$pvalue),]
summary(res)
#summary(resOrdered) #Is same, it is only reordering

#How many adjusted p-values were less than 0.1?
sum(res$padj <0.1, na.rm = TRUE)


#To change the adjusted p-value cut off change alpha
result0.01 <- results(dds, alpha=0.01)
summary(result0.01)
#How many adjusted p-values were less than 0.01?
sum(result0.01$padj <0.01, na.rm = TRUE)


#Specify comparison levels in results to get deferentially expressed results
#specificContrast <- results(dds3, contrast = c("drug", "Control", "128.10"))
#specificContrast
#Filter to find significant changes
#sigs <- na.omit (specificContrast)
#sigs <- sigs[sigs$padj <0.05,]


#Counts of any single Gene (here, gene with min padj value)gene=which.min(res$padj) OR "ENSG00000197355" 
#plotCount <- plotCounts(dds, gene=which.min(res$padj), intgroup = "drug", returnData = TRUE)
plotCount <- plotCounts(dds, gene="ENSG00000276600", intgroup = "drug", returnData = TRUE) #ENSG00000276600 is RAB7B Low expressed in drug
ggplot(plotCount, aes(x=drug, y=count, color =drug))+
  geom_point(position = position_jitter(w=0.1, h=0))+
  geom_text_repel(aes(label = rownames(plotCount)))+
  theme_bw()+
  ggtitle("UAP1L1 Gene Expression")+
  theme(plot.title = element_text(hjust = 0.5))

##DDS Calculation ----
# Apply transformation & estimate dispersion trend
# RLT: Regularized Log Transformation
# VST: Variance Stabilizing Transformation

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 2)
head(assay(rld), 2)

##Heatmaps ----
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#A71415"))(99)

# Extract Transformed values
#rld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformationrld<-vst(dds)
# Using rld transformed data for the analysis [rld <- rlog(dds)]
length(rld)

## creating distance matrix (Not Needed)
sampleDists_subset <- as.matrix(dist(t(assay(rld))))
hm<-pheatmap::pheatmap(as.matrix(sampleDists_subset),
                       annotation_col = colData, 
                       col=my_colors,
                       annotation_legend=TRUE)

#Plot Heatmap of Hierarchical Clustering
rld_mat <- assay(rld) #Extract the transformed matrix from the object
rld_cor <- cor(rld_mat) #Compute pairwise correlation values
head(rld_cor)
heat.colors <- RColorBrewer::brewer.pal(6, "BrBG")
pheatmap(rld_cor, annotation=colData, color = heat.colors, border_color = NA,
         fontsize = 10, fontsize_row = 10, height = 20)


##PCA Plots ----
## PCA Plot using VSD Transformation 
# Argument in plotPCA ntop=length(vsd)) to include all genes; 
#...default is ntop=500, to plot top 500 feature variance
plotPCAvst <- plotPCA(vsd, intgroup = "drug", returnData=FALSE, ntop=length(vsd)) 
plotPCAvstData <- plotPCA(vsd, intgroup="drug", returnData=TRUE, ntop=length(vsd))
pcaVST <- plotPCAvst+geom_label_repel(data=plotPCAvstData, aes(label=name))+
  ggtitle(label="PCA plot of All Cells: VST Transformed Data")+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

#PCA Plot using RLT Transformation #Almost a similar plot
plotPCArld <- plotPCA(rld, intgroup="drug", returnData=FALSE, ntop=length(rld)) 
plotPCArldData <- plotPCA(rld, intgroup="drug", returnData=TRUE, ntop=length(rld)) 
pcaRLDAll <- plotPCArld+geom_label_repel(data=plotPCArldData, 
                                      aes(label=name), 
                                      min.segment.length = 0.5)+
  ggtitle(label="PCA Plot of All Genes")+
  labs(caption = "Regularized Log Transformed (RLT) Count Data") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.caption = element_text(hjust=1, size = '8', color = 'grey', face = 'italic'))
print(pcaRLDAll)

## PCA Plot using RLD for default top 500 genes
plotPcaRLD500<-plotPCA(rld, intgroup="drug", returnData=FALSE)
plotPcaRLD500Data <- plotPCA(rld, intgroup="drug", returnData=TRUE)
pcaRLD500 <- plotPcaRLD500 + geom_label_repel(data=plotPcaRLD500Data, 
                                    aes(label=name),
                                    min.segment.length = 0.5,
                                    max.overlaps = 1)+
  ggtitle(label="PCA Plot of Top 500 Variable Genes")+
  labs(caption = "Regularized Log Transformed (RLT) Count Data") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5,  size = 14, face = "bold"),
        plot.caption = element_text(hjust=1, size = '8', color = 'grey', face = 'italic'),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #legend.position = "bottom",
        #strip.background = element_blank(),
        axis.text = element_text(size = 10)
        #panel.border = element_rect(color="black"),
        )
print(pcaRLD500)
#saveFigure(figure=pca_repel318, fileName = "PCAPlot_Repel", h=15, w=20)
#colData(dds)


##Plotting DEG ----
ComparisonColumn <- "drug"
factor1 <- "128-10"
factor2 <- "Control"
# resultsNames(dds)
e <- as.character(c(ComparisonColumn, factor1, factor2))
res <- lfcShrink(dds, contrast = e, type = "normal") #Calculates LFC and p-values

annotation_colors <- list(
  drug = c("128-10"="#9FD900", 
           "128-13"="#FAA800", 
              "130"="#ff5d8f", 
          "Control"="#6c757d"),
  
  cellLine =c("358"="#006E18", 
            "318" = "#832161")
)
icolors <- colorRampPalette(c("blue",
                              "white",
                              "red"))(99)

title <- paste(factor1, "vs", factor2)  

#TODO X and Y are not same. Investigate: "The assay function extracts the matrix of normalized values
# x <- as.data.frame(assay(rld)) #Transformed Normalized
# y <- as.data.frame(counts(dds, normalized=TRUE)) #Only Normalized
# 
# r <- as.data.frame(res) #Shrinked Result (L2FC, padj etc) based on "normal" algorithm


#==========AmCodeStart (Heatmap Based On Counts)----
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)),
                 by = "row.names",
                 sort = FALSE)
names(resdata)[1] <- "EnsembleID"
resdata <- resdata[order(resdata$padj),]
resdata <- resdata[complete.cases(resdata$padj),]
#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db

resdata$Gene<- mapIds(sp, keys=resdata$EnsembleID, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")

resdata <- resdata[complete.cases(resdata$Gene),]
resdata <- resdata[!duplicated(resdata$Gene),]
resdata <- resdata[order(resdata$log2FoldChange),]

#Heatmap of Top/Botttom 25 Log Fold Change Genes
Top25Down <- head(resdata, 25)
Top25Up <- tail(resdata, 25)
Top25 <- rbind(Top25Down, Top25Up)

A <- as.character(rownames(colData))
AG<-c(A, "Gene")
Top25<-Top25 %>% dplyr::select(all_of(AG))
Top25<-data.frame(Top25, check.names = FALSE)
rownames(Top25)<-Top25$Gene
Top25 <- Top25 %>% dplyr::select(all_of(A))

#Keep only the compared columns
# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, sep = "|"), names(Top25), value = TRUE)

# Subset the data frame to keep only those columns
Top25 <- Top25[, columns_to_keep] #%>% as.matrix(.)

# Filter the colData to keep only those rows
colData25 <- colData[columns_to_keep, ]

# Order the columns by 'drug'
ordered_indices <- order(colData25$drug, decreasing = TRUE)
Top25_ordered <- Top25[, ordered_indices]
colData25_ordered <- colData25[ordered_indices, ]

# Filter annotation colors to include only relevant levels
annotation_colors_filtered <- list(
  drug = annotation_colors$drug[names(annotation_colors$drug) %in% unique(colData25_ordered$drug)],
  cellLine = annotation_colors$cellLine[names(annotation_colors$cellLine) %in% unique(colData25_ordered$cellLine)]
)
pheatmap::pheatmap(Top25_ordered,
         scale="row",
         annotation=colData25_ordered,
         cluster_cols=F,
         cluster_rows=T,
         annotation_colors = annotation_colors_filtered,
         cellwidth=10,
         annotation_legend =TRUE,
         color= icolors,
         fontsize_row = 5,
         main=title)

## Complex Heat Map Code--
# library(ComplexHeatmap)
# library(circlize)

# col_fun <- colorRamp2(c(min(Top25), median(Top25), max(Top25)), c("blue", "white", "red"))
# Adjust color mapping to ensure proper visualization
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
title <- paste("Top 50 Fold Change:", factor1, "vs", factor2)


# Ensure the 'drug' column in colData is a factor with the correct order
#colData$drug <- factor(colData25$drug, levels = c("Control", "130"))
# Order the columns by 'drug'
ordered_indices <- order(colData25$drug, decreasing = TRUE)
Top25_ordered <- Top25[, ordered_indices]
colData25_ordered <- colData25[ordered_indices, ]

# Filter annotation colors to include only relevant levels
annotation_colors_filtered <- list(
  drug = annotation_colors$drug[names(annotation_colors$drug) %in% unique(colData25_ordered$drug)],
  cellLine = annotation_colors$cellLine[names(annotation_colors$cellLine) %in% unique(colData25_ordered$cellLine)]
)

# Create the combined column annotations with drug above cellLine
ha_combined <- HeatmapAnnotation(
  df = data.frame(drug = colData25_ordered$drug, cellLine = colData25_ordered$cellLine),
  col = annotation_colors_filtered,
  annotation_name_side = "right",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  show_legend = T # Disable legends in the annotation to custom order the legends using library(grid)
)

# Scale the rows
Top25_scaled <- t(scale(t(Top25_ordered)))


# Create the heatmap
hm25 <- Heatmap(
  Top25_scaled,
  name = "Expression",
  cluster_rows = T,
  cluster_columns = F,
  col = col_fun,
  top_annotation = ha_combined,
  column_title = title,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 10),
  #column_split = colData25$drug,
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "grey", lwd = 0.5))
  }
) 
# Optional: Customize the legend orders
# library(grid)
# Create custom legends for each annotation
legend_drug <- Legend(
  at = rev(names(annotation_colors_filtered$drug)),
  labels = rev(names(annotation_colors_filtered$drug)),
  title = "Drug",
  legend_gp = gpar(fill = rev(annotation_colors_filtered$drug))
)

legend_cellLine <- Legend(
  at = names(annotation_colors_filtered$cellLine),
  labels = names(annotation_colors_filtered$cellLine),
  title = "Cell Line",
  legend_gp = gpar(fill = annotation_colors_filtered$cellLine)
)

# Combine the legends in the desired order
legends <- packLegend(legend_drug, legend_cellLine)

# Draw the heatmap with custom legends
draw(hm25, heatmap_legend_side = "right", annotation_legend_side = "right", annotation_legend_list = legends)

#==========AmCodeEnd----
#==========CodesStarts (Based on RlogTransfrmed Data) ----
resdata_subset <- merge(as.data.frame(res), as.data.frame(assay(rld)), 
                        by="row.names", 
                        sort=FALSE)
# write.csv(resdata_subset, paste0(output_dir,"/","DifferentialExpressionAnalysis", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset)[1] <- "EnsembleID" #Renaming Row.names column as EnsembleID
# resdata_subset <- resdata_subset[order(resdata_subset$padj),]
# resdata_subset <- resdata_subset[!is.na(resdata_subset$padj),]
resdata_subset <- resdata_subset %>%
  filter(!is.na(padj)) %>% #Remove rows with padj as na
  arrange(padj) #Arrange by increasing padj (Decreasing significance)

#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db


resdata_subset$GeneSymbol<- mapIds(sp, keys=resdata_subset$EnsembleID, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset$EntrezID<- mapIds(sp, keys=resdata_subset$EnsembleID, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
# write.csv(resdata_subset, paste0(output_dir,"/","DifferentialExpressionAnalysis_cleaned", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

# resdata <- resdata[complete.cases(resdata$Gene),]
# resdata <- resdata[!duplicated(resdata$Gene),]
# resdata <- resdata[order(resdata$log2FoldChange),]

resdatagenes_subset <- resdata_subset[complete.cases(resdata_subset$GeneSymbol),] #Keep where gene symbol present
resdatagenes_subset <- resdatagenes_subset[!duplicated(resdatagenes_subset$GeneSymbol), ] #Remove duplicates #TODO-Decide Order
resdatagenes_subset <- resdatagenes_subset[order(resdatagenes_subset$log2FoldChange),] #Order ascending L2FC

#↑Resdatagenes last modification
### Top50 Fold Changes ----
Top50genesdown_subset<-head(resdatagenes_subset, 50)
Top50genesup_subset<-tail(resdatagenes_subset, 50)
Top_subset<-rbind(Top50genesdown_subset, Top50genesup_subset)


R <- as.character(rownames(colData))
R<-c(R, "GeneSymbol")
Top_subset<-Top_subset %>% dplyr::select(all_of(R))
TopD_subset<-data.frame(Top_subset, check.names = FALSE)
rownames(TopD_subset)<-Top_subset$GeneSymbol

#Keep only the compared columns
# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, sep = "|"), names(TopD_subset), value = TRUE)

# Subset the data frame to keep only those columns
Top50_f1vsf2 <- TopD_subset[, columns_to_keep]

# Filter the colData to keep only those rows
colData50 <- colData[columns_to_keep, ]

# Order the columns by 'drug'
ordered_indices <- order(colData50$drug, decreasing = TRUE)
Top50_f1vsf2_ordered <- Top50_f1vsf2[, ordered_indices]
colData50_ordered <- colData50[ordered_indices, ]


#Save for Venn Diagram Analysis
row_names <- rownames(Top50_f1vsf2_ordered)
# T50Control_130 <- data.frame(matrix(ncol = 1, nrow = length(row_names)))
# T50Control_130 <- data.frame(ColumnNames = column_names)

#T50Control_128.10 <- data.frame(Control_128.10 = row_names)
#T50Control_128.13 <- data.frame(Control_128.13 = row_names)
T50Control_130 <- data.frame(Control_130 = row_names)
# head(T50Control_130)

# Filter annotation colors to include only relevant levels
annotation_colors_filtered <- list(
  drug = annotation_colors$drug[names(annotation_colors$drug) %in% unique(colData50_ordered$drug)],
  cellLine = annotation_colors$cellLine[names(annotation_colors$cellLine) %in% unique(colData50_ordered$cellLine)]
)
## Top 50 fold changes
hmt50f1vsf2<-pheatmap::pheatmap(Top50_f1vsf2_ordered, scale="row", 
                              annotation_col=colData50_ordered,
                              annotation=colData50_ordered, 
                              cluster_cols = F,
                              cluster_rows = T,
                              cellwidth=10,
                              annotation_legend =TRUE, 
                              color= icolors,
                              annotation_colors = annotation_colors_filtered,
                              fontsize_row = 8,
                              main=paste("Top 50 Fold Change:", title))
# saveFigure(figure=hmt50,fileName="Top50FoldChange_heatmap_Control_128-13",h=12,w=12)


### Variable Genes ----
Ds<-resdatagenes_subset%>% dplyr::select(all_of(R), "GeneSymbol")
Ds<-Ds[!duplicated(Ds$GeneSymbol), ]
rownames(Ds)<-Ds$GeneSymbol
Ds<-Ds%>% dplyr::select(-"GeneSymbol")

#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds)),100)
mat <- Ds[topVarGenes, ]
mat <- mat - rowMeans(mat)

#Optionally Show Only Once Cell Line
cellLineToPlot <- "358"
# Find columns with "318" in their names
DsKeepCols <- grep(cellLineToPlot, names(Ds), value = TRUE)
# Subset the data frame to keep only those columns
DsFiltered <- Ds[, DsKeepCols]
topVarGenes <- head(order(-genefilter::rowVars(DsFiltered)),100)
mat <- DsFiltered[topVarGenes, ]
mat <- mat - rowMeans(mat)

#plot the variable genes in heatmap
hmVariable<-pheatmap::pheatmap(mat,
                      annotation_col=colData,
                      color= icolors, 
                      annotation_colors = annotation_colors,
                      annotation_legend =TRUE, 
                      scale="row", 
                      fontsize_row = 8, 
                      show_rownames=T, 
                      cluster_rows = T,
                      cluster_cols = F,
                      main="Top 100 Variable Genes - 358 Cell Line")
# saveFigure(figure=hmVariable,fileName="Top100VariableGenes All Cell Lines",h=12,w=12)

### Volcano Plot ----
vp<-EnhancedVolcano(resdata_subset,
                    lab = resdata_subset$GeneSymbol,
                    title = paste(title),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 10e-4, #Default P value cutoff is 10e-6
                    FCcutoff = 1, #Default log2FC cutoff is >|2|
                    pointSize = 3.0,
                    labSize = 6.0,
                    drawConnectors = TRUE)
print(vp)
########################################
#Pathway Analysis With Hallmark Geneset
########################################
#SKIP-STARTS=== Skip this section if part of the same code
  #GET-DDS
  ddsObject <- DESeqDataSetFromMatrix(countData = countData,
                                      colData = colData,
                                      design = ~ cellLine + drug)
  #Keep roows with >10 reads per smallest group of analysis
  smallestGroupSize <- 12
  #keep <-  rowSums(counts(dds2))>= 10
  keep <-  rowSums(counts(ddsObject)>= 10) >= smallestGroupSize
  ddsObject_filtered <- ddsObject[keep,]
  #Relevel the Drug to define Control
  ddsObject_filtered$drug <- relevel(ddsObject_filtered$drug, ref="Control")
  #Drop the level of drug which does not have samples (if any)
  ddsObject_filtered$drug <- droplevels(ddsObject_filtered$drug) #remove the levels (of drug) 
  dds <- DESeq(ddsObject_filtered)
  #GET-RLD
  # Apply transformation & estimate dispersion trend
  # RLT: Regularized Log Transformation
  # VST: Variance Stabilizing Transformation
  #vsd <- vst(dds, blind = FALSE)
  rld <- rlog(dds, blind=FALSE)


  
ComparisonColumn <- "drug"
factor1 <- "128-10"
factor2 <- "Control"
title <- paste(factor1, "vs", factor2)
# resultsNames(dds)
e <- as.character(c(ComparisonColumn, factor1, factor2))
#Shrink dds based on comparison data
res <- lfcShrink(dds, contrast = e, type = "normal")
#Combine shrunk results (with lfc, padj etc with normalized count data)
resdata <- merge(as.data.frame(res), #Get the counts data alongwith the comparison results
                 #as.data.frame(counts(dds, normalized=FALSE)), #Get original count data for each set
                 as.data.frame(assay(rld)), #Get transformed count values for each set
                 by = "row.names",
                 sort = FALSE)
  
names(resdata)[1] <- "EnsembleID" #Rename the first column (Row.names) as EnsembleID
#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db

#Add a column of Gene translated from EnsembleID
resdata$Gene<- mapIds(sp, 
                      keys=resdata$EnsembleID, 
                      column=c("SYMBOL"), 
                      keytype="ENSEMBL", 
                      multiVals="first")

resdata <- resdata[order(resdata$padj),] #Order by ascending p-adj significance;low p =significant at top
resdata <- resdata[complete.cases(resdata$padj),] #Keep rows which have p-adj data
resdata <- resdata[complete.cases(resdata$Gene),] #Remove rows that don't have assigned genes
resdata <- resdata[!duplicated(resdata$Gene),] #Remove rows with same gene, keeping significant ones (low padj)
resdata <- resdata[order(resdata$log2FoldChange),] #Now order per log2FoldChange
#Rename dataset to match codeS
resdatagenes_subset <- resdata
#Pick top and bottom most significant DEGs
#Skipped
#SKIP-ENDS=== Skip this section if part of the same code


#Pathway Analysis
human_hall_file<-paste0(input_dir,"/GSEA/h.all.v2023.2.Hs.symbols.gmt")
# human_hall_file<-paste0(input_dir,"/GSEA/h.all.v7.2.symbols.gmt")
genesets <- read.gmt(human_hall_file)

resdata_gsea <- resdatagenes_subset #For CodeS
# resdata_gsea <- resdata
resdata_gsea <- resdata_gsea[!is.na(resdata_gsea$log2FoldChange),] #Remove rows without log2fc data
ordered_data <- resdata_gsea[order(-resdata_gsea$log2FoldChange),] #Order by decreasing log2fc

logFC <- ordered_data$log2FoldChange #Take the log2FC column only
# names(logFC) <- ordered_data$Gene #Assign Row names as their corresponding Gene
names(logFC) <- ordered_data$GeneSymbol #CodeS #Assign Row names as their corresponding Gene

#If need to plot a known list of genes (intersection123) NOT POSSIBLE SINCE LIST IS TOO SMALL
#intersectionDEGs <- intersection123 #Intersection Genes among Top DEG between comparison groups
#logFC_filtered <- logFC[names(logFC) %in% intersectionDEGs] #keep rows where gene names are in intersectionDEGs

#print(logFC_filtered[1:10])

gsea_hallmark <- GSEA(logFC, 
                      TERM2GENE=genesets, 
                      verbose=TRUE, 
                      pvalueCutoff=1)

gsea_hallmark@result$Description <- gsub('HALLMARK_','',gsea_hallmark@result$Description)
gsea_hallmark@result$Description <- gsub('_',' ',gsea_hallmark@result$Description)

#Need LogFC in descending order along with the gene names
dotplot(gsea_hallmark,
        x = "NES",
        showCategory = 50,
        orderBy = "NES",
        #title = paste("Significant Differential Pathways:", title),
        color = "p.adjust", #Map 'p.adjust' to color
        font.size = 7) +
  #scale_color_gradient(low="red4", high = "tomato")
  scale_color_gradient2(limits = c(0, 0.9),
    breaks = c(0, 0.3),
    labels = c("0.0", "0.9"),
    low = "orange",
    mid = "red",
    high = "green",
    na.value = "blue") +
  ggtitle(paste("Significant Differential Pathways:", title)) +
  theme(plot.title = element_text(hjust = 0.5))
  
#NES is Normalized Enrichment Score

#Save Data for Venn Diagram
pathways <- gsea_hallmark@result$Description 
#pathways <- rownames(Top50_f1vsf2_ordered)
Path23Control_130 <- data.frame(Control_130 = pathways)


####OLD CODE FOLLOWS
#   
# human_hall_file<-paste0(input_dir,"/GSEA/h.all.v2023.2.Hs.symbols.gmt")
# hall <- read.gmt(human_hall_file)
# 
# resdatagenes_gsea <- resdatagenes_subset
# resdatagenes_gsea <- resdatagenes_gsea[!is.na(resdatagenes_gsea$log2FoldChange),]
# resdatagenes_gsea <- resdatagenes_gsea[order(-resdatagenes_gsea$log2FoldChange),]
# resdatagenes_gsea$FC_pval <- (resdatagenes_gsea$log2FoldChange)
# logFC.l2n <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$FC_pval
# names(logFC.l2n) <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$GeneSymbol
# 
# gsea.hall.l2n <- GSEA(logFC.l2n, TERM2GENE=hall, verbose=FALSE, pvalueCutoff=1)
# 
# gsea.hall.l2n@result$Description <- gsub('HALLMARK_', '', gsea.hall.l2n@result$Description)
# gsea.hall.l2n@result$Description <- gsub('_', ' ', gsea.hall.l2n@result$Description)
# 
# # #If Needed to save
# # gsea.hall.l2n.df <- as.data.frame(gsea.hall.l2n@result)
# # write.csv(gsea.hall.l2n.df, paste0(outsdir,"/","GSEA_output", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
# 
# dotplot(gsea.hall.l2n, 
#             x="NES", 
#             showCategory=50, 
#             orderBy= "NES",
#             #title = "Cohort1 vs Cohort2",
#             color="p.adjust",   # Map 'p.adjust' to color
#             font.size = 7) +
#   scale_color_gradient2(#limits = c(0, 0.05),
#                         #breaks = c(0, 0.05),
#                         #labels = c("0.0", "0.05"),
#                         low = "orange",
#                         mid = "red",
#                         high = "green",
#                         na.value = "blue") +
#   ggtitle(paste(factor1, "vs", factor2)) +
#   theme(plot.title = element_text(hjust = 0.5))

#NES is Normalized Enrichment Score
# saveFigure(figure=dp,fileName="HallmarkPathwayAnalysis")

### Go Pathway Analysis ----
library(clusterProfiler)
library(enrichplot)
# Use the example data set included with the package DOSE
data(geneList, package="DOSE")
head(geneList, 10)
# Set fold change > 2 as being DE genes
gene <- names(geneList)[abs(geneList)>2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

ggo <- groupGO(gene = gene,
             OrgDb = org.Hs.eg.db,
             ont = "CC",
             level = 3,
             readable = TRUE)

head(ggo)


ego <- enrichGO(gene = gene,
                universe = names(geneList),
                OrgDb = "org.Hs.eg.db", 
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable=TRUE)
head(ego)

median.FC.values <- runif(n = dim( ego@result)[1] , min=-15, max = 15)
names(median.FC.values) <- rownames(ego@result)

ego@result$medianFC <-  median.FC.values

head(ego@result)
head(as.data.frame(ego))

options(enrichplot.colors = c("pink", "blue"))
dotplot(ego, color="p.adjust") +
  scale_colour_gradient2(low="green", mid="white", high="red",
                         limits = c(-15, 15),
                         breaks = c(-15, -7.5, 0, 7.5, 15),
                         labels = c("-15down", "-7.5", "0", "7.5", "15up"),  
                         guide=guide_colorbar(reverse=TRUE) )  +
  labs(size="Count", colour="Median logFC")

############################################
## Venn Diagram & Upset Plot of Top 50 Genes
############################################
###Prepare data for Venn and Upset Plots----
T50Merged <- cbind(T50Control_128.10, T50Control_130, T50Control_128.13)
listInput <- list(
  Control_128.10 = T50Merged$Control_128.10,
  Control_130 = T50Merged$Control_130,
  Control_128.13 = T50Merged$Control_128.13
)
upsetData <- data.frame(listInput)
###Venn Diagram===----
#install.packages("ggVennDiagram")
library(ggVennDiagram)
# Create the Venn diagram #DONT CREATE LIST
# venn_data <- list(
#   "Control vs 128.10" = Control_128.10,
#   "Control vs 130" = Control_130,
#   "Control vs 128.13" = Control_128.13
# )

## Plot using ggVennDiagram===
#ggVennDiagram(venn_data, label_alpha = 0.5) +
ggVennDiagram::ggVennDiagram(listInput, label_alpha = 0.5) +
#ggVennDiagram::ggVennDiagram(upsetData, label_alpha = 0.5) + #Works too
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Venn Diagram of Top50 DEG Between Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )  # Change background color

# Find intersections
# intersection_12 <- intersect(Control_128.10, Control_130)
# intersection_13 <- intersect(Control_128.10, Control_128.13)
# intersection_123 <- Reduce(intersect, list(Control_128.10, Control_130, Control_128.13))

# intersection_all <- Reduce(intersect, list(set1, set2, set3, set4))
# Print intersections
print(paste("Intersection of Set1 and Set2:", toString(intersection_12)))
print(paste("Intersection of Set1, Set2, and Set3:", toString(intersection_123)))
print(paste("Intersection of all sets:", toString(intersection_all)))

intersection123 <- Reduce(intersect, listInput)
print(paste("Intersection", toString(intersection123)))

##Pathway Venn Diagram
PathwayMerged <- cbind(Path23Control_128.10, Path23Control_128.13, Path23Control_130)
pathwayVennInput <- list(
  Control_128.10 = PathwayMerged$Control_128.10,
  Control_130 = PathwayMerged$Control_130,
  Control_128.13 = PathwayMerged$Control_128.13 
)
upsetPathwayData <- data.frame(pathwayVennInput)
ggVennDiagram::ggVennDiagram(pathwayVennInput, label_alpha = 0.5) +
  #ggVennDiagram::ggVennDiagram(upsetData, label_alpha = 0.5) + #Works too
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Venn Diagram of Top50 Pathways Between Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )
###Upset Plot===----
# install.packages("UpSetR")
library(UpSetR)
#UpSetR::upset(as.data.frame(upsetData), 
UpSetR::upset(fromList(listInput), 
      nsets = 3,
      #sets = c("Control_128.10", "Control_130", "Control_128.13"),
      number.angles = 0,
      mb.ratio = c(0.65, 0.35),
      point.size = 5,
      line.size = 1.3,
      order.by = "freq",
      #group.by = "freq",
      mainbar.y.label = "Number of Common DEG",
      sets.x.label = "Number of Top DEG",
      main.bar.color = "#6666FF",
      sets.bar.color = "#7DD305" , #flouroscent green
      matrix.color = "#FF6699", #pink #FF7D58", #orange
      matrix.dot.alpha = 0.1,
      shade.color = "#999",
      empty.intersections = "on",
      text.scale = c(1.8, 2, 1.2, 1, 2, 2),
      #text.scale = 1.5,
      keep.order = TRUE)

###ComplexHeatmap::Upset===----
m <- make_comb_mat(listInput)
col_size = comb_size(m)
# Draw the UpSet plot
ComplexHeatmap::UpSet(m,
      top_annotation = upset_top_annotation(m),
      right_annotation = upset_right_annotation(m),
      left_annotation = NULL,
      row_names_side = "left",
      column_title = "Plot")

#Using ComplexUpset::upset==
#install.packages("ComplexUpset")
library(ComplexUpset)
head(upsetData)
upsetDataFrame <- UpSetR::upset(fromList(listInput))
ComplexUpsetData <- upsetDataFrame$New_data
intersections <- colnames(upsetDataFrame$New_data)[1:3]
# genres <- colnames(upsetData)
ComplexUpset::upset(data = ComplexUpsetData, 
                    intersect = intersections, 
                    name='Intersecting Groups', 
                    width_ratio = 0.1,
                    guides = 'over') +
  ggtitle("Intersection of DEGs Between Groups")
# intersection_subset = c("Control_128.10", "Control_130", "Control_128.13")
# arrangedVenn =  arrange_venn(ComplexUpsetData, sets=intersection_subset)