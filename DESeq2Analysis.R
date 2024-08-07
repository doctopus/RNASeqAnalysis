# Downstream Analysis of RNASeq

##### Define Functions ----
setupProject <- function(project) {
  #Create 'project' dir if not same as name of *.Rproj dir as root/proj/io, else
  #..set root as proj dir, create io as root/io; set input_dir & output_dir vars
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
install_and_load_packages <- function(cran_packages, bioc_packages) {
  # Install missing CRAN packages
  new_packages_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
  if (length(new_packages_cran) > 0) {install.packages(new_packages_cran)}
  # Install missing Bioconductor packages
  new_packages_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
  if (length(new_packages_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
    BiocManager::install(new_packages_bioc, update = FALSE)
  }
  # Load all packages
  all_packages <- c(cran_packages, bioc_packages)
  sapply(all_packages, require, character.only = TRUE)
}
sentence_case <- function(name) { 
  # Sentence case first word if not uppercase with/out numbers/"-" (eg.DN-A1)
  # Split the sentence into words
  words <- unlist(strsplit(name, " "))
  # Check if the first word should be converted
  first_word <- words[1]
  if (!grepl("^[A-Z0-9]+$", first_word) && !grepl("-", first_word)) {
    # Convert the first word to sentence case
    words[1] <- paste0(toupper(substring(first_word, 1, 1)), 
                       tolower(substring(first_word, 2)))
  }
  # Join the words back into a sentence
  return(paste(words, collapse = " "))
}
##### Setup Project ----
## Initiate project
setupProject("RNASeqAnalysis") ; print(paste0("Working dir is: ", getwd()))
# If any project specific override: output folder 1.1_HierarchicalCategory_Top100DEGs
# output_dir <- paste0(output_dir, "/1.1_HierarchicalCategory_Top100DEGs")

## Install & Load Packages
cran_packages <- c("annotate", "circlize", "clipr", "devtools", "EnhancedVolcano", "ggpubr", "ggrepel", "matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "viridis")
bioc_packages <- c("apeglm", "clusterProfiler", "ComplexHeatmap", "DESeq2", "DOSE", "enrichplot", "genefilter", "GSVA", "org.Hs.eg.db", "pathview", "xCell")
install_and_load_packages(cran_packages, bioc_packages)

##### Source & Process Input files ----
file_samplesheet<-paste0(input_dir,"/samplesheet.csv")
file_fc<-paste0(input_dir,"/featureCounts_0")
## Process input files
sampleData<-read.csv(file = file_samplesheet, 
                  header=TRUE, 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE,
                  row.names = "sample")
sampleData <- sampleData[,c("cellLine", "drug")]
sampleData <- sampleData %>% mutate(drug = gsub("-", ".", drug))
sampleData[, c("cellLine", "drug")] <- lapply(sampleData[, c("cellLine", "drug")], factor)
head(sampleData)

fc <- read.delim(file_fc, row.names = NULL, check.names = FALSE)
genes_to_keep <- rowSums(fc[, 8:ncol(fc)]) > 1        #Keep Genes which are expressed in >1 sample
desired_order <- rownames(sampleData)                    #Order of columns per sampleData

fc1 <- fc %>%
  filter(genes_to_keep) %>%                           # Keep genes with expression in >1 sample
  mutate(EnsembleID = gsub("\\..*$", "", Geneid)) %>% # Add column of removed dot of EnsembleID
  filter(!is.na(EnsembleID) & EnsembleID != "") %>%   # Remove empty or NA EnsembleID rows
  dplyr::select(7:ncol(.)) %>%                        # Remove the initial columns except length
  arrange(desc(Length)) %>%                           # Order the rows based on gene length
  filter(!duplicated(EnsembleID)) %>%                 # Remove duplicated EnsembleID keep 1st
  dplyr::select(-Length) %>%                          # Remove the Length column
  column_to_rownames("EnsembleID")                    # Set EnsembleID as row names

# Order the columns according to desired_order
fc2 <- fc1[, desired_order]
head(fc2)

##### Make DDSObject (INPUT NEEDED- Define if Doing Subset Analysis) ----
### ############################### ###
#subsetToAnalyze <- NULL #If not performing subset analysis
subsetToAnalyze <- "318" # Define 318 or 358 cellLine Subset to analyze 
### ############################### ###
perform_subset_analysis <- !is.null(subsetToAnalyze) && subsetToAnalyze != ""
# Perform conditional branching
if (perform_subset_analysis) {
  countData <- fc2 %>% dplyr::select(contains(subsetToAnalyze))
  colData <- sampleData %>% filter(cellLine == subsetToAnalyze)
  design_formula <- ~ drug
} else {
  countData <- fc2
  colData <- sampleData
  design_formula <- ~ cellLine + drug
}
# Check validity of data format
if (all(colnames(countData) == rownames(colData))) {
  message("Data ready for DESeq2")
} else {
  stop("Column names of countData do not match row names of colData")
}
#### Create DESeqDataSet object----
ddsObject <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = design_formula)

#Keeping rows that have at least 10 reads for a minimum number of samples
#Minimal number of samples is the smallest group size, eg here 12 of each cellLine
#..or minimal number of samples for which non-zero counts would be considered interesting; 3 replicates
if (perform_subset_analysis) {
  smallestGroupSize <- 3
} else {
  smallestGroupSize <- 12
}
#counts(ddsObject)
#keep <-  rowSums(counts(dds2))>= 10
keep <-  rowSums(counts(ddsObject)>= 10) >= smallestGroupSize

ddsObject_filtered <- ddsObject[keep,]

ddsObject_filtered$drug

#Since we are primarily comparing between different drugs, so our primary level
#of comparison is drugs, and here reference level is Control 
#(this only reorders, since the default comparison is with first in the list)
ddsObject_filtered$drug <- relevel(ddsObject_filtered$drug, ref="Control")

ddsObject_filtered$drug <- droplevels(ddsObject_filtered$drug) #remove the levels (of drug) 
# ...which do not have samples in the current data set. Here nothing removed

##### Run DESeq2 Analysis ----
dds <- DESeq(ddsObject_filtered)

##### DDS Quality Control (Optional) ----
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

#To find how many differentially expressed genes based on a p-value cutoff
table(res$padj < 0.05)

#To see the normalized counts from the dds object
normCounts <- counts(dds, normalized = TRUE)

##To define contrast pairs in the results argument as follows
resControlx128.10 <- results(dds, contrast = c("drug", "128.10", "Control"))
head(resControlx128.10)

#Alternatively pick a combination from 
resultsNames(dds)
resControlx128.13 <- results(dds, name="drug_128.13_vs_Control")
head(resControlx128.13)


##Plot DESeq2 dispersion re-estimation procedure
plotDispEsts(dds)

# Visualize the distribution of p-values
hist(res$pvalue, breaks=20, col="grey")
hist(res$pvalue, breaks=0:50/50, xlab="p value", col="grey", 
     main="Histogram of nominal p values")
hist(res$padj, breaks=20, col="grey")

# Two common visualization for DE analysis are MA-plot and Volcano Plot
# MA-plot: Visualize relationship between a genes’ mean expression..
#..and its fold-change between experimental conditions
plotMA(res, ylim = c(-3, 3))
plotMA(res)
#plotMA(resOrdered)

#Volcano-plot: DESeq2 does not provide a function, but use base R
#..Here red are highlighted to show genes that are DE with Padj<0.05
plot(res$log2FoldChange, -log10(res$pvalue), 
     xlab="log2 Fold-change",
     ylab="-log P-value",
     pch=20,
     cex=0.5)
points(res$log2FoldChange[res$padj<0.05], -log10(res$pvalue[res$padj<0.05]),
       col="red", pch=20, cex=0.5)
abline(v=0, h=-log10(0.05), lty="dashed", col="grey")




# Since we want rank and visualization, shrink data by either apeglm, ashr or "normal" algorithm
# We provide the name or number of the coefficient we want to shrink,
#..where the number refers to the order of the coefficient in resultsNames(dds)
# coef in the lfcShrink should be from the output of resultsNames(dds)

# See the comparison groups
resultsNames(dds)

resLFCapeglm <- lfcShrink(dds, coef = "drug_130_vs_Control", type = "apeglm")
resLFCnormal <- lfcShrink(dds, coef = "drug_130_vs_Control", type = "normal")


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

##### DDS Apply Transformation ----
# Apply transformation & estimate dispersion trend
# vsd <- vst(dds, blind = FALSE) # VST: Variance Stabilizing Transformation
rld <- rlog(dds, blind=FALSE) # RLT: Regularized Log Transformation (Selected for this analysis)
# head(assay(vsd), 2)
head(assay(rld), 2)

### Heatmaps ----
#my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
my_colors <- colorRampPalette(c("white", "#E32600"))(99)

# Extract Transformed values
#rld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformationrld<-vst(dds)
# Using rld transformed data for the analysis [rld <- rlog(dds)]
length(rld)

## creating distance matrix (Not Needed)
sampleDists_subset <- as.matrix(dist(t(assay(rld))))
hm<-pheatmap::pheatmap(as.matrix(sampleDists_subset),
                       annotation_col = colData, 
                       col=my_colors,
                       main = paste("Distance Matrix", subsetToAnalyze, "Cell Line"),
                       annotation_legend=TRUE)

#Plot Heatmap of Hierarchical Clustering
rld_mat <- assay(rld) #Extract the transformed matrix from the object
rld_cor <- cor(rld_mat) #Compute pairwise correlation values
head(rld_cor)
heat.colors <- RColorBrewer::brewer.pal(6, "BrBG")
pheatmap(rld_cor, 
         annotation=colData, 
         color = heat.colors, border_color = NA,
         main = paste("Hierarchical Correlation", subsetToAnalyze, "Cell Line"),
         fontsize = 10, fontsize_row = 10, height = 20)


### PCA Plots ----

## PCA Plot using VSD Transformation 
# Argument in plotPCA ntop=length(vsd)) to include all genes; 
#...default is ntop=500, to plot top 500 feature variance
plotPCAvst <- plotPCA(vsd, intgroup = "drug", returnData=FALSE, ntop=length(vsd)) 
plotPCAvstData <- plotPCA(vsd, intgroup="drug", returnData=TRUE, ntop=length(vsd))
pcaVST <- plotPCAvst+geom_label_repel(data=plotPCAvstData, aes(label=name))+
  ggtitle(label=paste("PCA plot of", subsetToAnalyze, "Cells Line: VST Transformed Data"))+
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

## PCA Plot using RLT Transformation #Almost a similar plot
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


### DEG Prepare Data for Plot (INPUT NEEDED- Define Comparison Groups)----
ComparisonColumn <- "drug"
factor1 <- "130" #Choose from 128.10, 128.13 & 130
factor2 <- "Control"
title <- paste(factor1, "vs", factor2)  
# resultsNames(dds)
e <- as.character(c(ComparisonColumn, factor1, factor2))
#Shrink dds based on comparison data
res <- lfcShrink(dds, contrast = e, type = "normal") #Shrinked Result (L2FC, padj etc); "normal" algorithm

annotation_colors <- list(
  drug = c("128.10"="#9FD900", 
           "128.13"="#FAA800", 
           "130"="#ff5d8f", 
           "Control"="#6c757d"),
  
  cellLine =c("358"="#006E18", 
              "318"="#832161")
)
icolors <- colorRampPalette(c("blue",
                              "white",
                              "red"))(99)
#Combine shrunk results (with lfc, padj etc with normalized count data)
resdata <- merge(as.data.frame(res), #Get the counts data alongwith the comparison results
                 as.data.frame(assay(rld)), #To get transformed count values
                 #as.data.frame(counts(dds, normalized=FALSE)), #To get original count data
                 by = "row.names",
                 sort = FALSE)

names(resdata)[1] <- "EnsembleID" #Rename the first column (Row.names) as EnsembleID
#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db

#Add a column of Gene translated from EnsembleID
resdata$Gene<- mapIds(sp, keys=resdata$EnsembleID, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata$EntrezID<- mapIds(sp, keys=resdata$EnsembleID, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")

resdata <- resdata[order(resdata$padj),] #Order by ascending p-adj significance;low p =significant at top
resdata <- resdata[complete.cases(resdata$padj),] #Keep rows which have p-adj data
resdata <- resdata[complete.cases(resdata$Gene),] #Remove rows that don't have assigned genes
resdata <- resdata[!duplicated(resdata$Gene),] #Remove rows with same gene, keeping significant ones (low padj)
resdata <- resdata[order(resdata$log2FoldChange),] #Now order per log2FoldChange :Ascending

####SKIP STARTS if !(Need only list of Top and Bottom as separate lists)----
TopN <- 500
Criteria <- "Down" # "Up" or "Down"

if (Criteria == "Up") {
  TopX <- head(resdata, TopN)
} else if (Criteria == "Down") {
  TopX <- tail(resdata, TopN)
} else {
  stop("Invalid value for Criteria. Select either 'Up' or 'Down'")
}

# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, "Gene", "log2FoldChange", "padj", sep = "|"), names(TopX), value = TRUE)
TopX <- TopX[, columns_to_keep]
gene_names <- mapIds(org.Hs.eg.db,
                     keys = TopX$Gene,
                     column = "GENENAME",
                     keytype = "SYMBOL",
                     multiVals = "first")

TopX$GeneName <- gene_names[TopX$Gene]
TopX <- TopX %>% mutate(GeneName = sapply(GeneName, sentence_case))
DetailsColName <- paste(subsetToAnalyze, "CellLine", factor1, "vs",factor2,  Criteria, "Regulated", sep = "_")
TopX[, DetailsColName] <- paste(TopX$Gene, TopX$GeneName, sep = ": ")

##MANUALLY Save this data in spreadsheet
df_temp <- TopX %>% dplyr::select(DetailsColName)
# install.packages("clipr")
# library(clipr)
# Export the data frame to clipboard
write_clip(df_temp)

####SKIP ENDS if !(Need only list of Top and Bottom as separate lists)----
#### Find Top Bottom DEG to Plot----
TopN <- 50
Top50Up<-head(resdata, TopN)
Top50Down<-tail(resdata, TopN)
Top100Rows<-rbind(Top50Up, Top50Down)
R <- as.character(rownames(colData))
# R<-c(R, "Gene")
Top100<-Top100Rows %>% dplyr::select(all_of(R), "Gene") #Keep selected columns
Top100<-data.frame(Top100, check.names = FALSE)
rownames(Top100)<-Top100$Gene

#Keep only the compared columns
# Find the column names that contain either factor1 or factor2
columns_to_keep <- grep(paste(factor1, factor2, sep = "|"), names(Top100), value = TRUE)

# Subset the data frame to keep only those columns
Top100_f1vsf2 <- Top100[, columns_to_keep]

# Filter the colData to keep only those rows
colData100_f1vsf2 <- colData[columns_to_keep, ]

# Order the columns by 'drug'
ordered_indices <- order(colData100_f1vsf2$drug, decreasing = TRUE)
Top100_f1vsf2_ordered <- Top100_f1vsf2[, ordered_indices]
colData100_f1vsf2_ordered <- colData100_f1vsf2[ordered_indices, ]

# Filter annotation colors to include only relevant levels
annotation_colors_filtered <- list(
  drug = annotation_colors$drug[names(annotation_colors$drug) %in% unique(colData100_f1vsf2_ordered$drug)],
  cellLine = annotation_colors$cellLine[names(annotation_colors$cellLine) %in% unique(colData100_f1vsf2_ordered$cellLine)]
)

#### Skip if !(Need to show Gene Description in Heatmap) ----
gene_symbols <- rownames(Top100_f1vsf2_ordered)
gene_names <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "GENENAME",
                     keytype = "SYMBOL",
                     multiVals = "first")
#Copy the data frame
new_df <- Top100_f1vsf2_ordered
# Assign gene names to GeneName column
new_df$GeneName <- gene_names[rownames(Top100_f1vsf2_ordered)]
# Convert GeneName to sentence case
new_df <- new_df %>% mutate(GeneName = sapply(GeneName, sentence_case))
# Append GeneName(Symbol) to The GeneName column
new_df$GeneName <- paste(rownames(Top100_f1vsf2_ordered), new_df$GeneName, sep = ": ")
##MANUALLY Save this data in spreadsheet
new_df_control_130 <- new_df %>% select("GeneName")
# install.packages("clipr")
library(clipr)
# Export the data frame to clipboard
write_clip(new_df_control_130)

# Append GeneName to row names separated by ":"
rownames(new_df) <- paste(rownames(Top100_f1vsf2_ordered), new_df$GeneName, sep = ": ")

##Remove the GeneName column if no longer needed
new_df$GeneName <- NULL
# Assign it back as the Top100_f1vsf2_ordered
Top100_f1vsf2_ordered <- new_df


#### Plot Top 100 DEG ----
hmt100f1vsf2<-pheatmap::pheatmap(Top100_f1vsf2_ordered, scale="row", 
              annotation_col=colData100_f1vsf2_ordered,
              annotation=colData100_f1vsf2_ordered, 
              cluster_cols = F,
              cluster_rows = T,
              cellwidth=10,
              annotation_legend =TRUE, 
              color= icolors,
              annotation_colors = annotation_colors_filtered,
              fontsize_row = 8,
              main=paste0("Top 300 Fold Change: ", subsetToAnalyze, " cell Line \n", "[", title, "]"))
# saveFigure(figure=hmt50,fileName="Top50FoldChange_heatmap_Control_128.13",h=12,w=12)

#######Save Gene Data for Venn Diagram Analysis----

###==VENN: Alternate *Top & Bottom* Flow to Save *Top & Bottom DEG* Separately ----
# Run after every new set of Top50 and Bottom50; this appends lists to listInputz
# rm(listInputz) #Remove the previously existing listInputz object if not want to append
# Initialize the listInputz if it doesn't exist
if (!exists("listInputz")) {
  listInputz <- list()
}

# Generate the list names based on the provided variables
list_name_up <- paste0(subsetToAnalyze, "_", factor2, "_", factor1, "_Up")
list_name_down <- paste0(subsetToAnalyze, "_", factor2, "_", factor1, "_Down")

# Extracting the values from the "Gene" column of Top50Up and Top50Down
genes_up <- Top50Up$Gene
genes_down <- Top50Down$Gene

# Create the new lists to be appended
new_lists <- list(
  setNames(list(as.list(genes_up)), list_name_up),
  setNames(list(as.list(genes_down)), list_name_down)
)

# Append new lists to listInputz if they don't already exist
for (new_list in new_lists) {
  new_list_name <- names(new_list)
  if (!new_list_name %in% names(listInputz)) {
    listInputz <- c(listInputz, new_list)
  }
}

# Print listInputz to verify the contents
print(listInputz)
##RUN TO PLOT: When listInputz has all sets of required data----
#VENN DIAGRAM: listInputz as input
ggVennDiagram::ggVennDiagram(listInputz, label_alpha = 0.5)+
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle(paste0(subsetToAnalyze, " cell Line \n", "Up and Down DEG Between Groups")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )
#UPSET PLOT: listInputz as input 
library(UpSetR)
#UpSetR::upset(as.data.frame(upsetData), 
UpSetR::upset(fromList(listInputz), 
              nsets = 6,
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
              sets.bar.color = "#30C38B" , #flouroscent green
              matrix.color = "#007DEF", #pink #FF6699", #orange
              matrix.dot.alpha = 0.1,
              shade.color = "#999",
              empty.intersections = "on",
              text.scale = c(1.8, 2, 1.2, 1, 2, 2),
              #text.scale = 1.5,
              keep.order = TRUE)

##=========Experiment: Attempting to extract each group of from the upset plot
UpsetPlotData <- unlist(listInputz, use.names = TRUE)
head(UpsetPlotData)

UpsetPlotData <- UpsetPlotData[!duplicated(UpsetPlotData)]
intersection123z <- Reduce(intersect, listInputz)
print(paste("Intersection", toString(intersection123z)))

###=Alternate Flow Skip Ends-----

###==VENN: Alternate *All DEG* Flow to Save *All DEG* for Venn Diagram Analysis----
# Initialize the listInput if it doesn't exist
if (!exists("listInput")) {
  listInput <- list()
}

# Generate the list name based on the provided variables
list_name <- paste0(subsetToAnalyze, "_", factor2, "_", factor1)

# Extract the row names from Top100_f1vsf2_ordered
genes <- rownames(Top100_f1vsf2_ordered)

# Create the new list to be appended
new_list <- setNames(list(as.list(genes)), list_name)

# Append the new list to listInput if it doesn't already exist
if (!list_name %in% names(listInput)) {
  listInput <- c(listInput, new_list)
}

# Print listInput to verify the contents
print(listInput)

##RUN TO PLOT: When listInput has all sets of required data----
#VENN DIAGRAM: listInput as input
ggVennDiagram::ggVennDiagram(listInput, label_alpha = 0.5)+
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle(paste0(subsetToAnalyze, " cell Line \n", "DEG Between Groups")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )
#UPSET PLOT: listInput as input 
library(UpSetR)
#UpSetR::upset(as.data.frame(upsetData), 
UpSetR::upset(fromList(listInput), 
              nsets = 3, #Modify to include total number of lists
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
              sets.bar.color = "#30C38B" , #flouroscent green
              matrix.color = "#007DEF", #pink #FF6699", #orange
              matrix.dot.alpha = 0.1,
              shade.color = "#999",
              empty.intersections = "on",
              text.scale = c(1.8, 2, 1.2, 1, 2, 2),
              #text.scale = 1.5,
              keep.order = TRUE)
###=Alternate *All DEG* Flow Ends-----

######VENN:  OLD WAY OF Manual List creation and Venn diagram ----
#Create individual lists each time
row_names <- rownames(Top100_f1vsf2_ordered)
# T50Control_130 <- data.frame(matrix(ncol = 1, nrow = length(row_names)))
# T50Control_130 <- data.frame(ColumnNames = column_names)

#T50Control_128.10 <- data.frame(Control_128.10 = row_names)
#T50Control_128.13 <- data.frame(Control_128.13 = row_names)
T50Control_130 <- data.frame(Control_130 = row_names)
# head(T50Control_130)
#### OLD WAY Venn Diagram for DEGs ----
###Prepare data for Venn and Upset Plots
T50Merged <- cbind(T50Control_128.10, T50Control_130, T50Control_128.13)
listInput <- list(
  Control_128.10 = T50Merged$Control_128.10,
  Control_130 = T50Merged$Control_130,
  Control_128.13 = T50Merged$Control_128.13
)
upsetData <- data.frame(listInput)
#install.packages("ggVennDiagram")
library(ggVennDiagram)
# Create the Venn diagram #DONT CREATE LIST
# venn_data <- list(
#   "Control vs 128.10" = Control_128.10,
#   "Control vs 130" = Control_130,
#   "Control vs 128.13" = Control_128.13
# )

## Plot using ggVennDiagram===
ggVennDiagram::ggVennDiagram(listInput, label_alpha = 0.5) +
  #ggVennDiagram::ggVennDiagram(upsetData, label_alpha = 0.5) + #Works too
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle(paste(subsetToAnalyze, "cell Line \n", "Top100 DEG Between Groups")) +
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


##### Variable Genes ----
#A copy of resdata which is already based on a comparison pair
resdataf1vsf2 <- resdata
R <- as.character(rownames(colData))
resdataf1vsf2 <- resdataf1vsf2 %>% dplyr::select(all_of(R), "Gene")
rownames(resdataf1vsf2) <- resdataf1vsf2$Gene
# Remove the Gene column
resdataf1vsf2$Gene <- NULL

#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(resdataf1vsf2)),100)
mat <- resdataf1vsf2[topVarGenes, ]
mat <- mat - rowMeans(mat)

#=SKIP-START If !(Show Only Once Cell Line) Now pipeline has branched logic, so no need for this---
cellLineToPlot <- "358"
# Find columns with "318" in their names
DsKeepCols <- grep(cellLineToPlot, names(resdataf1vsf2), value = TRUE)
# Subset the data frame to keep only those columns
resdataf1vsf2 <- resdataf1vsf2[, DsKeepCols]
topVarGenes <- head(order(-genefilter::rowVars(resdataf1vsf2)),100)
mat <- resdataf1vsf2[topVarGenes, ]
mat <- mat - rowMeans(mat)
#=SKIP-ENDS---
## Plot Variable Genes--
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

##### Volcano Plot ----
volcanoPlot<-EnhancedVolcano(resdata,
              lab = resdata$Gene,
              title = paste(subsetToAnalyze, "cell Line \n", "[", title, "]"),
              x = 'log2FoldChange',
              y = 'pvalue',
              pCutoff = 10e-4, #Default P value cutoff is 10e-6
              FCcutoff = 1, #Default log2FC cutoff is >|2|
              pointSize = 3.0,
              labSize = 6.0,
              drawConnectors = TRUE)
print(volcanoPlot)

##### HALLMARK: Pathway Analysis With Hallmark Geneset ----
########## Input: List of descending log2FoldChange with names as Genes
human_hall_file<-paste0(input_dir,"/GSEA/h.all.v2023.2.Hs.symbols.gmt")
# human_hall_file<-paste0(input_dir,"/GSEA/h.all.v7.2.symbols.gmt")
genesets <- read.gmt(human_hall_file)

#Copy Result(p.adj, L2FC result) with normalized count data
resdata_gsea <- resdata
resdata_gsea <- resdata_gsea[!is.na(resdata_gsea$log2FoldChange),] #Remove rows without log2fc data
ordered_data <- resdata_gsea[order(-resdata_gsea$log2FoldChange),] #Order by decreasing log2fc

logFC <- ordered_data$log2FoldChange #Take the log2FC column only
names(logFC) <- ordered_data$Gene #Assign Row names as their corresponding Gene

#head(logFC)
#If need to plot a known list of genes (intersection123) NOT POSSIBLE SINCE LIST IS TOO SMALL
#intersectionDEGs <- intersection123 #Intersection Genes among Top DEG between comparison groups
#logFC_filtered <- logFC[names(logFC) %in% intersectionDEGs] #keep rows where gene names are in intersectionDEGs

#print(logFC_filtered[1:10])

gsea_hallmark <- GSEA(logFC, 
                      TERM2GENE=genesets, 
                      verbose=TRUE, 
                      pvalueCutoff=0.05)

gsea_hallmark@result$Description <- gsub('HALLMARK_','',gsea_hallmark@result$Description)
gsea_hallmark@result$Description <- gsub('_',' ',gsea_hallmark@result$Description)

#Need LogFC in descending order along with the gene names
options(enrichplot.colours = c("#e5383b","#007DEF"))
dotplot(gsea_hallmark,
        x = "NES",
        showCategory = 50,
        #title = paste("Differential Hallmark Pathways:", title),
        orderBy = "NES",
        #color = "p.adjust", # Map 'p.adjust' to color
        font.size = 10) +
  ggtitle(paste("Differential Hallmark Pathways:", subsetToAnalyze, "cell Line \n", "[", title, "]")) +
  theme(plot.title = element_text(hjust = 0.5))
#NES is Normalized Enrichment Score

## Save Pathway Data for Venn Diagram ---
pathways <- gsea_hallmark@result$Description 
#pathways <- rownames(Top50_f1vsf2_ordered)

#Path23Control_128.10 <- data.frame(Control_128.10 = pathways)
#Path23Control_128.13 <- data.frame(Control_128.13 = pathways)
Path23Control_130 <- data.frame(Control_130 = pathways)
##== HALLMARK VENN :TODO Improve Flow for Venn Diagram for Pathway Analysis----
##Pathway Venn Diagram
# PathwayMerged <- cbind(Path23Control_128.10, Path23Control_128.13, Path23Control_130)
# pathwayVennInput <- list(
#   Control_128.10 = PathwayMerged$Control_128.10,
#   Control_130 = PathwayMerged$Control_130,
#   Control_128.13 = PathwayMerged$Control_128.13 
# )

pathwayVennInput <- list(
  Control_128.10 = as.list(as.vector(Path23Control_128.10[,1])),
  Control_130 = as.list(as.vector(Path23Control_130[,1])),
  Control_128.13 = as.list(as.vector(Path23Control_128.13[,1]))
)
upsetPathwayData <- data.frame(pathwayVennInput)
ggVennDiagram::ggVennDiagram(pathwayVennInput, label_alpha = 0.5) +
  #ggVennDiagram::ggVennDiagram(upsetData, label_alpha = 0.5) + #Works too
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle(paste(subsetToAnalyze, "cell Line \n", "Common Significant Hallmark Pathways Between Groups")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )
intersectionPathways <- Reduce(intersect, pathwayVennInput)
print(paste("Intersection", toString(intersectionPathways)))

##### GO: Pathway Analysis with GO (https://learn.gencore.bio.nyu.edu/rna-seq-analysis/deseq-2/)----
########## Input: List of descending log2FoldChange with names as Genes
# library(clusterProfiler)
# library(enrichplot)
# library("org.Hs.eg.db")

resdata_go <- resdata
head(resdata_go)
dim(resdata_go)
#Remove NA rows
resdata_go <- resdata_go[!is.na(resdata_go$Gene),] #Remove rows without Gene data
ordered_data <- resdata_go[order(-resdata_go$log2FoldChange),] #Order by decreasing log2fc

original_gene_list <- ordered_data$log2FoldChange #Take the log2FC column only

#Name the vector
names(original_gene_list) <- ordered_data$Gene #Assign Row names as their corresponding Gene
# view(original_gene_list)
#May use ↑ this Original_gene_list in KEGG Cluster Network Plot

gene_list = na.omit(original_gene_list) #Omit any NA values (Previously Done)
head(gene_list)

gene_list = sort(gene_list, decreasing = TRUE) #Sort in decreasing order of log2FoldChange (Already Done)
head(gene_list)

##GO comprises three orthogonal ontologies
#MF: Molecular Function, 
#BP: Biological Process,
#CC: Cellular Component,
#ALL: All Components

onto = "MF"
dictionary <- c(ALL = "ALL", CC="Cellular Component", BP="Biological Process", MF="Molecular Function")

# dictionary[[onto]]

gse <- gseGO(geneList = gene_list,
         ont = onto, # Determines number of ontology terms per feature (BP, MF, CC, ALL) to be displayed
         keyType = "SYMBOL", #Source of annotation, check options available by keytypes(org.Hs.eg.db)
         nPerm = 10000, #Higher number of permutations will get better result at the cost of longer analysis
         minGSSize = 3,
         maxGSSize = 800,
         pvalueCutoff = 0.05,
         verbose = TRUE,
         OrgDb = org.Hs.eg.db,
         pAdjustMethod = "none")

# library(DOSE)
# options(enrichplot.colours = c("#e5383b","#007DEF")) #Set colors of Dotplot
#### Dot Plot
dotplot(gse, 
        font.size = 15,
        label_format = 35,# Width of the labels
        title=paste0("GO Gene Set Enrichment Analysis \n", subsetToAnalyze, " Cell Line - ",
                    dictionary[[onto]], " - ", "[", title, "]"),
        showCategory=5, 
        split=".sign") +
  facet_grid(.~.sign)+
  theme(plot.title = element_text(size=15, face = "bold"))#+
  #ggtitle(paste("Go Enrichment Analysis:", subsetToAnalyze, "Cell Line -",
                # dictionary[[onto]], "-", "[", title, "]"))

#### Enrichment Map Plot (Similarity Matrix)
# Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. 
# ..In this way, mutually overlapping gene sets tend to cluster together, making it easy to identify functional modules.
x2 = pairwise_termsim(gse)

emapplot(x2, 
         shadowtext= TRUE, 
         font.size = 15,
         showCategory = 20) +
  ggtitle(paste0("GO Enrichment Map of GSE Analysis \n", subsetToAnalyze, " Cell Line - ",
                dictionary[[onto]], " - ", "[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold"))


#### Category Net Plot 
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, 
         categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = 3)+
  ggtitle(paste0("GO Enrichment Gene-Concept Network \n", subsetToAnalyze, " Cell Line - ",
                dictionary[[onto]], " - ", "[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Plot linkages of genes and enriched concepts in GO categories")

#### Ridge Plot
#Grouped by gene set, density plots are generated by using the frequency of fold change values per gene within each set. 
#Helpful to interpret up/down-regulated pathways.

ridgeplot(gse,
          label_format = 35, # Width of the labels)
          showCategory = 10) +
  ggplot2::labs(x = "Enrichment Distribution")+
  ggtitle(paste0("GO Enrichment Distribution \n", subsetToAnalyze, " Cell Line - ",
                dictionary[[onto]], " - ", "[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Ridgeline plot for GSEA result of Top 10 Categories")


#### GSEA Plot (Per Selected Pathway)---
#GSEA PLOT Use the `Gene Set` param for the index in the title, and as the value for geneSetId

#Running Enrichment Score Plot
#Plot of the Running Enrichment Score (green line) for a gene set as the analysis walks down the ranked gene list,
#..including the location of the maximum enrichment score (the red line). 
#..The black lines in the Running Enrichment Score show where the members of the gene set appear in the ranked list of genes, indicating the leading edge subset.
#..The Ranked list metric shows the value of the ranking metric (log2 fold change) as you move down the list of ranked genes. 
#..The ranking metric measures a gene’s correlation with a phenotype.

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

keytypes(org.Hs.eg.db)

#View Table: Can view pathway enrichment results by extracting result matrix from the kk2 object
#The enriched pathways are stored in the ID and Description columns.
go_result = as.matrix(gse@result) 
head(go_result)
go_result = go_result[, 1:10]
head(go_result)


#### Pathview (Per Selected Pathway)---
# Not for GO Data Set
 
##### KEGG: Pathway Analysis with KEGG Gene Set Enrichment Analysis----
########## Input: List of descending log2FoldChange with names as EntrezID
#Convert gene IDs for gseKEGG function
#Not all gene ID will be converted
ids = bitr(names(original_gene_list), fromType = "SYMBOL", toType="ENTREZID",
           OrgDb = org.Hs.eg.db)

#Remove duplicate IDS (Using SYMBOL here, but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]), ]

head(dedup_ids)
#Create a new dataframe entrezData which has the respective ENTREZ IDs for the gene symbols.
colnames(dedup_ids) = c("Gene", "Entrez_ID")
entrezData = merge(resdata_go, dedup_ids, by = "Gene")
head(entrezData)


#Create a vector of the gene universe
kegg_gene_list = entrezData$log2FoldChange

#Name vector with ENTREZ IDs
names(kegg_gene_list) = entrezData$Entrez_ID

#Omit any NA Values
kegg_gene_list = na.omit(kegg_gene_list)

#Sort in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
# view(kegg_gene_list)


kegg_organism = "hsa"
kk2 = gseKEGG(geneList = kegg_gene_list,
              organism = kegg_organism, #For Humans
              nPerm = 10000, #Higher number of permutations for more accuracy at the cost of length of analysis
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              pAdjustMethod = "none",
              keyType = "ncbi-geneid") 

#### Dot Plot
dotplot(kk2, showCategory = 5, 
        label_format = 35, # Width of the labels)
        font.size = 15,
        title=paste0("KEGG Gene Set Enrichment Analysis \n", subsetToAnalyze, " Cell Line - ", "[", title, "]"),
        split=".sign") +
  facet_grid(.~.sign) +
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Gene Set Enrichment Analysis of KEGG")

#### Enrichment Map Plot (Similarity Matrix)
# Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. 
# ..In this way, mutually overlapping gene sets tend to cluster together, making it easy to identify functional modules.
k2 = pairwise_termsim(kk2)
emapplot(k2, 
         shadowtext= TRUE,
         font.size = 15,
         showCategory = 20)+
  ggtitle(paste0("KEGG Enrichment Map of GSE Analysis \n", subsetToAnalyze, " Cell Line - ", "[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold"))


#### Category Net Plot 
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, 
         categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = 3)+
  ggtitle(paste0("KEGG Enrichment Gene-Concept Network \n", subsetToAnalyze, " Cell Line - ","[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Plot linkages of genes and enriched concepts in KEGG categories")

head(kk2@geneList)

#### Category Net Plot -with Gene Name
#Translate Entrez ID to Gene Symbol 
kk3 <- setReadable(kk2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cnetplot(kk3, 
         categorySize="pvalue", 
         foldChange=gene_list, 
         showCategory = 3)+
  ggtitle(paste0("KEGG Enrichment Gene-Concept Network \n", subsetToAnalyze, " Cell Line - ","[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Plot linkages of genes and enriched concepts in KEGG categories")

#### Ridge Plot
#Grouped by gene set, density plots are generated by using the frequency of fold change values per gene within each set. 
#Helpful to interpret up/down-regulated pathways.

ridgeplot(kk2,
          label_format = 35, # Width of the labels
          showCategory = 10) +
  ggplot2::labs(x = "Enrichment Distribution")+
  ggtitle(paste0("KEGG Enrichment Distribution \n", subsetToAnalyze, " Cell Line - ", "[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Ridgeline plot for GSEA result of Top 10 Categories")


#### GSEA Plot (Per Selected Pathway)---
#GSEA PLOT Use the `Gene Set` param for the index in the title, and as the value for geneSetId

#Running Enrichment Score Plot
#Plot of the Running Enrichment Score (green line) for a gene set as the analysis walks down the ranked gene list,
#..including the location of the maximum enrichment score (the red line). 
#..The black lines in the Running Enrichment Score show where the members of the gene set appear in the ranked list of genes, indicating the leading edge subset.
#..The Ranked list metric shows the value of the ranking metric (log2 fold change) as you move down the list of ranked genes. 
#..The ranking metric measures a gene’s correlation with a phenotype.

gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

#View Table: Can view pathway enrichment results by extracting result matrix from the kk2 object
#The enriched pathways are stored in the ID and Description columns.
kegg_result = as.matrix(kk2@result) 
head(kegg_result)
kegg_result = kegg_result[, 1:10]
head(kegg_result)


#### Pathview (Per Selected Pathway)---
#This will create a PNG and different PDF of the enriched KEGG pathway and xml files.
#gene.data: This is kegg_gene_list created above
#pathway.id: The user needs to enter this. Enriched pathways + the pathway ID are provided in the gseKEGG output table (above).
#species: Same as organism above in gseKEGG, which we defined as kegg_organism

# BiocManager::install("pathview")
# library(pathview)
#Get pathway.id from kk2$ID

r1 =pathview(gene.data = kegg_gene_list, pathway.id = "hsa04721", #Ids from gseKEGG output (kk2)
  species = kegg_organism)
# ↑ This saves the png, preview.png and an xlm file in the project dir

#==== REACTOME: Gene Set Enrichment Analysis of Reactome Pathway (NOT WORKING)----
# (https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html)
BiocManager::install("ReactomePA")
BiocManager::install("reactome.db")
library("reactome.db")
remove.packages("ReactomePA")
library("ReactomePA")
data("geneList")
head(geneList)
head(kegg_gene_list)
dim(kegg_gene_list)
de <- names(kegg_gene_list)[abs(kegg_gene_list) > 1.2] #Pick the highly foldchanged ones
length(de)
head(de)
x <- enrichPathway(gene=de, pvalueCutoff=0.05, readable=T)
y <- gsePathway(geneList = kegg_gene_list,
                minGSSize=120,
                pvalueCutoff=0.05,
                pAdjustMethod="BH",
                verbose=TRUE)
enrichmap(y)

# ↑ Error: ReactomePA depends on reactome.db, which is not installing


# Trying ReactomeGSA package
# (https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/using-reactomegsa.html)
BiocManager::install("ReactomeGSA")
library(ReactomeGSA)

available_methods <- get_reactome_methods(print_methods = FALSE, return_result = TRUE)

available_methods$name

reactome_result <- analyse_sc_clusters(kegg_gene_list, verbose = TRUE)

##=== GO Enrichment Analysis --Experimenting============
# (https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html)-
# Enrichment Analysis done for a set of genes (Here selecting the top log2FoldChange genes)
#Input: List of Genes with name as EntrezIDs

# Example: Use the example data set included with the package DOSE
# data(geneList, package="DOSE")
# head(geneList, 10)

#SKIP START--
# Set fold change > 1.2 as being DE genes Select a few genes (Directly get the set with Gene names)
gene <- names(kegg_gene_list)[abs(kegg_gene_list)>1.2]
length(gene) #10 Genes Selected

#Gene Name needed so can use the gene_list from GO analysis and skip next step
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

#gene_list is list with Gene Names
head(gene_list) #list with Gene Names
head(kegg_gene_list) #list with EntrezID and o
names(kegg_gene_list)

# SKIP ENDS --

#GO comprises three orthogonal ontologies
#MF: Molecular Function, 
#BP: Biological Process,
#CC: Cellular Component

# ggo <- groupGO(gene = names(kegg_gene_list),
#              OrgDb = org.Hs.eg.db,
#              ont = "BP",
#              level = 3,
#              readable = TRUE)
# 
# head(ggo)
# goplot(ggo)

#Select a few genes of interest (Here with high log2FoldChange)
gene <- names(kegg_gene_list)[abs(kegg_gene_list)>1]
length(gene) #10 Genes Selected

ego <- enrichGO(gene = gene,
                universe = names(kegg_gene_list), #
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", #Only Works for BP; Not for ALL/MP/CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable=TRUE)

head(ego)

#Induced GO DAG(Directed Acyclic Graph) of Significant Terms
goplot(ego)+
       # showCategory =50,#No effect
       # color = "p.adjust", #No effect
       # layout = "sugiyama", #No effect
       # geom = "text")+ #No effect
  ggtitle(paste0("GO Induced DAG Plot of Significant Terms \n", subsetToAnalyze, " Cell Line - ",
                dictionary["BP"], "-", "[", title, "]"))+
  theme(plot.title = element_text(size=15, face = "bold")) +
  labs(subtitle = "",
       caption = "Induced GO DAG (Directed Acyclic Graph) of Significant Terms")



# median.FC.values <- runif(n = dim( ego@result)[1] , min=-15, max = 15)
# names(median.FC.values) <- rownames(ego@result)
# 
# ego@result$medianFC <-  median.FC.values
# 
# head(ego@result)
# head(as.data.frame(ego))
# 
# options(enrichplot.colors = c("pink", "blue"))
# dotplot(ego, color="p.adjust") +
#   scale_colour_gradient2(low="green", mid="white", high="red",
#                          limits = c(-15, 15),
#                          breaks = c(-15, -7.5, 0, 7.5, 15),
#                          labels = c("-15down", "-7.5", "0", "7.5", "15up"),  
#                          guide=guide_colorbar(reverse=TRUE) )  +
#   labs(size="Count", colour="Median logFC")

############################################ --

###### Extra Goodies: Not Part of main code: May be removed ↓----

#### Get Gene Name from Gene Symbol ----
genelist <- data.frame(GeneSymbol = unlist(intersection123))
gene_symbols <- unique(genelist$GeneSymbol)
gene_names <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "GENENAME",
                     keytype = "SYMBOL",
                     multiVals = "first")

# Create a new column for gene names
genelist$GeneName <- gene_names[genelist$GeneSymbol]
#Convert to sentence case (Function defined at the top)
genelist <- genelist %>% mutate(GeneName = sapply(GeneName, sentence_case))

#### Export to clipboard----
install.packages("clipr")
library(clipr)
# Export the data frame to clipboard
write_clip(genelist)

#### Print table as image----
install.packages("gt")
library(gt)
gt_table <- gt(genelist) %>% tab_header (title = paste("Common DEG in ",subsetToAnalyze, "Cell Line"))
gtsave(gt_table, paste0(subsetToAnalyze, " Common DEG List 358", ".png"), path= input_dir)


#### Upset Plot===----
#Upset plot.. May be removed since it is moved to new code.
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

#### ComplexHeatmap::Upset===----
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




###=Alternative way to plot DEG Heatmap: Complex Heatmap Code----
## Complex Heat Map Code--
# library(ComplexHeatmap)
# library(circlize)

# col_fun <- colorRamp2(c(min(Top25), median(Top25), max(Top25)), c("blue", "white", "red"))
# Adjust color mapping to ensure proper visualization
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
title <- paste("Top 100 Fold Change:", factor2, "vs", factor1)

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
