# Define the project name, organize folder structure. Script in project root.
project <- "RNASeq-Kinase"
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
###
# Project specific output for V0.3 UnTrimmedNewIndex
output_dir <- "/Users/i/Dropbox/Clinic3.0/RStudio/RNASeq-Kinase/output/v0.3_UnTrimmedNewIndex"
###Select and Merge FC File----
# source_dir <- "W:/Shared/LRI/Labs/yuj2lab/Tapas"
# source_dir <- "/mnt/beegfs/beherat2/kinase_RNASeq/data"
#fc_dir <- paste0(input_dir, "/subread/featureCounts/")
# fc_dir <- paste0(source_dir, "/FeatureCountsIndex/")
#fc_dir <- paste0("../../../Documents/Clinic3.0/kinaseRNASeq/FeatureCounts/")

# Get list of bam files
file_list <- paste0(input_dir, "/bamFiles.list.txt")
p<- scan(file = file_list, what = "character")
samples=gsub(".Aligned.sortedByCoord.out.bam","",basename(p))
list<-as.numeric(c(1:length(samples)))

#Create a dataframe "assigned" to know which featureCounts have maximum assigned values
fcs<-c("featureCounts_0","featureCounts_s1","featureCounts_s2")
assigned<-data.frame(matrix(nrow=3,ncol=2))
colnames(assigned)<-c("file","assigned")
assigned$assigned<-0

for(f_file in 1:length(fcs)){
  assigned$file[f_file]<-fcs[f_file]
  for(i in 1:length(list)){
    f<-paste(fcs[f_file],".",i,".summary",sep='')
    x<-read.table(file = paste0(fc_dir,f),header = T, sep="\t",fill = T, stringsAsFactors = FALSE)
    
    assigned$assigned[f_file]<-assigned$assigned[f_file] + x[,2][which(x[,1]=="Assigned")]
  }
}

#Select which of the assigned dataframe has max assigned values
f_file<-which.max(assigned$assigned)

#### create matrix for the chosen counts file
f<-paste(fcs[f_file],".",list[1],sep='')
Mat<-read.table(file = paste0(fc_dir,f),header = T, sep="\t",fill = T, stringsAsFactors = FALSE)
Idx<-c(1:dim(Mat)[1]) #Total number of rows
Mat<-data.frame(Idx,Mat) #Create a column of row numbers called Idx
colnames(Mat)[8]<-samples[1] #Change the column name of counts to its bam file

for(i in 2:length(list)){
  f<-paste(fcs[f_file],".",i,sep='')
  tmp<-read.table(file = paste0(fc_dir,f),header = T, sep="\t",fill = T, stringsAsFactors = FALSE)
  tmp<-tmp[,c(1,7)]
  colnames(tmp)[2]<-samples[i]
  
  Mat<-merge(Mat,tmp,by.x = "Geneid",by.y = "Geneid",all.x = TRUE,all.y = TRUE)
}

#Order the rows in Mat file based on initial order Idx
Mat<-Mat[order(Mat$Idx),]
#Create the folder featureCountsCombined Else point to the featureCounts folder
#The combined featureCount folder as featureCounts_0 without any extension
write.table(Mat,file = file.path(input_dir, fcs[f_file]),append = FALSE,quote = FALSE,sep="\t",row.names=FALSE,col.names = TRUE)




### create DESeq Object ====
# BiocManager::install("DESeq2")
# BiocManager::install(c("DESeq2", "annotate", "org.Hs.eg.db", "AnnotationDbi"))
# install.packages(c("devtools"))
#packages1 <- c("tidyverse", "devtools", "annotate", "org.Hs.eg.db", "DESeq2")
#invisible(lapply(packages1, function(x) require(x, character.only = T, quietly = T)))
sapply(c("tidyverse", "devtools", "annotate", "org.Hs.eg.db", "DESeq2"), require, character.only = TRUE)
# suppressPackageStartupMessages({
#   library(tidyverse) #Loads dplyr
#   #library(dplyr)
#   library(devtools)
#   library(annotate) #Loads AnnotationDbi
#   library(org.Hs.eg.db)
#   #library(AnnotationDbi)
#   library(DESeq2)
# })

meta<-data.frame(Sample=c("318-Control_A1", "318-Control_A2", "318-Control_A3", 
                          "318-128-10_B1", "318-128-10_B2", "318-128-10_B3", 
                          "318-130_C1", "318-130_C2", "318-130_C3", 
                          "318-128-13_D1", "318-128-13_D2", "318-128-13_D3", 
                          "358-Control_E1", "358-Control_E2", "358-Control_E3", 
                          "358-128-10_F1", "358-128-10_F2", "358-128-10_F3", 
                          "358-130_G1", "358-130_G2", "358-130_G3", 
                          "358-128-13_H1", "358-128-13_H2", "358-128-13_H3"),
                 Group=c("318-Control", "318-Control", "318-Control", 
                            "318-128-10", "318-128-10", "318-128-10", 
                            "318-130", "318-130", "318-130", 
                            "318-128-13", "318-128-13", "318-128-13", 
                            "358-Control", "358-Control", "358-Control", 
                            "358-128-10", "358-128-10", "358-128-10", 
                            "358-130", "358-130", "358-130", 
                            "358-128-13","358-128-13", "358-128-13"))
print("metadata created")

#Read combined featureCounts File
fc_file<-paste0(input_dir,"/featureCounts_0")
fc<-read.delim(fc_file,row.names=1,check.names = FALSE)
#Keep rows where more than one genes are expressed
genes_to_keep <- rowSums(fc[,7:ncol(fc)]) > 1
fc_1<-fc[genes_to_keep,]
fc_1<-fc_1[,6:ncol(fc_1)] #Remove the initial columns
fc_1$EnsembleID<-gsub("\\..*$","",rownames(fc_1)) #New column with removed suffix of Ensemble ID with the dot
fc_1<-fc_1[order(fc_1[,2],decreasing = TRUE),] # Order the rows
fc_1<-fc_1[!duplicated(fc_1$EnsembleID),] #Remove duplicated rows based on EnsembleID keeping the first occurrence

#Get Bioconductor Annotation Database
sp <- org.Hs.eg.db

#TPM Normalization Transcripts per million
geneLengths <- as_vector(subset(fc_1, select = c(Length)))
rpk <- apply( subset(fc_1, select = c(-Length,-EnsembleID)), 2, #2 implies the function to be applied column-wise
              function(x) x/(geneLengths/1000)) #Reads per kilobase by dividing the reads by respective length(in Kb)
#normalize by sample size using rpk values
tpm <- as.data.frame(apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)) #Divide each value in rpk with sum of all values in that column(in Million)
colSums(tpm) #Sum of each column of tpm values for each sample, to be used for qc and normalization purpose
tpm$GeneID<-row.names(tpm) #Row names that is the original Ensemble ID as GeneID column
tpm$EnsemblID<-gsub("\\..*$", "", rownames(tpm)) #EnsembleID as stripped out EmsebleID without dot and the number
#Get Matching data from Bioconductor
tpm$GeneSymbol<- mapIds(sp, keys=tpm$EnsemblID, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
tpm$EntrezID<- mapIds(sp, keys=tpm$EnsemblID, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
rownames(tpm)<-tpm$EnsemblID #Change rownames to stripped out EnsembleID

# Write the tpm dataframe to a CSV file named tpm.csv in the output folder
write.csv(tpm, file = file.path(output_dir, "tpm.csv"))

#Count Matrix for DESeq; needs a count file only with GeneID without the rpk and tpm calculations
#fc_1$Length<-NULL
# Create a new data frame without the "Length" column
#fc_1_new <- subset(fc_1, select = -Length) 
fc_1_new <- fc_1 %>% mutate(Length = NULL) #Without removing the Length Column
## get entrz ID and Gene symbol
fc_1_new$GeneSymbol<-mapIds(sp,keys=fc_1_new$EnsembleID, column = c("SYMBOL"), keytype = c("ENSEMBL"),multiVals = "first")

saveRDS(fc_1_new,file = file.path(output_dir, "clean_counts_for_DESeq.rds"))
rownames(fc_1_new)<-fc_1_new$EnsembleID
write.csv(fc_1_new, file = file.path(output_dir, "raw_counts_matrix_clean.csv"))

#Create DESeq Object
#library(DESeq2)
comparison<-"Group" #Name of the independent variable
deseq_formula<-as.formula(paste("~",comparison,collapse="")) # ~ separates the Group(dependent) variable from the explanatory (independent) variable 
# collapse="" is not necessary but used above since paste0 adds a space between the pasted values
#fc_2<-fc_1[,1:6]
fc_2<-fc_1_new[,meta$Sample]
dds<-DESeq(DESeqDataSetFromMatrix(countData=fc_2,colData = meta,design = deseq_formula))
saveRDS(dds,file = file.path(output_dir, "DESeqObject.rds"))

#Normalize count data
dds_normalized<-as.data.frame(counts(dds, normalized=TRUE))
dds_normalized<-cbind(dds_normalized,fc_1_new[,c("EnsembleID","GeneSymbol")])
write.csv(dds_normalized, file = file.path(output_dir, "DESeqNormalizedDataset.csv"))

### save the normalized DESeq object as well and use that for further analysis
saveRDS(dds_normalized, file = file.path(output_dir, "DESeqNormalizedObject.rds"))

# help(package="DESeq2", help="html")
# vignette("DESeq2")

#DownStream Analysis----
# BiocManager::install(c("clusterProfiler", "EnhancedVolcano", "GSVA", "DOSE"))
# install.packages(c("devtools"))
# devtools::install_github('dviraran/xCell')
# BiocManager::install("genefilter")
# library(genefilter)
sapply(c("genefilter", "clusterProfiler", "EnhancedVolcano", "GSVA", "DOSE", "RColorBrewer", "pheatmap", "xCell"), require, character.only = TRUE)

# packages2 <- c("clusterProfiler", "EnhancedVolcano", "GSVA", "DOSE", "RColorBrewer", "pheatmap", "xCell")
# invisible(lapply(packages2, function(x) require(x, character.only = T, quietly = T)))

# suppressPackageStartupMessages({
#   library(c("clusterProfiler",
#             "EnhancedVolcano", 
#             "GSVA", 
#             "DOSE",
#     "RColorBrewer", 
#     "pheatmap", 
#             "xCell"))
#   })

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

deseqobj<-paste0(output_dir,"/DESeqObject.rds")
#meta defined above 
rownames(meta) <- meta$Sample #Change Rownames as Sample (instead of default serial numbers)
my_colors <- colorRampPalette(c("blue", "white", "red"))(99)
### change these values according to metadata and comparisons required
ComparisonColumn<-"Group"
factor1<-"318-Control"
factor2<-"318-128-10"
factor3<-"318-130"
factor4<-"318-128-13"
factor1<-"358-Control"
factor2<-"358-128-10"
factor3<-"358-130"
factor4<-"358-128-13"

## Use DESeq object
dds<-readRDS(deseqobj)
rld<-vst(dds) #estimate dispersion trend and apply a variance stabilizing transformation

## creating distance matrix
sampleDists_subset <- as.matrix(dist(t(assay(rld))))
hm<-pheatmap::pheatmap(as.matrix(sampleDists_subset),annotation_col = meta, col=my_colors,annotation_legend=TRUE)
saveFigure(figure=hm,fileName="DistanceMatrix",h=7,w=7)

## creating PCA plot
pca<-plotPCA(rld, intgroup="Group")
pca<-pca + geom_text(aes(label=name),vjust=2, size = 3)
saveFigure(figure=pca,fileName="PCAPlot",h=10,w=20)

#If needs to label samples using ggrepel
library(ggrepel)
zz <- plotPCA(rld, intgroup="Group", returnData=TRUE)
pca_repel <- pca+geom_label_repel(data=zz, aes(label=name))
saveFigure(figure=pca_repel, fileName = "PCAPlot_Repel", h=15, w=20)

### plotting DEG
e<-as.character(c(ComparisonColumn, factor1,factor2))
res<-lfcShrink(dds,contrast=e,type="normal") ## or use res<-results(dds,contrast=e)
#lfcShrrink is to shrink log2 fold changes
resdata_subset <- merge(as.data.frame(res), as.data.frame(assay(rld)), by="row.names", sort=FALSE)
write.csv(resdata_subset, paste0(output_dir,"/","DifferentialExpressionAnalysis", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset)[1] <- "Gene"
resdata_subset <- resdata_subset[order(resdata_subset$padj),]
resdata_subset <- resdata_subset[!is.na(resdata_subset$padj),]
#resdata_subset$EnsembleID <- gsub("\\..*$", "", resdata_subset$Gene)
resdata_subset$GeneSymbol<- mapIds(sp, keys=resdata_subset$Gene, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset$EntrezID<- mapIds(sp, keys=resdata_subset$Gene, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
write.csv(resdata_subset, paste0(output_dir,"/","DifferentialExpressionAnalysis_cleaned", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

resdatagenes_subset <- resdata_subset[complete.cases(resdata_subset$GeneSymbol),]
resdatagenes_subset <- resdatagenes_subset[!duplicated(resdatagenes_subset$GeneSymbol), ]

resdatagenes_subset <- resdatagenes_subset[order(resdatagenes_subset$log2FoldChange),]
Top50genesdown_subset<-head(resdatagenes_subset, 50)
Top50genesup_subset<-tail(resdatagenes_subset, 50)
Top_subset<-rbind(Top50genesdown_subset, Top50genesup_subset)
#subset the counts
L<-as.character(meta$Sample)
L<-c(L, "GeneSymbol")
Top_subset<-Top_subset%>% dplyr::select(all_of(L))
TopD_subset<-data.frame(Top_subset, check.names = FALSE)
#name the rows as genesymbols
rownames(TopD_subset)<-Top_subset$GeneSymbol

## top 50 fold changes
hmt50<-pheatmap::pheatmap(TopD_subset[1:(length(TopD_subset)-1)],scale="row", annotation_col=meta,
                       annotation_legend =TRUE, color= colorRampPalette(c("blue","white","red"))(99), fontsize_row = 8, main="Top50FoldChange")
saveFigure(figure=hmt50,fileName="Top50FoldChange_heatmap",h=12,w=12)






## variable genes
Ds2<-resdatagenes_subset%>% dplyr::select(all_of(L), "GeneSymbol")
Ds2<-Ds2[!duplicated(Ds2$GeneSymbol), ]
rownames(Ds2)<-Ds2$GeneSymbol
Ds2<-Ds2%>% dplyr::select(-"GeneSymbol")
#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds2)),100)
mat <- Ds2[topVarGenes, ]
mat <- mat - rowMeans(mat)

#plot the variable genes in heatmap
hmVariable<-pheatmap::pheatmap(mat,color= colorRampPalette(c("blue","white","red"))(99), annotation_col=meta,
                       annotation_legend =TRUE, scale="row", fontsize_row = 8, show_rownames=T, main="Top100variablegenes")
saveFigure(figure=hmVariable,fileName="Top100VariableGenes",h=12,w=12)


## volcano plot
vp<-EnhancedVolcano(resdata_subset,
                    lab = resdata_subset$GeneSymbol,x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 10e-4,
                    FCcutoff = 1)
saveFigure(figure=vp,fileName="VolcanoPlot",h=8,w=8)



### pathway analysis with hallmark geneset
### create a resources directory 
human_hall_file <- paste0(input_dir, "/GSEA/h.all.v2023.2.Hs.symbols.gmt")
hall<-read.gmt(human_hall_file)


resdatagenes_gsea <- resdatagenes_subset
resdatagenes_gsea <- resdatagenes_gsea[!is.na(resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea <- resdatagenes_gsea[order(-resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea$FC_pval <- (resdatagenes_gsea$log2FoldChange)
logFC.l2n <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$FC_pval
names(logFC.l2n) <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$GeneSymbol

gsea.hall.l2n <- GSEA(logFC.l2n, TERM2GENE=hall, verbose=FALSE, pvalueCutoff=1)
gsea.hall.l2n.df <- as.data.frame(gsea.hall.l2n@result)
write.csv(gsea.hall.l2n.df, paste0(output_dir,"/","GSEA_output", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

dp<-dotplot(gsea.hall.l2n, x="NES", showCategory=50, orderBy= "NES",font.size = 7) 
saveFigure(figure=dp,fileName="HallmarkPathwayAnalysis")



#==============================================
#==============================================
## Can be run from counts. Dont need DESeq object for this. Cannot be run on cell lines
counts<-paste0(output_dir,"/tpm.csv")
Df<-read.csv(counts, check.names = FALSE)
Df<-Df[!duplicated(Df$GeneSymbol),]
Df<-Df[!is.na(Df$GeneSymbol),]
rownames(Df)<-Df$GeneSymbol
Df<-Df[,meta$Sample]

### xCell

xOut<-xCellAnalysis(Df)
xOut<-data.frame(xOut, check.names = FALSE)

write.csv(xOut, paste0(output_dir, "/XCellScores.csv"),row.names = TRUE)

xOut$avg<-rowSums(xOut)/ncol(xOut)
xOut<-xOut[which(xOut$avg != xOut$`318-Control_A1`),]
xOut$avg<-NULL
rownames(meta)<-meta$Sample

hmxCell<-pheatmap::pheatmap(xOut,
                       annotation = meta,
                       show_rownames=T,
                       cluster_rows = T, 
                       cluster_cols = T, 
                       scale="row",
                       color=colorRampPalette(c("blue", "white", "red"))(99))
saveFigure(figure=hmxCell,fileName="xCell_heatmap")



