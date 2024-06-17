# Downstream analysis using EdgeR Package (Loads limma Package)
# Reference: https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# EdgeR stores data in a data object called DGEList 

library(edgeR)
# Create groupData as a factor variable from the row names of colData
groupData <- factor(rownames(colData))
str(countData)

DGEListObject <- DGEList(counts = countData,
             genes = rownames(countData))#,
             #group = groupData)

#### Classic Mode (Not Recommended) ----
#To run edgeR in classic mode, we need to perform 3 steps: 
#..calculate normalization factors, estimate dispersions, 
#..and finally perform the exact test for differential expression.
z <- calcNormFactors(y)
z <- estimateDisp(z)

et <- exactTest(z)

#### Pairwise comparison (https://rnnh.github.io/bioinfo-notebook/docs/DE_analysis_edgeR_script.html) ----
#To run in the more popular glm mode, there are two testing methods..
#..Likelihood ratio tests (recommended for singlecellRNASeq & datasets with no replicates)
#..and quasi-likelihood F-tests (recommended for bulkRNAsq)

dim(DGEListObject)

#Check if columns names of countData are matching the order of columns of colData
summary(colnames(countData) == rownames(colData))

#Add grouping information to DGEList object
DGEListObject$samples$group <- as.factor(colData$drug)

#Comparison later can be between groups, usually the factors are arranged alphabetically
#..however we can change the order by relevel(). This only changes the order
levels(DGEListObject$samples$group) #Before arranging

DGEListObject$samples$group <- relevel(DGEListObject$samples$group, ref = "Control")
levels(DGEListObject$samples$group) #After Arranging Control is the first group


dim(DGEListObject)

#Filtering Lowly Expressed Genes
#Create Object to filter genes with low expression
counts.keep <- filterByExpr(DGEListObject)
summary(counts.keep) #Use this to select rows of DGEListObject to keep the TRUES

DGEListObject <- DGEListObject[counts.keep, , keep.lib.sizes = FALSE] 
# keep.lib.sizes = FALSE; library sizes (Number of RNASeq reads of the samples are recalculated after filtering)
dim(DGEListObject)

#May remove counts.keep as it is not needed later
rm(counts.keep)

#Normalization of library size between the samples to minimize bias towards highly expressed genes
#Before Normalization
DGEListObject$samples$norm.factors
# Add Normalization by assigning calcNormFactors function directly to the DGEListObject
DGEListObject <- calcNormFactors(DGEListObject)
#After Normalization (data added)
DGEListObject$samples$norm.factors

#Estimate common dispersion and tagwise(genewise) dispersion
#By estimating gene dispersion, we are estimating the relative variability of true expression levels between replicates.

#Can define a design matrix from colData
condition_ <- colData$drug #Will be used in the estimateDisp() command as the design= argument
#Again estimateDisp() function can be added directly to DGEListObject

DGEListObject <- estimateDisp(DGEListObject,
                              design = model.matrix(~condition_))
#This will add many information to the DGEListObject including design

### Pairwise Testing ----
condition_
#The elements of this will be used as pair argument for the exactTest() function
#The results of the exactTest() will be assigned to their own .DGEExact object 
control_128.10.DGEExact <- exactTest(DGEListObject,
                                     pair = c("Control", "128-10"))
#exactTest compares first two group variables if not provided with the pair argument

#Extract the most differentially expressed genes for pair condition using topTags()
control_128.10.topTags <- topTags(control_128.10.DGEExact)

#To see the most differentially expressed genes
control_128.10.topTags #This adds logFC, logCPM PValue etc

#### GLM Approach from vignette(GLM Approach: https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) ----
design <- model.matrix(~0+group, data = DGEListObject$samples)
#The 0+ in the model formula to not include an intercept column and instead to include a column for each group
#The column names have a "group" prefix; change it back as per the defined groups
colnames(design) <- levels(DGEListObject$samples$group)

design


#Any groups can be compared using the contrast argument of the glmQLFTest or glmLRT function

#y <- normLibSizes(y)

fit <- glmQLFit(DGEListObject, design)
# qlf <- glmQLFTest(fit, coef = 2)

qlf <- glmQLFTest(fit, contrast = c(-1, 1, 0))
#If there are A, B, C in the comparison, above will compare B to A
#Meaning of contrast is to make comparison -1*A + 1*B + 0*C; which is simply B-A

# The contrast vector can be constructed using makeContrasts also, if it is convenient
#The same above comparison could have been made by
BvsA <- makeContrasts(B-A, levels = design)
qlf <- glmQLFTest(fit, contrast=BvsA)

#Can make three pairwise contrasts by
my.contrasts <- makeContrasts(BvsA=B-A, CvsB=C-B, CvsA=C-A, levels = design)

qlf.BvsA <- glmQLFTest(fit, contrast = my.contrast[,"BvsA"])
topTags(qlf.BvsA)
qlf.CvsA <- gmlQLFTest(fit, contrast = my.contrast[,"CvsA"])
topTags(qlf.CvsA)
