#Downstream analysis using EdgeR Package (Loads limma Package)
#EdgeR stores data in a data object called DGEList 

library(edgeR)
# Create groupData as a factor variable from the row names of colData
groupData <- factor(rownames(colData))

y <- DGEList(counts = countData,
             group = groupData)

keepEdge <- filterByExpr(y, group = groupData)



design <- model.matrix(~ cellLine + drug)

keepEdge <- 
y <- y[keepEdge, , keep.lib.sizes=FALSE]

y <- normLibSizes(y)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)