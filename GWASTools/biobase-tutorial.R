# ExpressionSet:
  # assayData = expression data from microarray
  # phenoData = metadata on samples
  # featureData = metadata about features on the microarray
  # annotation = metadata about technology used in experiment (alternative to featureData?)
  # protocolData = processing protocol information, usually from manufacturer files
  # experimentData = flexible to describe the experiment

library(Biobase)
# you may use other Bioconductor packages to read array data, in .CEL files or similar,
# into an ExpressionSet object

# 4 Building an ExpressionSet from scratch ---------
# 4.1 arrayData = F (features) x S (samples) matrix 
dataDirectory <- system.file("extdata", package="Biobase")
exprsFile <- file.path(dataDirectory, "exprsData.txt")
exprs <- as.matrix(read.table(exprsFile, header=T, sep="\t", row.names=1, as.is=T))
head(exprs)
class(exprs)
dim(exprs)
colnames(exprs)
minimalSet <- ExpressionSet(assayData=exprs) # minimal set without any metadata
print(minimalSet)
# 4.2 phenoData = S x V (covariates) table
pDataFile <- file.path(dataDirectory, "pData.txt")
pData <- read.table(pDataFile, row.names=1, header=T, sep="\t")
head(pData)
dim(pData)
rownames(pData)
summary(pData)
all(rownames(pData) == colnames(exprs)) # rows of phenoData = cols of arrayData, including names
sapply(pData, class) # find the class of each column of pData
# make a data frame to store information about the phenoData covariates
metadata <- data.frame(labelDescription=
                           c("Patient gender",
                             "Case/control status",
                             "Tumor progress on XYZ scale"),
                         row.names=c("gender", "type", "score"))
# labelDescription column is required
# now make the phenoData into an AnnotatedDataFrame
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)
phenoData
# some AnnotatedDataFrame functions:
pData(phenoData)
phenoData[c("A", "Z"), "gender"]
pData(phenoData[phenoData$score > 0.8, ])
# 4.3 now add other annotations and feature data
# because there are only a certain number of standard microarray chips, there are pre-packaged
# metadata objects for popular chips/instruments
# there are also metadata packages for GO and KEGG information
annotation <- "hgu95av2" # specify the name of a particular Affy chip
# you can also add custom metadata for features to the AnnotatedDataFrame, for example flagging 
# features of interest to your analysis
# 4.4 Experiment description using a MIAME object
# slots of MIAME: name, lab, contact, title, abstract, url, samples, hybridization, normControls, 
# preprocessing, pubMedIDs, other
experimentData <- new("MIAME", name="Pierre Fermat", lab="Francis Galton Lab", contact="pfermat@lab",
                      title="Smoking-Cancer Experiment", abstract="An example expressionSet",
                      url="lab.com", other=list(notes="Created from text files"))
print(experimentData)
# 4.5 finally assembling the expressionSet
exSet <- ExpressionSet(assayData = exprs,
                       phenoData = phenoData, 
                       experimentData = experimentData, 
                       annotation = "hgu95av2")
print(exSet)

# 5 ExpressionSet Basics ------------------
help("ExpressionSet-class")
exSet
exSet$gender
featureNames(exSet)[1:5]
sampleNames(exSet)[1:5]
varLabels(exSet)
mat <- exprs(exSet)
dim(mat)
# subset using [features, samples]
exSet[1:5, 1:3]
exprs(exSet[1:5, 1:3])
