my.path <- "~/Documents/gwastools-tutorial/"
setwd(my.path)
# netCDF format exploration ---------------------
summary.ncdf <-
  function( control){
    for ( i in 1:control$nvars ){
      
      vname = control$var[[i]]$name    # variable name
      ndims = control$var[[i]]$ndims   # number of dimensions for this variable
      
      dimstring = paste(vname,'( variable',as.character(i),') has shape')
      for (j in 1:ndims) {
        dimstring <- paste(dimstring, as.character(control$var[[i]]$dim[[j]]$len))
      }
      
      cat(dimstring, fill=TRUE)
    }
  }

library(ncdf)
# read in an example netCDF 
ex.nc <- open.ncdf("example.nc")
summary.ncdf(ex.nc)
print(ex.nc)
get.var.ncdf(ex.nc, "SN") # get one of the variables, e.g. for plotting

# make a netCDF object in R
data(volcano)
z <- 10 * volcano
x <- 100 * (1:nrow(z))
y <- 100 * (1:ncol(z))

# define the coordinate variables
dim1 <- dim.def.ncdf("EW", "meters", as.double(x))
dim2 <- dim.def.ncdf("SN", "meters", as.double(y))

# define empty elevation variable
varz <- var.def.ncdf("Elevation", "meters", list(dim1, dim2), -1, longname="The Classic R New Zealand Volcano")

# associate the variable with a netCDF file, put variable in the file, and close
nc.ex <- create.ncdf("volcano.nc", varz)
put.var.ncdf(nc.ex, varz, z)
close.ncdf(nc.ex)

##### GWAS DATA CLEANING tutorial --------------------
library(GWASTools)
library(GWASdata)

# 2.2 SNP Annotation object -----
# Load SNP annotation df
data(illumina_snp_annot)
head(illumina_snp_annot)
# make this into a SnpAnnotationDataFrame object (derived from AnnotatedDataFrame in Biobase)
library(Biobase)
help("AnnotatedDataFrame-class")
# SNP annotation object requires at least snpID, chromosome, and position arrays
snpAnnot <- SnpAnnotationDataFrame(illumina_snp_annot)
print(snpAnnot)
varLabels(snpAnnot)
varMetadata(snpAnnot)
head(pData(snpAnnot)) # look at the data
# add the metadata on the variables
meta <- varMetadata(snpAnnot)
meta[c("snpID", "chromosome", "position", "rsID", "alleleA", "alleleB",
        "BeadSetID", "IntensityOnly", "tAA", "tAB", "tBB", "rAA", "rAB", "rBB"),
       "labelDescription"] <- c("unique integer ID for SNPs",
                                   paste("integer code for chromosome: 1:22=autosomes,",
                                              "23=X, 24=pseudoautosomal, 25=Y, 26=Mitochondrial, 27=Unknown"),
                                    "base pair position on chromosome (build 36)",
                                       "RS identifier",
                                       "alelleA", "alleleB",
                                       "BeadSet ID from Illumina",
                                       "1=no genotypes were attempted for this assay",
                                       "mean theta for AA cluster",
                                       "mean theta for AB cluster",
                                       "mean theta for BB cluster",
                                       "mean R for AA cluster",
                                       "mean R for AB cluster",
                                       "mean R for BB cluster")
varMetadata(snpAnnot) <- meta
varMetadata(snpAnnot)
# get the data in several ways
snpAnnot$snpID
getSnpID(snpAnnot)

snpAnnot[["chromosome"]]
getChromosome(snpAnnot)
table(getChromosome(snpAnnot))

getPosition(snpAnnot)
getVariable(snpAnnot, "rsID")

# 2.3 Scan Annotation object -----
# required variables: scanID, sex
data(illumina_scan_annot)
scanAnnot <- ScanAnnotationDataFrame(illumina_scan_annot)
print(scanAnnot)
varLabels(scanAnnot)
head(pData(scanAnnot))
# Add metadata to describe the columns
meta <- varMetadata(scanAnnot)
meta[c("scanID", "subjectID", "family", "father", "mother",
       "CoriellID", "race", "sex", "status", "genoRunID", "plate",
       "batch", "file"), "labelDescription"] <-
  c("unique ID for scans",
    "subject identifier (may have multiple scans)",
    "family identifier",
    "father identifier as subjectID",
    "mother identifier as subjectID",
    "Coriell subject identifier",
    "HapMap population group",
    "sex coded as M=male and F=female",
    "simulated case/control status" ,
    "genotyping instance identifier",
    "plate containing samples processed together for genotyping chemistry",
    "simulated genotyping batch",
    "raw data file")
varMetadata(scanAnnot) <- meta
varMetadata(scanAnnot)

# 2.4 Creating GDS and netCDF data files -----
path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
geno.file <- "tmp.geno.gds" #paste0(my.path, "tmp.geno.gds")
# only using first 3 samples
scan.annotation <- illumina_scan_annot[1:3, c("scanID", "genoRunID", "file")]
names(scan.annotation)[2] <- "scanName"

snp.annotation <- illumina_snp_annot[,c("snpID", "rsID", "chromosome", "position")]
# indicate which column of SNP annotation is referenced in data files
names(snp.annotation)[2] <-  "snpName"
head(snp.annotation)

# Making Genotype files -----
col.nums <- as.integer(c(1,2,12,13))
names(col.nums) <- c("snp", "sample", "a1", "a2")
col.nums
# works fine with either different setwd or a different path in filename
# Uknown issue left the file handle open last time
diag.geno.file <- "diag.geno.RData" #paste0(my.path, "diag.geno.RData")
diag.geno <- createDataFile(path = path,
                            filename = geno.file,
                            file.type = "gds",
                            variables = "genotype",
                            snp.annotation = snp.annotation,
                            scan.annotation = scan.annotation,
                            sep.type = ",",
                            skip.num = 11,
                            col.total = 21,
                            col.nums = col.nums,
                            scan.name.in.file = 1,
                            diagnostics.filename = diag.geno.file,
                            verbose = FALSE)
# Look at the values included in the "diag.geno" object which holds
#   all output from the function call
names(diag.geno)
# `read.file' is a vector indicating whether (1) or not (0) each file
#   specified in the `files' argument was read successfully
table(diag.geno$read.file)
# `row.num' is a vector of the number of rows read from each file
table(diag.geno$row.num)
# `sample.match' is a vector indicating whether (1) or not (0)
#   the sample name inside the raw text file matches that in the
#   sample annotation data.frame
table(diag.geno$sample.match)
# `snp.chk' is a vector indicating whether (1) or not (0)
#   the raw text file has the expected set of SNP names
table(diag.geno$snp.chk)
# `chk' is a vector indicating whether (1) or not (0) all previous
#   checks were successful and the data were written to the data file
table(diag.geno$chk)

# Run the function checkGenotypeFile to check that the GDS file 
# contains the same data as the raw data files.
check.geno.file <- "check.geno.RData"
check.geno <- checkGenotypeFile(path = path,
                                filename = geno.file,
                                file.type = "gds",
                                snp.annotation = snp.annotation,
                                scan.annotation = scan.annotation,
                                sep.type = ",",
                                skip.num = 11,
                                col.total = 21,
                                col.nums = col.nums,
                                scan.name.in.file = 1,
                                check.scan.index = 1:3,
                                n.scans.loaded = 3,
                                diagnostics.filename = check.geno.file,
                                verbose = FALSE)
# Look at the values included in the "check.geno" object which holds
#   all output from the function call
names(check.geno)
# 'geno.chk' is a vector indicating whether (1) or not (0) the genotypes
#   match the text file
table(check.geno$geno.chk)

# Now that you've created a GDS file, you can create a GdsGenotypeReader object, which you
# can use for analyses involving the whole dataset
gds <- GdsGenotypeReader(geno.file)
nscan(gds)
nsnp(gds)
head(getScanID(gds))
head(getSnpID(gds))
head(getChromosome(gds))
head(getPosition(gds))
# genotypes for the first 3 samples and  the first 5 SNPs
getGenotype(gds, snp=c(1,5), scan=c(1,3))
# and when you're done with the reader, you need to close the object, just like a file handle
close(gds)
# GdsGenotypeReader requires snp.id, snp.chromosome, snp.position, sample.id, and genotype variables
# you could also just make a more general GdsReader object to find out which variables the file has
gds <- GdsReader(geno.file)
getVariableNames(gds)
close(gds)

# Making Intensity files -----
qxy.file <- "tmp.qxy.gds"
col.nums <- as.integer(c(1,2,5,16,17))
names(col.nums) <- c("snp", "sample", "quality", "X", "Y")
diag.qxy.file <- "diag.qxy.RData"
diag.qxy <- createDataFile(path = path,
                           filename = qxy.file,
                           file.type = "gds",
                           variables = c("quality","X","Y"),
                           snp.annotation = snp.annotation,
                           scan.annotation = scan.annotation,
                           sep.type = ",",
                           skip.num = 11,
                           col.total = 21,
                           col.nums = col.nums,
                           scan.name.in.file = 1,
                           diagnostics.filename = diag.qxy.file,
                           verbose = FALSE)
diag.qxy
# now check to make sure the file matches the individual data and the annotation data
check.qxy.file <- "check.qxy.RData"
check.qxy <- checkIntensityFile(path = path,
                                filename = qxy.file,
                                file.type = "gds",
                                snp.annotation = snp.annotation,
                                scan.annotation = scan.annotation,
                                sep.type = ",",
                                skip.num = 11,
                                col.total = 21,
                                col.nums = col.nums,
                                scan.name.in.file = 1,
                                check.scan.index = 1:3,
                                n.scans.loaded = 3,
                                diagnostics.filename = check.qxy.file,
                                verbose = FALSE)
check.qxy
# And now you can use the GdsIntensityReader class
gds <- GdsIntensityReader(qxy.file)
nscan(gds)
nsnp(gds)
head(getScanID(gds))
head(getSnpID(gds))
head(getChromosome(gds))
head(getPosition(gds))
# quality score for the first 3 samples and  the first 5 SNPs
getQuality(gds, snp=c(1,5), scan=c(1,3))
# X intensity for the first 3 samples and  the first 5 SNPs
getX(gds, snp=c(1,5), scan=c(1,3))
# Y intensity for the first 3 samples and  the first 5 SNPs
getY(gds, snp=c(1,5), scan=c(1,3))
close(gds)
# GdsIntensityReader requires these variables: snp, chromosome, position, sampleID and
# one or more of (quality, X, Y, BAlleleFreq, LogRRatio)

# B Allele Freq. and LogRR files -----
baflrr.file <- "tmp.bl.gds"

col.nums <- as.integer(c(1,2,20,21))
names(col.nums) <- c("snp", "sample", "BAlleleFreq", "LogRRatio")

diag.baflrr.file <- "diag.baflrr.RData"
diag.baflrr <- createDataFile(path = path,
                          filename = baflrr.file,
                          file.type = "gds",
                          variables = c("BAlleleFreq","LogRRatio"),
                          snp.annotation = snp.annotation,
                          scan.annotation = scan.annotation,
                          sep.type = ",",
                          skip.num = 11,
                          col.total = 21,
                          col.nums = col.nums,
                          scan.name.in.file = 1,
                          diagnostics.filename = diag.baflrr.file,
                          verbose = FALSE)
# and again, we can use GdsIntensityReader to read this file
gds <- GdsIntensityReader(baflrr.file)
getBAlleleFreq(gds, snp=c(1,5), scan=c(1,3))
getLogRRatio(gds, snp=c(1,5), scan=c(1,3))
close(gds)

file.remove(geno.file, qxy.file, baflrr.file)
file.remove(diag.geno.file, diag.qxy.file, diag.baflrr.file)
file.remove(check.geno.file, check.qxy.file)

# 2.5 Combine data files with SNP and scan annotations -----
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(genofile)
# only the GDS file
genoData <- GenotypeData(gds)
# add the annotations
genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
# show off the functions of the GenotypeData object
nsnp(genoData)
nscan(genoData)
range(getScanID(genoData))
hasSex(genoData)
table(getSex(genoData))
hasScanVariable(genoData, "subjectID")
head(getScanVariable(genoData, "subjectID"))
getScanVariableNames(genoData)
range(getSnpID(genoData))
table(getChromosome(genoData, char=TRUE))
table(getChromosome(genoData, char=FALSE))
head(getPosition(genoData))
hasSnpVariable(genoData, "rsID")
head(getSnpVariable(genoData, "rsID"))
getSnpVariableNames(genoData)
getGenotype(genoData, snp=c(1,5), scan=c(1,5))
close(genoData)
# IntensityData objects are similar, but take a GdsIntensityReader object
# here's the Quality, X, Y data 
qxyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
gds <- GdsIntensityReader(qxyfile)
qxyData <- IntensityData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
getQuality(qxyData, snp=c(1,5), scan=c(1,5))
getX(qxyData, snp=c(1,5), scan=c(1,5))
getY(qxyData, snp=c(1,5), scan=c(1,5))
close(qxyData)
# here's the BAF, LRR data
baflrrfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
gds <- GdsIntensityReader(baflrrfile)
baflrrData <- IntensityData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
getBAlleleFreq(baflrrData, snp=c(1,5), scan=c(1,5))
getX(baflrrData, snp=c(1,5), scan=c(1,5)) # Gives an error because X isn't a variable found in this file
getLogRRatio(baflrrData, snp=c(1,5), scan=c(1,5))
close(baflrrData)

# 3 Batch Quailty Checks --------------------
# open the GDS file and create a GenotypeData object
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gds <- GdsGenotypeReader(gdsfile)
# sex is required for this function, so we need the scan annotation
genoData <-  GenotypeData(gds, scanAnnot=scanAnnot)
# 3.1 Missing call rates for samples and SNPs -----
## missing.n1: per SNP missing rate over all samples
# calculate missing.n1 separately for each sex
miss <- missingGenotypeBySnpSex(genoData)
names(miss)
head(miss$missing.counts)
miss$scans.per.sex
head(miss$missing.fraction)
allequal(snpAnnot$snpID, as.numeric(names(miss$missing.fraction))) # check that order matches
snpAnnot$missing.n1 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples",
  "except that females are excluded for Y chr SNPs")
summary(snpAnnot$missing.n1)
# plot missingness per SNP to look for outliers
hist(snpAnnot$missing.n1, ylim=c(0, 100), 
     col = "gray50", 
     xlab="SNP missing call rate",
     main="Missing Call Rate for All Probes")
# number of SNPs missing for every sample
length(snpAnnot$missing.n1[snpAnnot$missing.n1 == 1])
# fraction of autosomal SNPs with missing < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome < 23]
length(x[x < 0.05]) / length(x)
# same fraction for X chromosome
x <- snpAnnot$missing.n1[snpAnnot$chromosome == 23]
length(x[x < 0.05]) / length(x)
# same fraction for Y chromosome
x <- snpAnnot$missing.n1[snpAnnot$chromosome == 25]
length(x[x < 0.05]) / length(x)
## missing.e1: per sample missing rate over all SNPs, with bad SNPs (missing.n1 == 1) excluded
# find SNPs with 100% missing rate to exclude
sum(snpAnnot$missing.n1 == 1)
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1]
length(snpexcl)
# missing call rate per sample with snps excluded
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
names(miss)
head(miss$missing.counts)
head(miss$snps.per.chr)
# check on the number of excluded SNPs
sum(miss$snps.per.chr)
nrow(snpAnnot) - sum(miss$snps.per.chr)
head(miss$missing.fraction)
# check that the ordering is correct
allequal(names(miss$missing.fraction), scanAnnot$scanID)
scanAnnot$missing.e1 <- miss$missing.fraction
["missing.e1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n1<1",
  "except that Y chr SNPs are excluded for females")
summary(scanAnnot$missing.e1)
summary(scanAnnot$missing.e1[scanAnnot$sex == "M"]) # missing.e1 in males
summary(scanAnnot$missing.e1[scanAnnot$sex == "F"]) # missing.e1 in females
sum(scanAnnot$missing.e1 > 0.05) # Number of samples missing in > 5% of SNPs
# this number of bad samples is 0 because the HapMap data is high quality, but we'll pretend there were some bad ones anyway
auto <- colnames(miss$missing.counts) %in% 1:22 # autosome columns
missa <- rowSums(miss$missing.counts[ , auto]) / sum(miss$snps.per.chr[auto]) # missing rate in autosomes
summary(missa)
missx <- miss$missing.counts[ , "X"] / miss$snps.per.chr["X"] # missing rate in X chromosome
summary(missx)
# check that they match scan annotation
allequal(names(missa), scanAnnot$scanID)
allequal(names(missx), scanAnnot$scanID)
# and add these to the scan annotation
scanAnnot$miss.e1.auto <- missa
scanAnnot$miss.e1.xchr <- missx
# order scanAnnot by missing.e1 so duplicate subjectIDs with a higher missing rate are 
# marked as duplicates
scanAnnot <- scanAnnot[order(scanAnnot$subjectID, scanAnnot$missing.e1), ]
scanAnnot$duplicated <- duplicated(scanAnnot$subjectID) # duplicated() marks all repeats, but doesn't mark the first occurrence
table(scanAnnot$duplicated)
# VERY IMPORTANT! return scanAnnot to scanID order
scanAnnot <- scanAnnot[order(scanAnnot$scanID), ]
allequal(scanAnnot$scanID, sort(scanAnnot$scanID)) # a check on the fixed ordering
varMetadata(scanAnnot)["duplicated", "labelDescription"] <-
  "TRUE for duplicate scan with higher missing.e1"

## missing.n2: per SNP missing rate over all samples, with missing.e1 < 0.05 excluded
# find scanIDs for samples with missing.e1 > 0.05
scan.exclude <- scanAnnot$scanID[scanAnnot$missing.e1 > 0.05]
miss <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)
snpAnnot$missing.n2 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples with missing.e1<0.05",
  "except that females are excluded for Y chr SNPs")
summary(snpAnnot$missing.n2)

## missing.e2: per sample missing rate over all SNPs, with missing.n2 < 0.05 excluded
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n2 >= 0.05] # get snpIDs to exclude
length(snpexcl)
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
# now add these missing rates to the scan Annotation
scanAnnot$missing.e2 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n2<0.05",
  "except that Y chr SNPs are excluded for females")
summary(scanAnnot$missing.e2)
# and plot it in another histogram
hist(scanAnnot$missing.e2, col="gray50", 
     xlab="Fraction of missing calls over all probes with missing call rate < 0.05",
     main="Histogram of Sample Missing Call Rate for All Samples")

# 3.2 Calculate missing call rates by batch -----
# In the example HapMap data, each sample was in a different batch, as controls for another experiment
# so we're using a simulated batch variable
varLabels(scanAnnot)
getVariable(scanAnnot, "plate")
getVariable(scanAnnot, "batch")
length(unique(scanAnnot$batch))
table(scanAnnot$batch, useNA="ifany")
# Plot distribution of samples per batch
barplot(table(scanAnnot$batch), xlab="Batch", 
        ylab="Number of Samples", 
        main = "Distribution of Samples per Batch")
# now look at mean missing call rate per batch for all snps
batches <- unique(scanAnnot$batch)
bmiss <- rep(NA, length(batches))
names(bmiss) <- batches
bn <- rep(NA, length(batches))
names(bn) <- batches
for (i in 1:length(batches))
{
  missing.in.batch <- scanAnnot$missing.e1[is.element(scanAnnot$batch, batches[i])]
  bmiss[i] <- mean(missing.in.batch)
  bn[i] <- length(missing.in.batch)
}
# Now we can do an ANOVA to get a regression line of mean missing in batch vs. samples in batch
# in cases of higher batch numbers, we could use this to find outlier batches
# you can also do this separately for autosomes and X chromosome (use imagination)
batch.lm <- lm(bmiss ~ bn)
anova(batch.lm)
plot(bn, bmiss, pch=20,
     xlab = "Number of samples per batch", 
     ylab = "Mean missing call rate")
abline(batch.lm$coefficients)

# 3.3 Chi-Square test of allele frequency differences in batches -----
# chi-square test for differences in allele frequency between each batch and all other pooled batches
batch.chi <- batchChisqTest(genoData, batchVar="batch", return.by.snp=TRUE)
close(genoData)
dim(batch.chi$chisq)
batch.chi$lambda # genomic inflation factor
batch.chi$mean.chisq # mean chi-square value for each of the batches
# now test for association between batches and population groups using chi-square 
pop.table <- table(scanAnnot$race, useNA="ifany")
pop.table[1] / sum(pop.table)
pop.table[2] / sum(pop.table)
pop.table <- table (scanAnnot$race, scanAnnot$batch)
pop.batch.chi <- chisq.test(pop.table) # test for "ethnic effects"
pop.batch.chi$p.value
# fraction of CEU in each batch
ceu <- rep(NA, length(batches))
names(ceu) <- sort(batches)
for (i in 1:length(batches))
{
  this.batch <- scanAnnot$race[is.element(scanAnnot$batch, batches[i])]
  this.batch.ceu <- length(this.batch[this.batch == "CEU"])
  ceu[i] <- this.batch.ceu / length(this.batch)
}
allequal(names(ceu), names(batch.chi$mean.chisq)) # check that the batches match
# look for batches that differ from other batches of a similar ethnic composition (indicates possible genotyping artifact)
# plot mean chi-square vs fraction of CEU samples in batch
plot(ceu, batch.chi$mean.chisq, pch = 20,
     xlab = "Fraction of CEU samples per batch",
     ylab = "Mean Chi-square value",
     main = "Fraction of CEU samples per batch vs. Mean chi-square value")
abline(v=mean(ceu), lty = 2, col="red", lwd = 1.5)
# this plot is not interpretable because it's not a large enough dataset

# 4 Sample quality checks ----------------------
# 4.1 sample genotype quality scores -----
# calculate mean and median quality scores per sample
# needs IntensityData (with quality scores) and GenotypeData (with missingness information)
# open the data files
qxyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
qualGDS <- GdsIntensityReader(qxyfile)
qualData <- IntensityData(qualGDS, scanAnnot=scanAnnot)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genoGDS <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(genoGDS, scanAnnot=scanAnnot)
# look for scans with low quality
qual.results <- qualityScoreByScan(qualData, genoData)
close(qualData)
head(qual.results)
# plot the hist of median quality scores
hist(qual.results[ , "median.quality"], main = "Median Genotype Quality Scores of Samples", 
     xlab = "Median Quality", col = "gray50", xlim = c(0, 1.0))
# all of these are of very high quality

# 4.2 B allele frequency variance analysis -----
blfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
blGDS <- GdsIntensityReader(blfile)
blData <- IntensityData(blGDS, scanAnnot=scanAnnot)
nbins <- rep(12, 3)
slidingBAF12 <- sdByScanChromWindow(blData, genoData, nbins=nbins)
names(slidingBAF12)
dim(slidingBAF12[["21"]])
sds.chr <- meanSdByChromWindow(slidingBAF12, scanAnnot$sex)
sds.chr[["21"]]
sds.chr[["X"]]
# find windows with sample-chromosome pairs that have very high sd(BAF) compared to same window in other samples
res12bin4sd <- findBAFvariance(sds.chr, slidingBAF12, scanAnnot$sex, sd.threshold=4)
table(res12bin4sd[ , "chromosome"])
# plot BAF of all SNPs in identified chromosome-sample pairs vs. position
scanID <- as.integer(res12bin4sd[, "scanID"])
chrom <- as.integer(res12bin4sd[, "chromosome"])
chrom[res12bin4sd[, "chromosome"] == "X"] <- 23
bincode <- paste("Bin", res12bin4sd[, "bin"], sep = " ")
chromIntensityPlot(blData, scanID, chrom, info=bincode, ideogram=FALSE)
close(blGDS)

# 4.3 Missingness and heterozygosity within samples -----
# looking for:
  # samples with relatively high heterozygosity for all chromosomes (possible contaminated samples)
  # high missingness in specific chromosomal regions
miss <- missingGenotypeByScanChrom(genoData)
miss.rate <- t(apply(miss$missing.counts, 1, function(x) {
  x / miss$snps.per.chr}))
miss.rate <- as.data.frame(miss.rate)

cols <- names(miss.rate) %in% c(1:22, "X", "XY")
boxplot(miss.rate[,cols], main="Missingness by Chromosome",
        ylab="Proportion Missing", xlab="Chromosome")

boxplot(miss.rate$X ~ scanAnnot$sex,
        main="X Chromosome Missingness by Sex",
        ylab="Proportion Missing")

# Calculate heterozygosity by scan by chromosome
het.results <- hetByScanChrom(genoData)
close(genoData)
# Ensure heterozygosity results are ordered correctly
allequal(scanAnnot$scanID, rownames(het.results))
# Write autosomal and X chr heterozygosity to sample annot
scanAnnot$het.A <- het.results[,"A"]
scanAnnot$het.X <- het.results[,"X"]
varMetadata(scanAnnot)["het.A", "labelDescription"] <-
  "fraction of heterozygotes for autosomal SNPs"
varMetadata(scanAnnot)["het.X", "labelDescription"] <-
  "fraction of heterozygotes for X chromosome SNPs"

boxplot(scanAnnot$het.A ~ scanAnnot$race,
        main="Autosomal Heterozygosity")

female <- scanAnnot$sex == "F"
boxplot(scanAnnot$het.X[female] ~ scanAnnot$race[female],
        main="X Chromosome Heterozygosity in Females")

# 5 Sample Identity Checks ---------------------
# 5.1 Misannotated sex check -----
qxyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
intenGDS <- GdsIntensityReader(qxyfile)
inten.by.chrom <- meanIntensityByScanChrom(intenGDS)
names(inten.by.chrom)
close(intenGDS)
# use mean intensities by sample to identify misannotated sexes or sex chr aneuploidies
mninten <- inten.by.chrom[[1]] # mean intensities
dim(mninten)
allequal(scanAnnot$scanID, rownames(mninten)) # check that sample order is correct
# assign colors to each sex
xcol <- rep(NA, nrow(scanAnnot))
xcol[scanAnnot$sex == "M"] <- "blue"
xcol[scanAnnot$sex == "F"] <- "green"
nx <- sum(snpAnnot$chromosome == 23)
ny <- sum(snpAnnot$chromosome == 25)
# plot 1: All intensities
x1 <- mninten[ , "X"]
y1 <- mninten[ , "Y"]
main1 <- "Mean X vs \nMean Y Chromosome Intensity"
# plot 2: X Heterozygosity vs. X intensity
x2 <- mninten[ , "X"]
y2 <- scanAnnot$het.X
main2 <- "Mean X chromosome intensity vs\nmean X chromosome heterozygosity"
# plot 3: Heterozygosity on X vs Y intensity
y3 <- mninten[ , "Y"]
x3 <- scanAnnot$het.X
main3 <- "Mean X chromosome heterozygosity vs\nMean Y chromosome intensity"
# plot 4: X vs Autosomal heterozygosity
x4 <- scanAnnot$het.A[scanAnnot$sex == "F"]
y4 <- scanAnnot$het.X[scanAnnot$sex == "F"]
main4 <- "Mean autosomal heterozygosity vs\nmean X chromosome heterozygosity"
cols <- c("blue", "green")
mf <- c("male", "female")
xintenlab <- paste("X intensity (n=", nx, ")", sep="")
yintenlab <- paste("Y intensity (n=", ny, ")", sep="")
# plotting!
pdf("DataCleaning-sex.pdf")
par(mfrow=c(2,2))
plot(x1, y1, xlab=xintenlab, ylab=yintenlab,
     main=main1, col=xcol, cex.main=0.8)
legend("topright",mf,col=cols,pch=c(1,1))
plot(x2, y2, col=xcol, xlab=xintenlab,
     ylab="X heterozygosity", main=main2, cex.main=0.8)
plot(x3, y3, col=xcol, ylab=yintenlab,
     xlab="X heterozygosity", main=main3, cex.main=0.8)
plot(x4,y4, col="green", xlab="Autosomal heterozygosity",
     ylab="X heterozygosity", main=main4, cex.main=0.8)
dev.off()

# 5.2 Relatedness and IBD -----
library(SNPRelate)
gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
gdsobj <- snpgdsOpen(gdsfile)
ibdobj <- snpgdsIBDKING(gdsobj)
snpgdsClose(gdsobj)
names(ibdobj)
dim(ibdobj$kinship)
ibdobj$kinship[1:5, 1:5]
# get the pedigree from the scanAnnotation data and check it for internal consistency errors
# "The function pedigreeCheck finds any of a number of possible errors and inconsistencies within 
# pedigree data. If no problems are encountered, the output is NULL. If problems are encountered, 
# output contains information for the errors encountered"
ped <- pData(scanAnnot)[ , c("family", "subjectID", "father", "mother", "sex")]
dim(ped)
head(ped)
names(ped) <- c("family", "individ", "father", "mother", "sex") # rename columns for compatibility
head(ped)
(chk <- pedigreeCheck(ped)) # you'll get an error message 
dups <- chk$duplicates # first try removing duplicates
uni.ped <- pedigreeDeleteDuplicates(ped, dups)
(chk <- pedigreeCheck(uni.ped)) # now try adding a row for a parent who is missing an individual entry
ni <- chk$parent.no.individ.entry
parent <- data.frame(family=ni$family, individ=ni$parentID,
                     father=0, mother=0, sex="F",
                     stringsAsFactors=FALSE)
ped.complete <- rbind(uni.ped, parent)
(chk <- pedigreeCheck(ped.complete)) # I'm not sure what exactly the identified subfamilies mean
# "There are multiple subfamilies identified, so we will need to assign new family IDs to the 
# subfamilies. One subfamily has two unrelated people (likely founders), so we remove this 
# family from the pedigree."
ped.complete <- ped.complete[ped.complete$family != 58,]
subf <- chk$subfamilies.ident
table(subf$family)
subf.ids <- subf$individ[subf$subfamily == 2]
newfam <- ped.complete$individ %in% subf.ids
ped.complete$family[newfam] <- paste0(ped.complete$family[newfam], "-2")
table(ped.complete$family)
pedigreeCheck(ped.complete) # and now the pedigreeCheck passes
# "Now from the verified sample list excluding duplicates, we can calculate the expected
# relationships among the samples..."
rels <- pedigreePairwiseRelatedness(ped.complete)
length(rels$inbred.fam)
relprs <- rels$relativeprs
relprs[1:5, ]
table(relprs$relation)
# prepare to plot IBD estimates color coded by expected relationships from pedigree
samp <- pData(scanAnnot)[ , c("scanID", "subjectID")] # start a dataframe
samp <- samp[match(ibdobj$sample.id, samp$scanID), ] # find rows for ibdobj corresponding to rows in new df
names(samp) <- c("scanID", "Individ")
ibd <- snpgdsIBDSelection(ibdobj, kinship.cutoff=1/32) # get an ibd df for relationships > 1/32 kinship coef
head(ibd)
ibd <- merge(ibd, samp, by.x="ID1", by.y="scanID")
head(ibd)
ibd <- merge(ibd, samp, by.x="ID2", by.y="scanID", suffixes=c("1","2"))
head(ibd)
ibd$ii <- pasteSorted(ibd$Individ1, ibd$Individ2)
head(ibd)
# 
relprs$ii <- pasteSorted(relprs$Individ1, relprs$Individ2)
ibd <- merge(ibd, relprs[,c("ii","relation")], all.x=TRUE)
head(ibd)
names(ibd)[names(ibd) == "relation"] <- "exp.rel"
ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"
ibd$exp.rel[is.na(ibd$exp.rel)] <- "U"
head(ibd)
table(ibd$exp.rel, useNA="ifany")

# assign observed relationships
ibd$obs.rel <- ibdAssignRelatednessKing(ibd$IBS0, ibd$kinship)
head(ibd)
table(ibd$obs.rel, useNA="ifany")

# now I can plot it!
## thresholds for assigning relationships using kinship coefficients
## in table 1 of Manichaikul (2010)
cut.dup <- 1/(2^(3/2))
cut.deg1 <- 1/(2^(5/2))
cut.deg2 <- 1/(2^(7/2))
cut.deg3 <- 1/(2^(9/2))
cols <- c(Dup="magenta", PO="cyan", U="black")

plot(ibd$IBS0, ibd$kinship, col=cols[ibd$exp.rel], pch=20, 
     xlab="Fraction of IBS=0", ylab="Kinship coefficient")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=names(cols), col=cols, pch=1)

# 5.3 PCA for population structure -----
filt <- get(data(pcaSnpFilters.hg18)) # regions to filter out due to known population correlation
chrom <- getChromosome(snpAnnot) 
pos <- getPosition(snpAnnot)
snpID <- getSnpID(snpAnnot)
# make an array with T/F values for filtering out snps above
snp.filt <- rep(TRUE, length(snpID))
for (f in 1:nrow(filt)) {
  snp.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos
           & pos < filt$end.base[f]] <- FALSE
}
snp.sel <- snpID[snp.filt]
length(snp.sel) # number of included snps

sample.sel <- scanAnnot$scanID[!scanAnnot$duplicated] # filter out duplicate samples
length(sample.sel)

gdsobj <- snpgdsOpen(gdsfile)
# prune the snps based on ld (just like in PLINK)
snpset <- snpgdsLDpruning(gdsobj, sample.id=sample.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1))
snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned) # total number of snps to use
# now do the PCA - citation looks like it's using the same method as EigenSTRAT
pca <- snpgdsPCA(gdsobj, sample.id=sample.sel, snp.id=snp.pruned)
names(pca)
length(pca$eigenval)
dim(pca$eigenvect)
# Calculate percentage of variance explained by each PC 
pc.frac <- pca$eigenval/sum(pca$eigenval)
lbls <- paste("EV", 1:4, "\n", format(pc.frac[1:4], digits=2), sep="")
samp <- pData(scanAnnot)[match(pca$sample.id, scanAnnot$scanID), ]
cols <- rep(NA, nrow(samp))
cols[samp$race == "CEU"] <- "blue"
cols[samp$race == "YRI"] <- "green"
# Plot the first four PCs in pairs
pairs(pca$eigenvect[, 1:4], col=cols, labels=lbls,
      main = "CEU: blue, YRI: green", pch=20)
# Now make a parallel coordinates plot
par.coord <- pca$eigenvect
rangel <- apply(par.coord, 2, function(x) range(x)[1])
rangeh <- apply(par.coord, 2, function(x) range(x)[2])
std.coord <- par.coord
for (i in 1:14)
  std.coord[,i] <- (par.coord[,i] - rangel[i])/(rangeh[i]-rangel[i])
plot(c(0,15), c(0,1), type = 'n', axes = FALSE, ylab = "", xlab = "",
     main = "Parallel Coordinates Plot
     CEU: blue, YRI: green")
for (j in 1:13)
  for (i in sample(1:nrow(std.coord)) )
    lines(c(j,j+1), std.coord[i,c(j,j+1)], col=cols[i], lwd=0.5)
axis(1, at = 1:14, labels = paste("PC",1:14, sep = "."))
# "he genetic diversity in the YRI group is apparent in the later 
# eigenvectors, while the remaining groups remain in clusters throughout."
# Look for correlations between SNP regions and specific eigenvectors
corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:4)
snpgdsClose(gdsobj)
snp <- snpAnnot[match(corr$snp.id, snpID), ]
chrom <- getChromosome(snp, char=TRUE)
pdf("DataCleaning-corr.pdf")
par(mfrow=c(4,1))
for (i in 1:4) {
  snpCorrelationPlot(abs(corr$snpcorr[i,]), chrom,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

# 6 Case-Control Confounding -------------------
# Check principal components and missing call rate for correlation with case status

# 6.1 Principal Components Differences -----
# make a df of principal components labeled with scanID, case control status, and "race"
princomp <- as.data.frame(pca$eigenvect)
head(princomp)
samples.nodup <- pData(scanAnnot)[!scanAnnot$duplicated,]
princomp$scanID <- as.factor(samples.nodup$scanID)
princomp$case.ctrl.status <- as.factor(samples.nodup$status)
princomp$race <- as.factor(samples.nodup$race)
head(princomp)
# prepare for plotting
pc.percent <- 100 * pca$eigenval[1:32]/sum(pca$eigenval)
pc.percent
lbls <- paste("EV", 1:3, "\n", format(pc.percent[1:3], digits=2), "%", sep="")
table(samples.nodup$status)
cols <- rep(NA, nrow(samples.nodup))
cols[samples.nodup$status == 1] <- "green"
cols[samples.nodup$status == 0] <- "magenta"
# 
pairs(pca$eigenvect[,1:3], col=cols, labels=lbls, pch=20,
      main = "First Three EVs by Case-Control Status")
# No eigenvectors are entireley cases or entirely controls, so this is good
# this bunch of boxplots would be better faceted in ggplot
boxplot(princomp[, 1] ~ princomp$case.ctrl.status,
        ylab = "PC1", main = "PC1 vs. Case-control Status")
boxplot(princomp[, 2] ~ princomp$case.ctrl.status,
        ylab = "PC2", main = "PC2 vs. Case-control Status")
boxplot(princomp[, 3] ~ princomp$case.ctrl.status,
        ylab = "PC3", main = "PC3 vs. Case-control Status")
# look for a relationship with race
aov.p1 <- aov(princomp[,1] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p1)
aov.p2 <- aov(princomp[,2] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p2)
aov.p3 <- aov(princomp[,3] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p3)
# No strong relationship bt first 3 PCs and race -- looks good!

# 6.2 Missing call rate difference -----
# Is there a difference in missingness between cases and controls?
lm.all <- lm(scanAnnot$missing.e1 ~ scanAnnot$status)
summary(aov(lm.all))
boxplot(scanAnnot$missing.e1 ~ scanAnnot$status, ylab = "Mean missing call rate", 
        main = "Mean missing call rate by case status")

# 7 Chromosome anomaly detection -----------------------

# 7.1 BAF filtering -----
blfile <- system.file("extdata", "illumina_bl.gds", package="GWASdata")
blgds <- GdsIntensityReader(blfile)
blData <-  IntensityData(blgds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genogds <- GdsGenotypeReader(genofile)
genoData <-  GenotypeData(genogds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# identify low quality samples from sd(BAF)
baf.sd <- sdByScanChromWindow(blData, genoData, var="BAlleleFreq")
names(baf.sd)
med.baf.sd <- medianSdOverAutosomes(baf.sd)
low.qual.ids <- med.baf.sd$scanID[med.baf.sd$med.sd > 0.05]
# exclude some snps from known problematic regions
chrom <- snpAnnot$chromosome
pos <- snpAnnot$position
hla.df <- get(data(HLA.hg18))
hla <- chrom == "6" & pos >= hla.df$start.base & pos <= hla.df$end.base
xtr.df <- get(data(pseudoautosomal.hg18))
xtr <- chrom == "X" & pos >= xtr.df["X.XTR", "start.base"] &
  pos <= xtr.df["X.XTR", "end.base"]
centromeres <- get(data(centromeres.hg18))
gap <- rep(FALSE, nrow(snpAnnot))
for (i in 1:nrow(centromeres)) {
  ingap <- chrom == centromeres$chrom[i] & pos > centromeres$left.base[i] &
    pos < centromeres$right.base[i]
  gap <- gap | ingap
}
ignore <- snpAnnot$missing.n1 == 1 #ignore includes intensity-only and failed snps
snp.exclude <- ignore | hla | xtr | gap
snp.ok <- snpAnnot$snpID[!snp.exclude]
# anomSegmentBAF uses circular binary segmentation to find change points in BAF 
scan.ids <- scanAnnot$scanID[1:10]
chrom.ids <- 21:23
baf.seg <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
                          chrom.ids=chrom.ids, snp.ids=snp.ok, verbose=FALSE)
head(baf.seg)
# now filter the segments for chromosome anomalies, treating low quality samples differently
baf.anom <- anomFilterBAF(blData, genoData, segments=baf.seg,
                          snp.ids=snp.ok, centromere=centromeres, low.qual.ids=low.qual.ids,
                          verbose=FALSE)
names(baf.anom)
baf.filt <- baf.anom$filtered
head(baf.filt)

# 7.2 Loss of heterozygosity -----
# find LOH regions by looking for homozygous runs with changes in LRR
# (this excludes known anomalies detected by BAF)
# anomDetectLOH looks for LRR change points using circular binary segmentation
loh.anom <- anomDetectLOH(blData, genoData, scan.ids=scan.ids, verbose=FALSE,
                          chrom.ids=chrom.ids, snp.ids=snp.ok, known.anoms=baf.filt)
names(loh.anom)
loh.filt <- loh.anom$filtered
head(loh.filt)

# 7.3 Statistics within anomaolous segments -----
# create a data frame to use here
baf.filt$method <- "BAF"
if (!is.null(loh.filt)) 
{
  loh.filt$method <- "LOH"
  cols <- intersect(names(baf.filt), names(loh.filt))
  anoms <- rbind(baf.filt[ , cols], loh.filt[ , cols])
} else {
  anoms <- baf.filt
}
anoms$anom.id <- 1:nrow(anoms)
# anomSegStats calculates several interesting statistics 
# "Some of these are basic statistics for the characteristics of the anomaly and for 
# measuring deviation of LRR or BAF from expected. Other statistics are used in downstrean 
# quality control analysis, including detecting terminal anomalies and investigating 
# centromere-spanning anomalies."
stats <- anomSegStats(blData, genoData, snp.ids=snp.ok, anom=anoms, centromere=centromeres)
names(stats)
# plot these anomalies
snp.not.ok <- snpAnnot$snpID[snp.exclude]
anomStatsPlot(blData, genoData, anom.stats=stats[1, ], 
              snp.ineligible=snp.not.ok, centromere=centromeres, cex.leg=1)

# 7.4 Identify low quality samples -----
# look for regions with high sd(LRR), including all snps
lrr.sd <- sdByScanChromWindow(blData, var="LogRRatio", incl.hom=TRUE)
med.lrr.sd <- medianSdOverAutosomes(lrr.sd)
# we need to know the number of segments found in anomaly detection
baf.seg.info <- baf.anom$seg.info
loh.seg.info <- loh.anom$base.info[ , c("scanID", "chromosome", "num.segs")]
# add an "eligible" column to the snp annotation df
snpAnnot$eligible <- !snp.exclude
# use different thresholds for BAF and LOH low quality filtering
baf.low.qual <- anomIdentifyLowQuality(snpAnnot, med.baf.sd, baf.seg.info,
                                       sd.thresh=0.1, sng.seg.thresh=0.0008, auto.seg.thresh=0.0001)
loh.low.qual <- anomIdentifyLowQuality(snpAnnot, med.lrr.sd, loh.seg.info,
                                       sd.thresh=0.25, sng.seg.thresh=0.0048, auto.seg.thresh=0.0006)
close(blData)
close(genoData)

# 7.5 Filter out the anomalies -----
# use set MissingGenotypes to filter out these low quality regions
# anomalies to filter
anom.filt <- stats[ , c("scanID", "chromosome", "left.base", "right.base")]
# the required whole.chrom column can be used to filter out sex chromosome anomalies
anom.filt$whole.chrom <- FALSE
# select unique subjects
subj <- scanAnnot$scanID[!scanAnnot$duplicated]
subj.filt.file <- "subj_filt.gds"
# this creates a new GDS file with the chosen genotypes set to missing
setMissingGenotypes(genofile, subj.filt.file, anom.filt, file.type="gds",
                    sample.include=subj, verbose=FALSE)
(gds <- GdsGenotypeReader(subj.filt.file))
close(gds)

# 8 SNP Quality Checks -------------------------
# looking for SNPs subject to genotyping artifacts
# 8.1 Duplicate Sample Discordance -----
snp.excl <- scanAnnot$scanID[scanAnnot$missing.e1 >= 0.05] # exclude snps & samples with high missingness
length(snp.excl)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genoGDS <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(genoGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# calculate duplicate discordance
dupdisc <- duplicateDiscordance(genoData, subjName.col="subjectID", 
                                scan.exclude=scan.excl, snp.exclude=snp.excl)
names(dupdisc)
head(dupdisc$discordance.by.snp)
length(dupdisc$discordance.by.subject)
dupdisc$discordance.by.subject[[2]]
# do some data manipulation because you only care about one value in these symmetric matrices
npair <- length(dupdisc$discordance.by.subject)
disc.subj <- rep(NA, npair)
subjID <- rep(NA, npair)
race <- rep(NA, npair)
for (i in 1:npair) {
  disc.subj[i] <- dupdisc$discordance.by.subject[[i]][1,2]
  subjID[i] <- names(dupdisc$discordance.by.subject)[i]
  race[i] <- scanAnnot$race[scanAnnot$subjectID == subjID[i]][1]
}
dat <- data.frame(subjID=subjID, disc=disc.subj, pop=race,
                  stringsAsFactors=FALSE)
summary(dat$disc)
# Assign colors for the duplicate samples based on population group.
dat$col <- NA
dat$col[dat$pop == "CEU"] <- "blue"
dat$col[dat$pop == "YRI"] <- "green"
dat <- dat[order(dat$disc),]
dat$rank <- 1:npair
# Plot the sample discordance rates color coded by race.
plot(dat$disc, dat$rank, col=dat$col, ylab="rank", pch=20,
     xlab="Discordance rate between duplicate samples",
     main="Duplicate Sample Discordance by Continental Ancestry")
legend("bottomright", unique(dat$pop), pch=rep(20,2), col=unique(dat$col))
# Use the binomial distribution to calculate the prob. of observing > x discordant genotypes
# in a total of n pairs of duplicates
duplicateDiscordanceProbability(npair)

# 8.2 Mendelian error checking -----
## look for parent-offspring genotype combinations that violate Mendelian inheritance
# make a mendelList object
men.list <- with(pData(scanAnnot), mendelList(family, subjectID,
                                              father, mother, sex, scanID))
res <- mendelListAsDataFrame(men.list)
head(res)
dim(res)
# Only want to use SNPs with missing.n1 < 0.05
snp.excl <- snpAnnot$snpID[snpAnnot$missing.n1 >= 0.05]
length(snp.excl)
# identify the mendelian errors
mend <- mendelErr(genoData, men.list, snp.exclude=snp.excl)
names(mend)
head(mend$trios)
names(mend$snp)
# Calculate the error rate per SNP
err <- mend$snp$error.cnt / mend$snp$check.cnt
table(err == 0, useNA="ifany")
#
plot(err, rank(err), xlab="Error Rate (fraction)", pch=20,
     ylab="rank", main="Mendelian Error Rate per SNP, ranked")
# look at errors per family 
fam <- mend$snp$error.cnt
n <- mend$snp$check.cnt
summary(fam)
# SNPs with errors
length(fam[n > 0 & fam > 0])
# SNPs for which more than one family has an error
length(fam[n > 0 & fam > 1])
# Get the SNPs with valid trios for error detection
val <- length(fam[n > 0])
noerr <- length(fam[n > 0 & fam == 0])
# Divide to get fraction with no errors
noerr / val
# add mendelian error values to snp Annotation
snp.sel <- match(names(mend$snp$error.cnt), snpAnnot$snpID)
snpAnnot$mendel.err.count[snp.sel] <- mend$snp$error.cnt
# number of valid families for checking is save as mendel.err.sampsize
snpAnnot$mendel.err.sampsize[snp.sel] <- mend$snp$check.cnt
allequal(snpAnnot$snpID, sort(snpAnnot$snpID))
# The high number of NA values is due to the filtering out of SNPs
#    before the Mendelian error rate calculation
sum(is.na(snpAnnot$mendel.err.count))
sum(is.na(snpAnnot$mendel.err.sampsize))
varMetadata(snpAnnot)["mendel.err.count", "labelDescription"] <-
  paste("number of Mendelian errors detected in trios averaged over",
        "multiple combinations of replicate genotyping instances")
varMetadata(snpAnnot)["mendel.err.sampsize", "labelDescription"] <-
  "number of opportunities to detect Mendelian error in trios"
# For high Mendelian error snps, look at genotype cluster plots to investigate further
# we expect a lack of defined genotype clusters (causes poor call rate)
# Get a vector of SNPs to check
snp <- pData(snpAnnot)
snp$err.rate <- snp$mendel.err.count /
  snp$mendel.err.sampsize
snp <- snp[order(snp$err.rate, decreasing=TRUE),]
snp <- snp[1:9,]
xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
xyGDS <- GdsIntensityReader(xyfile)
xyData <- IntensityData(xyGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
pdf(file="DataCleaning-mendel.pdf")
par(mfrow = c(3,3))
mtxt <- paste("SNP", snp$rsID, "\nMendelian Error Rate",
              format(snp$err.rate, digits=5))
genoClusterPlot(xyData, genoData, snpID=snp$snpID, main.txt=mtxt,
                cex.main=0.9)
dev.off()
close(xyData)
## Mendelian errors per family
# Calculate the fraction of SNPs with an error for each trio
trios <- mend$trios
trios$Mend.err <- trios$Men.err.cnt/trios$Men.cnt
summary(trios$Mend.err)
# Start by pulling out the vectors needed from `trios'
tmp <- trios[, c("fam.id", "Mend.err")]; dim(tmp)
# Change fam.id to match the sample annotation column name
names(tmp) <- c("family", "Mend.err.rate.fam")
# Merge the variables into the sample annotation file
scanAnnot$mend.err.rate.fam <- NA
for (i in 1:nrow(tmp)) {
  ind <- which(is.element(scanAnnot$family, tmp$family[i]))
  scanAnnot$mend.err.rate.fam[ind] <- tmp$Mend.err.rate.fam[i]
}
head(scanAnnot$mend.err.rate.fam)
varMetadata(scanAnnot)["mend.err.rate.fam", "labelDescription"] <-
  "Mendelian error rate per family"
# Look for SNPs that may not be accurately called across all populations for minor allele
# Get the families that have non-NA values for the family
#     Mendelian error rate
fams <- pData(scanAnnot)[!is.na(scanAnnot$mend.err.rate.fam) &
                           !duplicated(scanAnnot$family), c("family",
                                                            "mend.err.rate.fam", "race")]
dim(fams)
table(fams$race, useNA="ifany")
# Assign colors for the different ethnicities in these families
pcol <- rep(NA, nrow(fams))
pcol[fams$race == "CEU"] <- "blue"
pcol[fams$race == "YRI"] <- "green"
#
plot(fams$mend.err.rate.fam*100, rank(fams$mend.err.rate.fam),
     main="Mendelian Error rate per Family, ranked",
     xlab="Mendelian error rate per family (percent)",
     ylab="rank", col=pcol, pch=20)
legend("bottomright", c("CEU", "YRI"), pch=c(20, 20), col=c("blue", "green"))

# 8.3 Hardy-Weinberg equilibrium testing -----
# use Fisher's exact test to detect departure from HWE in each SNP
# for this test, filter out duplicates and non-founders, and run for each population separately
# also, it excludes X in males, Y, and MT snps
head(pData(scanAnnot)[,c("father", "mother")])
nonfounders <- scanAnnot$father != 0 & scanAnnot$mother != 0
table(nonfounders)
scan.excl <- scanAnnot$scanID[scanAnnot$race != "CEU" | nonfounders | scanAnnot$duplicated]
length(scan.excl)
hwe <- gwasExactHW(genoData, scan.exclude=scan.excl)
close(genoData)
#
names(hwe)
dim(hwe)
# Check on sample sizes for autosomes and X chromosome
hwe$N <- hwe$nAA + hwe$nAB + hwe$nBB
summary(hwe$N[is.element(hwe$chromosome,1:22)])
summary(hwe$N[is.element(hwe$chromosome,23)])
hwe$p.value[1:10]
sum(is.na(hwe$p.value[hwe$chromosome == 24])) # XY
sum(is.na(hwe$p.value[hwe$chromosome == 23])) # X
hwe$MAF[1:10]
hwe[1:3, c("nAA", "nAB", "nBB")]
# estimate inbreeding coefficient per SNP
summary(hwe$f)
#
hist(hwe$f, main="Histogram of the Inbreeding Coefficient
     For CEU Samples", xlab="Inbreeding Coefficient")
# histogram centered around 0 indicates hidden population substructure is unlikely
# Check the MAF of those SNPs with f=1
chkf <- hwe[!is.na(hwe$f) & hwe$f==1,]; dim(chkf)
summary(chkf$MAF)
# QQplots visualize at which value SNPs begin to deviate from HWE expectation
hwe.0 <- hwe[hwe$MAF > 0,]; dim(hwe.0) # exclude MAF=0 SNPs
# Only keep the autosomal SNPs for first plot
pval <- hwe.0$p.value[is.element(hwe.0$chromosome, 1:22)]
length(pval)
pval <- pval[!is.na(pval)]
length(pval)
# X chromosome SNPs for plot 2
pval.x <- hwe.0$p.value[is.element(hwe.0$chromosome, 23)]
length(pval.x)
pval.x <- pval.x[!is.na(pval.x)]
length(pval.x)
#
pdf(file = "DataCleaning-hwe.pdf")
par(mfrow=c(2,2))
qqPlot(pval=pval, truncate = FALSE, main="Autosomes, all")
qqPlot(pval=pval, truncate = TRUE, main="Autosomes, truncated")
qqPlot(pval=pval.x, truncate = FALSE, main="X chromosome, all")
qqPlot(pval=pval.x, truncate = TRUE, main="X chromosome, truncated")
dev.off()
# look for a pattern in pvalue vs MAF
plot(hwe.0$MAF, -log10(hwe.0$p.value),
     xlab="Minor Allele Frequency", ylab="-log(p-value)",
     main="Minor Allele Frequency vs\nP-value", pch=20)

# 9 Preliminary Association Tests -------------------
# 9.1 Association test -----
genoGDS <- GdsGenotypeReader(subj.filt.file)
subjAnnot <- scanAnnot[scanAnnot$scanID %in% getScanID(genoGDS),]
genoData <- GenotypeData(genoGDS, scanAnnot=subjAnnot)
# use the filtered GDS file with unique subjects only (from 7.5)
assoc <- assocTestRegression(genoData, outcome="status", covar.list=c(""),
                             ivar.list=c(""), model.type="logistic", robust=TRUE,
                             chromosome.set=c(24:26))
close(genoData)

# 9.2 QQ Plots -----
qqPlot(pval=assoc$model.1.additive.LR.pval.G,
       truncate=TRUE, main="QQ Plot of Likelihood Ratio Test p-values")
# no outliers visible - not surprising since the case-control status was made up

# 9.3 Manhattan plots of p values -----
chrom <- getChromosome(snpAnnot, char=TRUE)
chrom.sel <- chrom %in% c("XY", "Y", "M")
manhattanPlot(assoc$model.1.additive.LR.pval.G,
              chromosome=chrom[chrom.sel])

# 9.4 SNP cluster plots -----
# for the strongest hits, look at the genotype clustering plots to make sure things look right
# Identify SNPs with lowest p-values
snp <- pData(snpAnnot)[chrom.sel, c("snpID", "rsID")]
allequal(snp$snpID, assoc$snpID)
snp$pval <- assoc$model.1.additive.LR.pval.G
snp <- snp[order(snp$pval),]
snp <- snp[1:9,]
xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
xyGDS <- GdsIntensityReader(xyfile)
xyData <- IntensityData(xyGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
genoGDS <- GdsGenotypeReader(genofile)
genoData <- GenotypeData(genoGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
pdf(file="DataCleaning-cluster.pdf")
par(mfrow = c(3,3))
mtxt <- paste("SNP", snp$rsID, "\np =", format(snp$pval, digits=4))
genoClusterPlot(xyData, genoData, snpID=snp$snpID, main.txt=mtxt)
dev.off()
close(xyData)
close(genoData)


unlink(subj.filt.file)


