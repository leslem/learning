### R code from vignette source 'vignettes/GWASTools/inst/doc/Affymetrix.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Affymetrix.Rnw:68-88
###################################################
library(GWASTools)
library(GWASdata)
# Load the SNP annotation (simple data frame)
data(affy_snp_annot)
# Create a SnpAnnotationDataFrame
snpAnnot <- SnpAnnotationDataFrame(affy_snp_annot)
# names of columns
varLabels(snpAnnot)
# data
head(pData(snpAnnot))
# Add metadata to describe the columns
meta <- varMetadata(snpAnnot)
meta[c("snpID", "chromosome", "position", "rsID", "probeID"), 
  "labelDescription"] <- c("unique integer ID for SNPs",                      
  paste("integer code for chromosome: 1:22=autosomes,",
   "23=X, 24=pseudoautosomal, 25=Y, 26=Mitochondrial, 27=Unknown"), 
  "base pair position on chromosome (build 36)",
  "RS identifier",
  "unique ID from Affymetrix")
varMetadata(snpAnnot) <- meta


###################################################
### code chunk number 2: Affymetrix.Rnw:110-137
###################################################
# Load the scan annotation (simple data frame)
data(affy_scan_annot)
# Create a ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(affy_scan_annot)
# names of columns
varLabels(scanAnnot)
# data
head(pData(scanAnnot))
# Add metadata to describe the columns
meta <- varMetadata(scanAnnot)
meta[c("scanID", "subjectID", "family", "father", "mother",
  "CoriellID", "race", "sex", "status", "genoRunID", "plate",
  "alleleFile", "chpFile"), "labelDescription"] <- 
   c("unique integer ID for scans", 
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
  "data file with intensities",
  "data file with genotypes and quality scores")
varMetadata(scanAnnot) <- meta


###################################################
### code chunk number 3: Affymetrix.Rnw:160-240
###################################################
geno.file <- "tmp.geno.gds"

# first 3 samples only
scan.annotation <- affy_scan_annot[1:3, c("scanID", "genoRunID", "chpFile")]
names(scan.annotation) <- c("scanID", "scanName", "file")

# indicate which column of SNP annotation is referenced in data files
snp.annotation <- affy_snp_annot[,c("snpID", "probeID", "chromosome", "position")]
names(snp.annotation)[2] <- "snpName"

# col.nums is an integer vector indicating which columns of the raw text file 
# contain variables for input
col.nums <- as.integer(c(2,3))
names(col.nums) <- c("snp", "geno")

# Define a path to the raw data CHP text files
path <- system.file("extdata", "affy_raw_data", package="GWASdata")

diag.geno.file <- "diag.geno.RData"
diag.geno <- createDataFile(path = path,
  filename = geno.file, 
  file.type = "gds",
  variables = "genotype",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation, 
  sep.type = "\t",
  skip.num = 1,
  col.total = 6,
  col.nums = col.nums,
  scan.name.in.file = -1,
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
#   checks were successful and the data were written to the file
table(diag.geno$chk)

# open the file we just created
(genofile <- GdsGenotypeReader(geno.file))
# Take out genotype data for the first 3 samples and 
#   the first 5 SNPs
(genos <- getGenotype(genofile, snp=c(1,5), scan=c(1,3)))
close(genofile)

# check that the data file matches the raw text files
check.geno.file <- "check.geno.RData"
check.geno <- checkGenotypeFile(path = path,
  filename = geno.file, 
  file.type = "gds",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation, 
  sep.type = "\t",
  skip.num = 1,
  col.total = 6,
  col.nums = col.nums,
  scan.name.in.file = -1,
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


###################################################
### code chunk number 4: Affymetrix.Rnw:258-311
###################################################
qual.file <- "tmp.qual.gds"

scan.annotation <- affy_scan_annot[1:3, c("scanID", "genoRunID", "chpFile")]
names(scan.annotation) <- c("scanID", "scanName", "file")

snp.annotation <- affy_snp_annot[,c("snpID", "probeID", "chromosome", "position")]
names(snp.annotation)[2] <- "snpName"

col.nums <- as.integer(c(2,4))
names(col.nums) <- c("snp", "quality")

path <- system.file("extdata", "affy_raw_data", package="GWASdata")

diag.qual.file <- "diag.qual.RData"
diag.qual <- createDataFile(path = path,
  filename = qual.file, 
  file.type = "gds",
  variables = "quality",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation, 
  sep.type = "\t",
  skip.num = 1,
  col.total = 6,
  col.nums = col.nums,
  scan.name.in.file = -1,
  diagnostics.filename = diag.qual.file,
  verbose = FALSE)

# Open the GDS file we just created
(qualfile <- GdsIntensityReader(qual.file))
# Take out the quality scores for the first 
#     5 SNPs for the first 3 samples
(qual <- getQuality(qualfile, snp=c(1,5), scan=c(1,3)))
close(qualfile)

# check
check.qual.file <- "check.qual.RData"
check.qual <- checkIntensityFile(path = path,
  filename = qual.file, 
  file.type = "gds",
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation, 
  sep.type = "\t",
  skip.num = 1,
  col.total = 6,
  col.nums = col.nums,
  scan.name.in.file = -1,
  check.scan.index = 1:3,
  n.scans.loaded = 3,
  affy.inten = FALSE,
  diagnostics.filename = check.qual.file,
  verbose = FALSE)



###################################################
### code chunk number 5: Affymetrix.Rnw:317-359
###################################################
xy.file <- "tmp.xy.gds"

# we need the allele files instead of CHP files
scan.annotation <- affy_scan_annot[1:3, c("scanID", "genoRunID", "alleleFile")]
names(scan.annotation) <- c("scanID", "scanName", "file")

snp.annotation <- affy_snp_annot[,c("snpID", "probeID", "chromosome", "position")]
names(snp.annotation)[2] <- "snpName"

diag.xy.file <- "diag.xy.RData"
diag.xy <- createAffyIntensityFile(path = path,
  filename = xy.file, 
  file.type = "gds",
  snp.annotation = snp.annotation, 
  scan.annotation = scan.annotation, 
  diagnostics.filename = diag.xy.file,
  verbose = FALSE)

# Open the file we just created
(intenfile <- GdsIntensityReader(xy.file))
# Take out the normalized X intensity values for the first 
#     5 SNPs for the first 3 samples
(xinten <- getX(intenfile, snp=c(1,5), scan=c(1,3)))
close(intenfile)

# check the intensity files - note affy.inten=TRUE
check.xy.file <- "check.xy.RData"
check.xy <- checkIntensityFile(path = path,
  filename = xy.file,
  file.type = "gds", 
  snp.annotation = snp.annotation,
  scan.annotation = scan.annotation, 
  sep.type = "\t",
  skip.num = 1,
  col.total = 2,
  col.nums = setNames(as.integer(c(1,2,2)), c("snp", "X", "Y")),
  scan.name.in.file = -1,
  check.scan.index = 1:3,
  n.scans.loaded = 3,
  affy.inten = TRUE,
  diagnostics.filename = check.xy.file,
  verbose = FALSE)


###################################################
### code chunk number 6: Affymetrix.Rnw:445-463
###################################################
bl.file <- "tmp.bl.gds"
xyData <- IntensityData(GdsIntensityReader(xy.file),
                        snpAnnot=snpAnnot)
genoData <- GenotypeData(GdsGenotypeReader(geno.file),
                         snpAnnot=snpAnnot)
BAFfromGenotypes(xyData, genoData,
                 filename = bl.file, 
                 file.type = "gds",
                 min.n.genotypes = 0,
                 call.method ="by.study",
                 verbose = FALSE)
close(xyData)
close(genoData)
# Open the file we just created
(blfile <- GdsIntensityReader(bl.file))
# Look at the BAlleleFrequency values for the first 5 SNPs
(baf <- getBAlleleFreq(blfile, snp=c(1,5), scan=c(1,3)))
close(blfile)


###################################################
### code chunk number 7: Affymetrix.Rnw:468-471
###################################################
file.remove(geno.file, qual.file, xy.file, bl.file)
file.remove(diag.geno.file, diag.qual.file, diag.xy.file)
file.remove(check.geno.file, check.qual.file, check.xy.file)


