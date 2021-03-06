## Dada2 and phyloseq demo for Ganda lab
# Emily Bean, 11/2020 - updated for Capsicum by Stephanie Clouser 2/2022

# ---- install and read data ----

### dada2 and phyloseq are both Bioconductor packages

# BiocManager was installed manually in packages tab

# load package
require(BiocManager)

# install dada2
BiocManager::install("dada2")

# install from sources
no

# load the package
require(dada2)

# load our environment
load("Capsicum-Environment.RData")


# path to directory that contains fastq files
# Github doesn't handle file storage well so these are on OneDrive
# PATH = "~/The Pennsylvania State University/Ganda, Erika - Shared-Trello-Projects/Mouse_Opioid/Mouse_Opioid_Reads/"
PATH <- "~/Desktop/16s-workshop/"

# test this variable
PATH

# let's see if anything is in this folder
ourlist <- list.files(PATH)

# how many fastq files are in our folder
length(ourlist)

#  set pathway to the database you prefer - we will use Greengenes
## NOTE:Github doesn't handle fastq.gz files well so download these locally
# dada2 maintains a list of databases: https://benjjneb.github.io/dada2/training.html
DB = "~/Downloads/silva_nr99_v138.1_train_set.fa"


# paired end characterization; most Illumina files are sample names + "_R1_001.fastq" for forward reads
# however, sequences downloaded from NCBI have patterns: "_1.fastq" for forward and "_2.fastq" for reverse
PATTERNF = "_1.fq.gz"
PATTERNR = "_2.fq.gz"

## ---- getFiles----

# get forward and reverse files
## for this demo, subset to the first 2 forward and 2 reverse files
fwdFiles <- list.files(PATH, pattern = PATTERNF, full.names = TRUE)[1:10]
revFiles <- list.files(PATH, pattern = PATTERNR, full.names = TRUE)[1:10]

# check to make sure that the lengths of both files are the same
if(length(fwdFiles) != length(revFiles)) {
  
  stop("There is an unequal number of forward and reverse files")
}

# get sample names
fwdNames <- sapply(strsplit(basename(fwdFiles), PATTERNF), `[`, 1)
revNames <- sapply(strsplit(basename(revFiles), PATTERNR), `[`, 1)

## NOTE: DEFAULT CODE ASSUMES FWD AND REV FILES ARE ORDERED
# error catch if unordered
if(any(!fwdNames %in% revNames)) {
  
  stop("forward and reverse files are out of order")
  
}

### ---- filterAndTrim ----

# create subdirectory for filtered files
filtForward <- file.path("~/Desktop/16s-workshop/filtered", paste0(fwdNames, "_F_filt.fastq.gz"))
filtReverse <- file.path("~/Desktop/16s-workshop/filtered", paste0(revNames, "_R_filt.fastq.gz"))

## Dada2 can plot Phred qualities but FastQC is much better
plotQualityProfile(fwdFiles[1])

# white blocks in the plot are because it is binned quality scores...additional steps are needed.

# filter and trim
cleaned <- filterAndTrim(
  # set forward and reverse paths
  fwd = fwdFiles, rev = revFiles,
  # set paths for the filtered files that will be created
  filt = filtForward, filt.rev = filtReverse,
  # add any necessary filtering parameters, dada2 does not allow ambiguous bases
  maxN = 0,
  
  # set maxEE to 2 and minimum length to 100bp
  maxEE = 2,
  minLen = 100,
  
  # MAC ONLY: multithread ability
  multithread = TRUE, 
  verbose = TRUE
)

# path to filtered and cleaned reads
CLEANEDPATH = "~/Desktop/16s-workshop/filtered/"

# visualize the quality of filtered data on a forward & reverse read
# in reality: do this step in FastQC as well
plotQualityProfile(CLEANEDPATH[1])

# good spot to do another fastqc check in terminal on the trimmed data


## -----------------dada2 algorithm------------------

### NOTE: some patterns are re-done because this script can also be run
# after filtering with Trimmomatic or another program like QIIME

# pattern that specifies which reads are forward or reverse
# if single-read pairs, only specifyforward 
FILTEREDF = "_F_filt.fastq.gz"
FILTEREDR = "_R_filt.fastq.gz"

# get forward and reverse reads
forward <- sort(list.files(CLEANEDPATH, pattern = FILTEREDF, full.names = TRUE))
reverse <- sort(list.files(CLEANEDPATH, pattern = FILTEREDR, full.names = TRUE))

# check to make sure that the lengths of both files are the same and that they match names
fwdNames <- sapply(strsplit(basename(forward), FILTEREDF), `[`, 1)
revNames <- sapply(strsplit(basename(reverse), FILTEREDR), `[`, 1)

# error catch

if(length(fwdNames) != length(revNames)) {
  stop("The number of forward and reverse files do not match.")
} else {
  
  if(any(!fwdNames%in% revNames)) {
    
    stop("Forward and reverse reads are out of order.")
  }
}

# perform error learning
errF <- learnErrors(forward, 
                    multithread = TRUE,
                    verbose = TRUE)
errR <- learnErrors(reverse, 
                    multithread = TRUE,
                    verbose = TRUE)

# visualize error plots - with binned quality score, it will look bad: we want the black lines
# to match up with the red lines
plotErrors(errF, nominalQ = TRUE)

# visualize reverse plots
plotErrors(errR, nominalQ = TRUE)

# perform denoising on forward and reverse reads (sample interference)
dadaForward <- dada(derep = forward, 
                    err = errF, 
                    multithread = TRUE)
dadaReverse <- dada(derep = reverse, 
                    err = errR, 
                    multithread = TRUE)

# merge paired reads
mergers <- mergePairs(dadaF = dadaForward,
                      derepF = forward,
                      dadaR = dadaReverse,
                      derepR = reverse,
                      verbose = TRUE)

# construct sequence table of ASVs
seqtab <- makeSequenceTable(samples = mergers)

# view the dimensions of the table
dim(seqtab)

# inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(unqs = seqtab, 
                                    method = "consensus",
                                    multithread = TRUE,
                                    verbose = TRUE)

# view dimensions of seqence table with no chimeras
dim(seqtab.nochim)

# assign taxonomy using the Silva database
tax <- assignTaxonomy(seqs = seqtab.nochim, 
                      refFasta = DB, 
                      multithread = TRUE,
                      verbose = TRUE)

# OPTIONAL: add species to table using Silva database
taxa_species <- addSpecies(tax, "~/Downloads/silva_species_assignment_v138.1.fa",
                           verbose = TRUE)

# export the dataframe to csv file 
write.csv(tax, "PATH\\taxa10.csv", row.names = FALSE)

# export the datafram with species to csv file
write.csv(taxa_species, "PATH\\taxa_species10.csv", row.names = FALSE)

# end of dada2 - export objects
save(seqtab.nochim, file = "asv-table.RData")
save(tax, file = "taxonomy-table.RData")

## This is the end of the dada2 algorithm

## --------phyloseq analysis-------------------

# install phyloseq package
BiocManager::install("phyloseq")
no

# activate package
require(phyloseq)

### ---- phyloseq ----

# create demo sample data - in reality, you would have metadata
# associated with your dataset
sampleDF <- data.frame(sample = sapply(strsplit(rownames(seqtab.nochim), "_F"), `[`, 1),
                       treatment = c("Treatment", "Control"),
                       sex = c("Male", "Female"),
                       row.names = rownames(seqtab.nochim))

# make phyloseq object
## NOTE: Phyloseq has wrappers for import directly from QIIME among others
# check phyloseq Github for more info

## Phyloseq objects need 3 things: ASV table, taxonomy table, and sample data
# OPTIONAL: add a phylogenetic tree
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(sampleDF), 
               tax_table(tax))

# default settings are messy; replace strings of DNA with an ASV sequence 
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# examine the phyloseq object
ps

# to examine parts of phyloseq object, use "@"
ps@sam_data

# to get columns within the sample data, use "@" followed by "$"
ps@sam_data$treatment

# what does the ASV table look like?
head(ps@otu_table)

# the taxonomy table serves as a lookup table for ASVs
head(ps@tax_table)

# how many Phyla were classified?
get_taxa_unique(physeq = ps,
                taxonomic.rank = "Phylum")

## NOTE: phyloseq plots are wrappers for ggplot and take gg commands
#Abundance barplot by Sample
plot_bar(physeq = ps, 
         x = "sample", 
         y = "Abundance", 
         fill = "Phylum", 
         title = "Abundance Barplot by Sample")

# Alpha diversity plot
plot_richness(ps, x="treatment", measures=c("Shannon", "Simpson"), color="treatment")
