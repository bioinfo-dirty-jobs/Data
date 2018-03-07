#!/usr/bin/env Rscript

# Code modified from 
# https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

## TO DO:
## - manage multiple conditions (e.g. batch + treatment): the current implementation works with one condition only

suppressMessages(library("sleuth"))

kal_folder <- snakemake@output[["dir"]]
labels     <- unlist(strsplit(snakemake@params[["classes"]],","))
c1         <- unlist(strsplit(snakemake@params[["class1"]],","))
c2         <- unlist(strsplit(snakemake@params[["class2"]],","))

# direct output to a file
o2f <- file(snakemake@log[[1]], open = "wt")
sink(o2f, type="message", append=FALSE, split=FALSE)

# Get transcript to gene information
# We assume that:
# - kallisto used the ensembl transcript identifiers

if (Biobase::testBioCConnection()){
	
	# Connect to biomaRt
	ensembl <- biomaRt::useEnsembl(biomart=snakemake@params[["database"]],
                                   dataset=snakemake@params[["dataset"]],
                                   version=snakemake@params[["version"]])

	# Get transcript to gene info
	t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
										 "ensembl_gene_id",
										 "external_gene_name",
										 "ensembl_transcript_id_version"),
										 mart = ensembl)
	# Rename columns									 
	t2g <- dplyr::rename(t2g,
						 target_id = ensembl_transcript_id_version,
						 ens_gene = ensembl_gene_id,
						 ext_gene = external_gene_name)

} else {
	stop("transcript2gene file is missing, and there is no internet connection", call.=FALSE)
}

# Get sample names from classes
sample_id <- c(c1,c2)

# Path to kallisto results
kal_dirs <- file.path(kal_folder,sample_id)

# Conditions
sample_cond <- c(rep(labels[1],length(c1)),
				 rep(labels[2],length(c2))
				)

# Build sample to condition dataframe
s2c <- data.frame(sample=sample_id,
				  condition=sample_cond,
				  path=kal_dirs,
				  stringsAsFactors=FALSE)

# Load the kallisto processed data into the object
so <- sleuth_prep(s2c, target_mapping = t2g)

# Estimate parameters for the sleuth response error measurement (full) model
so <- sleuth_fit(so, ~condition, 'full')

# Perform test
which_beta <- paste("condition",labels[2],sep="")
so <- sleuth_wt(so,  which_beta, 'full')

# Extract results
sleuth_table <- sleuth_results(so, test = which_beta, which_model = 'full', show_all = FALSE)

# Create gene table
gene_table <- sleuth_gene_table(so, test = which_beta, test_type = "wt", which_model = 'full')

# return output to the terminal
sink(type = "message")
sink()

## Write results
# Gene table
write.table(gene_table, paste(snakemake@output[["gene_table"]]), sep="\t")

# Diff expression results
write.table(sleuth_table, paste(snakemake@output[["sleuth_table"]]), sep="\t")

# Sleuth results as R object
save(so, file=paste(snakemake@output[["sleuth_object"]]))

