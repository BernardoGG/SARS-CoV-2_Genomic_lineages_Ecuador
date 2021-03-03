#############################################################################################################################
######### SARS-CoV-2 GISAID downsampler #####################################################################################
#############################################################################################################################

######### Darlan da Silva Candido & Bernardo Gutierrez ######################################################################

library(seqinr)
library(dplyr)
library(tidyverse)
library(data.table)
library(plyr)
library(lubridate)

######### Read metadata and sequence data ###################################################################################

# Read files
meta1 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_15_(1).tsv", sep = "\t", strip.white = TRUE)
meta2 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_15_(2).tsv", sep = "\t", strip.white = TRUE)
meta3 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_15_(3).tsv", sep = "\t", strip.white = TRUE)
meta4 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(4).tsv", sep = "\t", strip.white = TRUE)
meta5 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(5).tsv", sep = "\t", strip.white = TRUE)
meta6 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(6).tsv", sep = "\t", strip.white = TRUE)
meta7 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(7).tsv", sep = "\t", strip.white = TRUE)
meta8 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(8).tsv", sep = "\t", strip.white = TRUE)
meta9 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(9).tsv", sep = "\t", strip.white = TRUE)
meta10 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(10).tsv", sep = "\t", strip.white = TRUE)
meta11 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(11).tsv", sep = "\t", strip.white = TRUE)
meta12 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(12).tsv", sep = "\t", strip.white = TRUE)
meta13 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(13).tsv", sep = "\t", strip.white = TRUE)
meta14 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(14).tsv", sep = "\t", strip.white = TRUE)
meta15 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(15).tsv", sep = "\t", strip.white = TRUE)
meta16 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(16).tsv", sep = "\t", strip.white = TRUE)
meta17 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(17).tsv", sep = "\t", strip.white = TRUE)
meta18 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(18).tsv", sep = "\t", strip.white = TRUE)
meta19 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(19).tsv", sep = "\t", strip.white = TRUE)
meta20 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(20).tsv", sep = "\t", strip.white = TRUE)
meta21 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(21).tsv", sep = "\t", strip.white = TRUE)
meta22 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(22).tsv", sep = "\t", strip.white = TRUE)
meta23 <- read.csv("Data/genetic_data/GISAID_complete_highcoverage_20210107/gisaid_hcov-19_2021_01_07_16_(23).tsv", sep = "\t", strip.white = TRUE)

meta <- do.call("rbind", list(meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10,
                              meta11, meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19,
                              meta20, meta21, meta22, meta23)) # <-- creates a single metadata file

meta_ec <- read.csv("Data/genetic_data/GISAID_hCoV-19_20210107_EC_all.tsv", sep = "\t", strip.white = TRUE)

rm(meta1, meta2, meta3, meta4, meta5, meta6, meta7, meta8, meta9, meta10,
   meta11, meta12, meta13, meta14, meta15, meta16, meta17, meta18, meta19,
   meta20, meta21, meta22, meta23)

seqs <- read.fasta(file = "Data/genetic_data/GISAID_hCoV-19_20210107_full.fasta", seqtype = "DNA", as.string = TRUE,
                  set.attributes = FALSE) # <-- import sequences file (w/o sequences from Ecuador)

seqs_ec <- read.fasta(file = "Data/genetic_data/GISAID_hCoV-19_20210107_EC.fasta", seqtype = "DNA", as.string = TRUE,
                      set.attributes = FALSE) # <-- import Ecuador sequences file

# All sequences (not from Ecuador) must be longer than 22250 (75% of the hCoV-19/Wuhan/IVDC-HB-01/2019|EPI_ISL_402119|2019-12-30,
# 29891)

######### Clean and filter metadata #########################################################################################

# Separate locations into Continent, Country, Region and City
meta <- meta %>%
  separate(Location, c("Continent", "Country", "Region", "City"), sep = "/") # <-- expect warning for added 'NA's

meta_ec <- meta_ec %>%
  separate(Location, c("Continent", "Country", "Region", "City"), sep = "/") # <-- expect warning for added 'NA's

# Remove white spaces
meta_op <- as.data.frame(apply(meta, 2, function(x) gsub('\\s+', '',x)))

# Clean date format
meta_op$Collection.date <- ymd(meta_op$Collection.date) # <-- no "NNNN failed to parse" message expected as all sequences 
                                                        #     should contain a collection date

# Check country names
unique(meta_op$Country) # <-- 135 countries expected

######### Create data frame with column that equals the fasta file sequence ID field ########################################

# Merge virus name, country and collection date to match sequence names in fasta file
meta_op2 <- unite(meta_op, Sequence_id, Virus.name, Accession.ID, Collection.date, sep = "|")
meta_op2 <- cbind(meta_op2, meta_op$Collection.date)
colnames(meta_op2)[16] <- "Collection.date"
meta_op2 <- meta_op2[meta_op2$Country!="Ecuador",] # <-- filter out sequences from Ecuador that were not originally excluded from the metadata

meta_ec2 <- unite(meta_ec, Sequence_id, Virus.name, Accession.ID, Collection.date, sep = "|")
meta_ec2 <- cbind(meta_ec2, meta_ec$Collection.date)
colnames(meta_ec2)[16] <- "Collection.date"

######### Subset sequences based on the filtering applied to the metadata ###################################################

## Downsample metadata to include one sequence per country per day (systematic downsample)
meta_down <-  plyr::ddply(meta_op2,.(Country, Collection.date),function(x) x[sample(nrow(x),1),]) # <-- systematic downsampling

meta_down_rate <- meta_down[sample(nrow(meta_down), round(nrow(meta_down)*0.1)),] # <-- 10% of systematic downsampling set for rate estimation

## Downsample metadata to include 5559 random sequences
meta_down_rand1 <-  meta_op2[sample(nrow(meta_op2), nrow(meta_down)), ] # <-- random subsampling no.1
meta_down_rand2 <-  meta_op2[sample(nrow(meta_op2), nrow(meta_down)), ] # <-- random subsampling no.2
meta_down_rand3 <-  meta_op2[sample(nrow(meta_op2), nrow(meta_down)), ] # <-- random subsampling no.3

## Downsample metadata to include the global lineages observed in Ecuador (one data set per lineage)
setDT(meta_op2)
meta_B.1.1 <- meta_op2[meta_op2$Lineage=="B.1.1",]
meta_B.1 <- meta_op2[meta_op2$Lineage=="B.1",]
meta_B.1.5.6 <- meta_op2[meta_op2$Lineage=="B.1.5.6",]
meta_B.1.1.1 <- meta_op2[meta_op2$Lineage=="B.1.1.1",]
meta_B.1.102 <- meta_op2[meta_op2$Lineage=="B.1.102",]
meta_B.1.67 <- meta_op2[meta_op2$Lineage=="B.1.67",]
meta_B.1.29 <- meta_op2[meta_op2$Lineage=="B.1.29",]

#meta_down_B.1 <-  plyr::ddply(meta_B.1,.(Country, Collection.date),function(x) x[sample(nrow(x),1),]) # <-- B.1 systematic downsampling
#meta_down_B.1.1 <-  plyr::ddply(meta_B.1.1,.(Country, Collection.date),function(x) x[sample(nrow(x),1),]) # <-- B.1.1 systematic downsampling


######### Generate sequence sets with full sampling dates ###################################################################

## Systematically downsampled data set
subset_seqs_syst <- seqs[names(seqs) %in% meta_down$Sequence_id]
#write.fasta(sequences = c(subset_seqs_syst, seqs_ec), names = names(c(subset_seqs_syst, seqs_ec)), nbchar = 300,
#            file.out = "Data/genetic_data/GISAID_hCoV-19_20210107_syst_subset.fasta")

subset_seqs_syst_rate <- seqs[names(seqs) %in% meta_down_rate$Sequence_id]
#write.fasta(sequences = subset_seqs_syst_rate, names = names(subset_seqs_syst_rate), nbchar = 300,
#            file.out = "Data/genetic_data/GISAID_hCoV-19_20210107_syst_subset_rates.fasta") # <-- random 10% of systematically downsampled
#                                                                                            # data set for rate estimation

## Randomly downsampled data sets
subset_seqs_rand1 <- seqs[names(seqs) %in% meta_down_rand1$Sequence_id]
#write.fasta(sequences = c(subset_seqs_rand1, seqs_ec), names = names(c(subset_seqs_rand1, seqs_ec)), nbchar = 300,
#            file.out = "Data/genetic_data/GISAID_hCoV-19_20210107_random_subset_1.fasta")

subset_seqs_rand2 <- seqs[names(seqs) %in% meta_down_rand2$Sequence_id]
#write.fasta(sequences = c(subset_seqs_rand2, seqs_ec), names = names(c(subset_seqs_rand2, seqs_ec)), nbchar = 300,
#            file.out = "Data/genetic_data/GISAID_hCoV-19_20210107_random_subset_2.fasta")

subset_seqs_rand3 <- seqs[names(seqs) %in% meta_down_rand3$Sequence_id]
#write.fasta(sequences = c(subset_seqs_rand3, seqs_ec), names = names(c(subset_seqs_rand3, seqs_ec)), nbchar = 300,
#            file.out = "Data/genetic_data/GISAID_hCoV-19_20210107_random_subset_3.fasta")

## Global lineage data sets
B.1.1_seqs <- seqs[names(seqs) %in% meta_B.1.1$Sequence_id]
#write.fasta(sequences = c(B.1.1_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.1"]]),
#            names = names(c(B.1.1_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.1"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1.1_EC.fasta")

B.1_seqs <- seqs[names(seqs) %in% meta_B.1$Sequence_id]
#write.fasta(sequences = c(B.1_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1"]]),
#            names = names(c(B.1_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1_EC.fasta")

B.1.5.6_seqs <- seqs[names(seqs) %in% meta_B.1.5.6$Sequence_id]
#write.fasta(sequences = c(B.1.5.6_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.5.6"]]),
#            names = names(c(B.1.5.6_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.5.6"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1.5.6_EC.fasta")

B.1.1.1_seqs <- seqs[names(seqs) %in% meta_B.1.1.1$Sequence_id]
#write.fasta(sequences = c(B.1.1.1_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.1.1"]]),
#            names = names(c(B.1.1.1_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.1.1"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1.1.1_EC.fasta")

B.1.102_seqs <- seqs[names(seqs) %in% meta_B.1.102$Sequence_id]
#write.fasta(sequences = c(B.1.102_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.102"]]),
#            names = names(c(B.1.102_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.102"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1.102_EC.fasta")

B.1.67_seqs <- seqs[names(seqs) %in% meta_B.1.67$Sequence_id]
#write.fasta(sequences = c(B.1.67_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.67"]]),
#            names = names(c(B.1.67_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.67"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1.67_EC.fasta")

B.1.29_seqs <- seqs[names(seqs) %in% meta_B.1.29$Sequence_id]
#write.fasta(sequences = c(B.1.29_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.29"]]),
#            names = names(c(B.1.29_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.29"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_B.1.29_EC.fasta")

B.1_down_seqs <- seqs[names(seqs) %in% meta_down_B.1$Sequence_id]
#write.fasta(sequences = c(B.1_down_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1"]]),
#            names = names(c(B.1_down_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_syst_subset_B.1_EC.fasta")

B.1.1_down_seqs <- seqs[names(seqs) %in% meta_down_B.1.1$Sequence_id]
#write.fasta(sequences = c(B.1.1_down_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.1"]]),
#            names = names(c(B.1.1_down_seqs, seqs_ec[names(seqs_ec) %in% meta_ec2$Sequence_id[meta_ec2$Lineage=="B.1.1"]])),
#            nbchar = 300, file.out = "GISAID_hCoV-19_2020-10_30-18_syst_subset_B.1.1_EC.fasta")

######### Cleanup and data export #########################################################################################

write_csv(meta_ec2, "metadata_genomes_EC.csv")
rm(c(sources_continent, sources_country, seqs, seqs_ec, meta, meta_op, meta_ec))
