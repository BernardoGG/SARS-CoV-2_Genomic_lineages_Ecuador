#############################################################################################################################
######### SARS-CoV-2 phylogenetics metadata preparation #####################################################################
#############################################################################################################################

######### Bernardo Gutierrez ################################################################################################

source("GISAID_downsampler.R")

######### Create data frames ################################################################################################

## Geographical locations
# Create data frame for sequence location (Ecuador vs non-Ecuador)

intros <- rbind(data.frame(seq = meta_op2$Sequence_id,
                           loc1 = rep("nonEC", nrow(meta_op2))),
                data.frame(seq = meta_ec2$Sequence_id,
                           loc1 = rep("EC", nrow(meta_ec2))))
colnames(intros) <- c("seq", "source")

# Create data frame for sequence location (continent)

sources_continent <- rbind(data.frame(seq = meta_op2$Sequence_id,
                                      loc1 = meta_op2$Continent),
                           data.frame(seq = meta_ec2$Sequence_id,
                                      loc1 = rep("Ecuador", nrow(meta_ec2))))
colnames(sources_continent) <- c("seq", "continent")

# Create data frame for sequence location (country)

sources_country <- rbind(data.frame(seq = meta_op2$Sequence_id,
                                    loc1 = meta_op2$Country),
                         data.frame(seq = meta_ec2$Sequence_id,
                                    loc1 = meta_ec2$Country))
colnames(sources_country) <- c("seq", "country")

## Temporal information
# Create data frame for sequence collection dates

dates <- rbind(data.frame(seq = meta_op2$Sequence_id,
                          date = meta_op2$Collection.date),
               data.frame(seq = genomic_EC$Sequence_id,
                          date = genomic_EC$Collection.date))

## Ecuador full metadata (for descriptive plots)
genomic_EC_manual <- read.csv("Data/genetic_data/metadata/metadata_genomes_manualinput_EC.csv", sep = ",", strip.white = TRUE) %>% arrange(seq_id)

#genomic_USFQ <- read.csv("DB_Genomes_SARS-CoV2.csv", sep = ",") %>%
#  separate(seq_id, c("Country_name", "Seq", "Year"), sep = "/") %>%
#  select(-Country_name, -Year, -lineage, -sex, -age) %>% arrange(Seq)

genomic_EC <- meta_ec2 %>% select(-Region, -Continent, -Country, -City, -Host, -Passage, -Specimen, -Clade,
                                  -Additional.host.information, -Additional.location.information, -Collection.date) %>% arrange(Sequence_id)

genomic_EC$Gender <- as.factor(genomic_EC$Gender)
genomic_EC$Patient.age <- as.numeric(genomic_EC$Patient.age) # <-- expect NAs introduced by coercion

genomic_EC$Lineage <- genomic_EC_manual$lineage %>% as.factor()
genomic_EC$Collection.date <- genomic_EC_manual$Collection.date %>% as.Date()
genomic_EC$Cluster <- genomic_EC_manual$cluster %>% as.factor()
genomic_EC$Region <- genomic_EC_manual$region %>% as.factor()
genomic_EC$Province <- genomic_EC_manual$admin_1 %>% as.factor()
genomic_EC$Lat1 <- genomic_EC_manual$lat_1
genomic_EC$Long1 <- genomic_EC_manual$long_1
genomic_EC$City <- genomic_EC_manual$admin_2 %>% as.factor()
genomic_EC$Lat2 <- genomic_EC_manual$lat_2
genomic_EC$Long2 <- genomic_EC_manual$long_2
#genomic_EC$Admin3 <- genomic_USFQ$admin_3 %>% as.factor()
#genomic_EC$Admin4 <- genomic_USFQ$admin_4 %>% as.factor()
#genomic_EC$Latf <- genomic_USFQ$lat_f
#genomic_EC$Longf <- genomic_USFQ$long_f
#genomic_EC$Symptom.onset.date <- genomic_USFQ$symptom_onset_date %>% dmy()
#genomic_EC$Travel.history <- genomic_USFQ$travel.history
#genomic_EC$Travel.history.dest <- genomic_USFQ$travel_history_destination
#genomic_EC$Notes <- genomic_USFQ$epidemiological_link
#genomic_EC$Indigenous.comm <- genomic_USFQ$indigenous_community


######### Downsample metadata data frames based on GISAID downsampler dictionaries ###########################################

# Systematic downsampling
intros_syst_subset <- intros[intros$seq %in% names(subset_seqs_syst),]
sources_syst_continent_subset <- sources_continent[sources_continent$seq %in% names(subset_seqs_syst),]
sources_syst_country_subset <- sources_country[sources_country$seq %in% names(subset_seqs_syst),]
dates_syst_subset <- dates[dates$seq %in% names(subset_seqs_syst),]

# Random downsampling
intros_rand_subset1 <- intros[intros$seq %in% names(subset_seqs_rand1),]
sources_rand_continent_subset1 <- sources_continent[sources_continent$seq %in% names(subset_seqs_rand1),]
sources_rand_country_subset1 <- sources_country[sources_country$seq %in% names(subset_seqs_rand1),]
dates_rand_subset1 <- dates[dates$seq %in% names(subset_seqs_rand1),]

intros_rand_subset2 <- intros[intros$seq %in% names(subset_seqs_rand2),]
sources_rand_continent_subset2 <- sources_continent[sources_continent$seq %in% names(subset_seqs_rand2),]
sources_rand_country_subset <- sources_country[sources_country$seq %in% names(subset_seqs_rand2),]
dates_rand_subset <- dates[dates$seq %in% names(subset_seqs_rand2),]

intros_rand_subset3 <- intros[intros$seq %in% names(subset_seqs_rand3),]
sources_rand_continent_subset3 <- sources_continent[sources_continent$seq %in% names(subset_seqs_rand3),]
sources_rand_country_subset3 <- sources_country[sources_country$seq %in% names(subset_seqs_rand3),]
dates_rand_subset3 <- dates[dates$seq %in% names(subset_seqs_rand3),]

######### Cleanup and data export #########################################################################################

write_tsv(rbind(intros_syst_subset, intros[grep("Ecuador", intros$seq),]), "Data/genetic_data/metadata/GISAID_hCoV-19_20210107_syst_subset_importation.tsv")
write_tsv(rbind(sources_syst_continent_subset, sources_continent[grep("Ecuador", sources_continent$seq),]), "Data/genetic_data/metadata/GISAID_hCoV-19_20210107_syst_subset_continent.tsv")
write_tsv(rbind(sources_syst_country_subset, sources_country[grep("Ecuador", sources_country$seq),]), "Data/genetic_data/metadata/GISAID_hCoV-19_20210107_syst_subset_country.tsv")
write_tsv(rbind(dates_syst_subset, dates[grep("Ecuador", dates$seq),]), "Data/genetic_data/metadata/GISAID_hCoV-19_20210107_syst_subset_dates.tsv")

rm(genomic_EC_manual)
