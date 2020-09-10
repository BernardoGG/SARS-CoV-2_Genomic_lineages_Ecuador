###################################
#### SARS-CoV-2 sequence stats ####
###################################

####### Bernardo Gutierrez ########

library("ape")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("stringr")

## Function to clean sequence location metadata "cleanlocs"
# x = Sequence metadata data frame (must contain a "location" column)
# y = Dictionary dataframe containing old and new "location" levels
locvector <- vector()
cleanlocs <- function(x, y){
  for(i in 1:length(x$location)){
    for(j in 1:length(y$old))
      if (x$location[i]==y$old[j]) {
        locvector[i] = as.character(y$new[j])
        break
      }
  }
  print(locvector)
}

##################
###### Data ######
##################

### Sequences at two data points
## Criteria: complete (>29kbp), high coverage (<1% Ns; <0.05% unique AA mutations not seen in other seqs in db)
## and low coverage excluded (>5% Ns).

## 1. Sequences overlapping with current travel history data
seqs_may <- read.FASTA("gisaid_hcov-19_up_to_2020_05_01.fasta") # <-- GISAID seqs collected up to 2020/05/01

labs_may <- names(seqs_may)
labs_may <- gsub("/", "|", labs_may) %>% str_split_fixed(pattern = fixed("|"), 6) %>% as.data.frame()
colnames(labs_may) <- c("virus", "location", "lab_id", "year", "EPI_ISL", "collection_date")

## Clean up sequence locations
temp <- as.vector(levels(labs_may$location))

# Create dictionary dataframe
lch <- c(grep("Wuhan", temp), grep("Beijing", temp), grep("Changzhou", temp), grep("Chongqing", temp), grep("Foshan", temp),
         grep("Fujian", temp), grep("Fuyang", temp), grep("Fuzhou", temp), grep("Ganzhou", temp), grep("Guangdong", temp),
         grep("Guangzhou", temp), grep("Hangzhou", temp), grep("Harbin", temp), grep("Hefei", temp), grep("Henan", temp),
         grep("Jian", temp), grep("Jiangsu", temp), grep("Jiangxi", temp), grep("Jingzhou", temp), grep("Jiujiang", temp),
         grep("Lishui", temp), grep("Nanchang", temp), grep("NanChang", temp), grep("Pingxiang", temp), grep("Shandong", temp),
         grep("Shanghai", temp), grep("Shangrao", temp), grep("Shaoxing", temp), grep("Shenzhen", temp), grep("Sichuan", temp),
         grep("Tianmen", temp), grep("Weifang", temp), grep("Xinyu", temp), grep("Yichun", temp), grep("Yingtan", temp),
         grep("Yunnan", temp), grep("Zhejiang", temp))
temp[lch] <- "China"
lro <- c(grep("Romania ", temp), grep("Bucuresti", temp))
temp[lro] <- "Romania"
luk <- c(grep("England", temp), grep("Scotland", temp), grep("Northern Ireland", temp), grep("Wales", temp))
temp[luk] <- "UK"
lsa <- grep("SaudiArabia", temp)
temp[lsa] <- "Saudi Arabia"
lba <- grep("Bahrein", temp)
temp[lba] <- "Bahrain"
lex <- c(grep("bat", temp), grep("env", temp), grep("Felis catus", temp), grep("mink", temp),
         grep("pangolin", temp), grep("tiger", temp), grep("canine", temp))
temp[lex] <- "EXCLUDE"

ldf <- data.frame(old = as.vector(levels(labs_may$location)), new = temp) # <-- dictionary data frame
labs_may$location <- cleanlocs(labs_may, ldf) # <- apply custom function to clean "location" column (long runtime)
locations_may <- as.data.frame(table(labs_may$location)) %>% arrange(Var1) # <-- summarise sequence locations

## 2. All high quality sequences (as defined by the "Criteria")
seqs_sept <- read.FASTA("gisaid_hcov-19_2020_09_08.fasta") # <-- GISAID seqs submitted up to 2020/09/08

labs_sept <- names(seqs_sept)
labs_sept <- gsub("/", "|", labs_sept) %>% str_split_fixed(pattern = fixed("|"), 6) %>% as.data.frame()
colnames(labs_sept) <- c("virus", "location", "lab_id", "year", "EPI_ISL", "collection_date")

## Clean up sequence locations
temp <- as.vector(levels(labs_sept$location))

# Create dictionary dataframe
lch <- c(grep("Wuhan", temp), grep("Beijing", temp), grep("Changzhou", temp), grep("Chongqing", temp), grep("Foshan", temp),
         grep("Fujian", temp), grep("Fuyang", temp), grep("Fuzhou", temp), grep("Ganzhou", temp), grep("Guangdong", temp),
         grep("Guangzhou", temp), grep("Hangzhou", temp), grep("Harbin", temp), grep("Hefei", temp), grep("Henan", temp),
         grep("Jian", temp), grep("Jiangsu", temp), grep("Jiangxi", temp), grep("Jingzhou", temp), grep("Jiujiang", temp),
         grep("Lishui", temp), grep("Nanchang", temp), grep("NanChang", temp), grep("Pingxiang", temp), grep("Shandong", temp),
         grep("Shanghai", temp), grep("Shangrao", temp), grep("Shaoxing", temp), grep("Shenzhen", temp), grep("Sichuan", temp),
         grep("Tianmen", temp), grep("Weifang", temp), grep("Xinyu", temp), grep("Yichun", temp), grep("Yingtan", temp),
         grep("Yunnan", temp), grep("Zhejiang", temp), grep("Liaoning", temp))
temp[lch] <- "China"
lro <- c(grep("Romania ", temp), grep("Bucuresti", temp))
temp[lro] <- "Romania"
luk <- c(grep("England", temp), grep("Scotland", temp), grep("Northern Ireland", temp), grep("Wales", temp))
temp[luk] <- "UK"
lsa <- grep("SaudiArabia", temp)
temp[lsa] <- "Saudi Arabia"
lba <- grep("Bahrein", temp)
temp[lba] <- "Bahrain"
lex <- c(grep("bat", temp), grep("env", temp), grep("Felis catus", temp), grep("mink", temp),
         grep("pangolin", temp), grep("tiger", temp), grep("canine", temp), grep("cat", temp))
temp[lex] <- "EXCLUDE"

ldf <- data.frame(old = as.vector(levels(labs_sept$location)), new = temp) # <-- dictionary data frame
labs_sept$location <- cleanlocs(labs_sept, ldf) # <- apply custom function to clean "location" column (long runtime)
locations_sept <- as.data.frame(table(labs_sept$location)) %>% arrange(Var1) # <-- summarise sequence locations

### Cumulative number of cases at two data points
## Data taken from Our World in Data, includes cases reported between 2020-05-01 and 2020-09-08.
cases_owid <- read.csv("owid-covid-data_may_sept.csv")
cases_owid$location <- as.character(cases_owid$location)
cases_owid$location[cases_owid$location=="Democratic Republic of Congo"] <- "DRC"
cases_owid$location[cases_owid$location=="Faeroe Islands"] <- "Faroe Islands"
cases_owid$location[cases_owid$location=="Macedonia"] <- "North Macedonia"
cases_owid$location[cases_owid$location=="Timor"] <- "Timor-Leste"
cases_owid$location[cases_owid$location=="United Kingdom"] <- "UK"
cases_owid$location[cases_owid$location=="United States"] <- "USA"

cases_may <- cases_owid[cases_owid$date=="2020-05-01",] %>% arrange(location) %>%
  select(location, date, total_cases) %>% filter(location %in% levels(locations_may$Var1))

cases_sept <- cases_owid[cases_owid$date=="2020-09-08",] %>% arrange(location) %>%
  select(location, date, total_cases) %>% filter(location %in% levels(locations_sept$Var1))

### Cases versus sequences data frame
counts_may <- data.frame(country = cases_may$location,
                          may_cases = cases_may$total_cases,
                          may_seqs = select(filter(locations_may, Var1 %in% cases_may$location), Freq))
colnames(counts_may)[3] <- "may_seqs"

counts_sept <- data.frame(country = cases_sept$location,
                         sept_cases = cases_sept$total_cases,
                         sept_seqs = select(filter(locations_sept, Var1 %in% cases_sept$location), Freq))
colnames(counts_sept)[3] <- "sept_seqs"

######################
###### Analyses ######
######################

### Plot number of cases versus number of sequences on both dates, colourdes by travek history to Ecuador
## Countries reported in travel histories of Ecuadorian patients
travel_histories_EC <- c("Argentina", "Belgium", "Brazil", "Chile", "China", "Colombia", "Costa Rica", "Cuba",
                         "France", "Egypt", "Germany", "Guatemala", "Italy", "Mexico", "Netherlands", "Panama",
                         "Peru", "Puerto Rico", "Spain", "Sweden", "USA", "UK", "Venezuela")

## Generate plots
ggplot() +
  geom_point(data = counts_may, aes(x = log(may_cases), y = log(may_seqs), color = "No travel history reported")) +
  geom_point(data = filter(counts_may, country %in% travel_histories_EC),
             aes(x = log(may_cases), y = log(may_seqs), color = country)) +
  theme_light() + labs(title = "Cases versus sequences - May 1", color = "Reported travel history to Ecuador",
                       x = "Number of cases (log)", y = "Number of sequences (log)") +
  geom_vline(xintercept = log(counts_may$may_cases[counts_may$country=="Ecuador"]), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = log(counts_may$may_seqs[counts_may$country=="Ecuador"]), linetype = "dashed", color = "gray")

ggplot() +
  geom_point(data = counts_sept, aes(x = log(sept_cases), y = log(sept_seqs), color = "No travel history reported")) +
  geom_point(data = filter(counts_sept, country %in% travel_histories_EC),
             aes(x = log(sept_cases), y = log(sept_seqs), color = country)) +
  theme_light() + labs(title = "Cases versus sequences - September 8", color = "Reported travel history to Ecuador",
                       x = "Number of cases (log)", y = "Number of sequences (log)") +
  geom_vline(xintercept = log(counts_sept$sept_cases[counts_sept$country=="Ecuador"]), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = log(counts_sept$sept_seqs[counts_sept$country=="Ecuador"]), linetype = "dashed", color = "gray")

### Create tables for sequences/cases ratios and expected taxa given specific threshold ratio values
## May dataset
may_taxa_est <- data.frame(country = counts_may$country,
                           ratio = counts_may$may_seqs/counts_may$may_cases,
                           sequences = counts_may$may_seqs,
                           target_taxa_0.005 = round(0.005*counts_may$may_cases),
                           diff_0.005 = round(0.005*counts_may$may_cases) - counts_may$may_seqs,
                           target_taxa_0.008 = round(0.008*counts_may$may_cases),
                           diff_0.008 = round(0.008*counts_may$may_cases) - counts_may$may_seqs) %>%
  arrange(-ratio)

nrow(labs_may) + sum(may_taxa_est$diff_0.005) #<-- dataset size at seq/case ratio of 0.005 
nrow(labs_may) + sum(may_taxa_est$diff_0.008) #<-- dataset size at seq/case ratio of 0.008

## September dataset
sept_taxa_est <- data.frame(country = counts_sept$country,
                           ratio = counts_sept$sept_seqs/counts_sept$sept_cases,
                           sequences = counts_sept$sept_seqs,
                           target_taxa_0.005 = round(0.005*counts_sept$sept_cases),
                           diff_0.005 = round(0.005*counts_sept$sept_cases) - counts_sept$sept_seqs,
                           target_taxa_0.008 = round(0.008*counts_sept$sept_cases),
                           diff_0.008 = round(0.008*counts_sept$sept_cases) - counts_sept$sept_seqs) %>%
  arrange(-ratio)

nrow(labs_sept) + sum(sept_taxa_est$diff_0.005) #<-- dataset size at seq/case ratio of 0.005 
nrow(labs_sept) + sum(sept_taxa_est$diff_0.008) #<-- dataset size at seq/case ratio of 0.008
