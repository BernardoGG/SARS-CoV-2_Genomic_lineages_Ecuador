#############################################################################################################################
######### SARS-CoV-2 Ecuador lineages plots #################################################################################
#############################################################################################################################

######### Bernardo Gutierrez ################################################################################################

source("GISAID_downsampler.R")
source("SARS2_phylo_metadata.R")
#load("COVID19_Ecuador_spatialdynamics_200420.RData") # <-- data from MSP report
#load("COVID19_Ecuador_spatialdynamicsplots_200420.RData") # <-- plots from MSP report
library('ggplot2')
library('ggpubr')
library('mapplots')
library('GADMTools')
library('sp')
library('RColorBrewer')
library('scatterpie')
library('patchwork')


######### Various plotting objects ##########################################################################################

## Generate maps for Ecuador & the Galapagos
map <- gadm_sp_loadCountries(fileNames = 'ECU', level = 1, basefile = './shapefiles', simplify = 0.01)

# Mainland Ecuador
map_ec <- gadm_subset(map, level = 1, regions = c("Carchi", "Esmeraldas", "Imbabura", "Sucumbios", "Pichincha",
                                                  "Santo Domingo de los Tsachilas", "Napo", "Manabi", "Orellana",
                                                  "Cotopaxi", "Tungurahua", "Los Rios", "Bolivar", "Pastaza", "Chimborazo",
                                                  "Guayas", "Santa Elena", "Cañar", "Morona Santiago", "Azuay", "El Oro",
                                                  "Loja", "Zamora Chinchipe"))
map_ec <- fortify(map_ec[[2]])

# Galapagos
map_gal <- gadm_subset(map, level = 1, regions = "Galápagos")
map_gal <- fortify(map_gal[[2]])


## Generate data frames
# Ecuacovid summarised epidemiological data extracted from MSP reports
ecuacovid <- read.csv("Data/MSP_epi_data/ecuacovid.csv", header = TRUE) # <-- file cloned on 2021-02-12
ecuacovid$created_at <- dmy(ecuacovid$created_at)

# OWID and GISAID cases and sequencing summary
latamseq <- read.csv("Data/OWID_GISAID_regional_data_20201210.csv", header = TRUE) # <-- data cutoff 2020-12-10
latamseq$country <- as.factor(latamseq$country)
latamseq$country <- factor(latamseq$country, levels = c("Suriname", "Guyana", "Uruguay", "Paraguay", "Venezuela",
                                                        "Bolivia", "Panama", "Ecuador", "Chile", "Peru", "Colombia",
                                                        "Argentina", "Brazil")) # <-- re-order

latamseq <- latamseq[order(latamseq$cases),]
latamseq$log_seqs_per_1k_deaths <- log(latamseq$seqs_per_1k_deaths)
latamseq$log_seqs_per_1k_deaths[latamseq$log_seqs_per_1k_deaths=="-Inf"] <- NA
conc <- vector()
seqs_per_deaths <- vector()
for(i in 1:nrow(latamseq)){
  conc <- rep(levels(latamseq$country)[i], round(latamseq$seqs_per_1k_deaths[i]))
  seqs_per_deaths <- c(seqs_per_deaths, conc)
}
rm(conc)
seqs_per_deaths <- as.data.frame(seqs_per_deaths)
seqs_per_deaths$seqs_per_deaths <- factor(seqs_per_deaths$seqs_per_deaths,
                                          levels = c("Suriname", "Guyana", "Uruguay", "Paraguay", "Venezuela",
                                                        "Bolivia", "Panama", "Ecuador", "Chile", "Peru", "Colombia",
                                                        "Argentina", "Brazil")) # <-- re-order

# Create transmission lineage summaries
genomic_EC$Cluster[genomic_EC$Cluster=="outgroup/singleton"] <- "singletons"
genomic_EC$Cluster <- factor(genomic_EC$Cluster, levels = c("V", "U", "T", "S", "R", "Q", "P", "O", "N", "M", "L",
                                                            "K", "J", "I", "H*", "G", "F", "E", "D*", "C", "B", "A", "singletons"))

trans_lineages_ec <- genomic_EC[genomic_EC$Province!="Galapagos",] %>% group_by(Cluster, Lat1, Long1) %>% dplyr::summarise(n = n()) %>%
  spread(Cluster, n, fill = 0)
trans_lineages_ec <- cbind(as.vector(rowSums(trans_lineages_ec[3:ncol(trans_lineages_ec)])), trans_lineages_ec)
colnames(trans_lineages_ec)[1] <- "Total"

trans_lineages_gal <- genomic_EC[genomic_EC$Province=="Galapagos",] %>% group_by(Cluster, Lat1, Long1) %>% dplyr::summarise(n = n()) %>%
  spread(Cluster, n, fill = 0)
trans_lineages_gal <- cbind(as.vector(rowSums(trans_lineages_gal[3:ncol(trans_lineages_gal)])), trans_lineages_gal)
colnames(trans_lineages_gal)[1] <- "Total"

# Number of sequences versus cumulative number of deaths
casenumbers_dec_ec <- data.frame(Province = as.vector(levels(genomic_EC$Province)),
                                 Cases = c(7732, 4935, 26080, 4707, 14061, 2420, 5744, 12670, 2463, 2576, 3399, 2986,
                                           5274, 5695, 7252, 72305, 6971, 3422, 1605, 2100, 2360, 2934, 1628, 861),
                                 Deaths = c(553+199, 237+54, 1851+1684, 361+242, 1257+1059, 383+280, 371+116, 219+20,
                                            71+14, 99+6, 114+2, 344+123, 303+65, 194+15, 245+47, 2015+273, 339+297,
                                            35+0, 75+2, 53+17, 61+17, 98+3, 55+1, 4+2)) # <- data taken from @Salud_Ec
                                                                                        #    Twitter feed from 2020-12-10

seqnumbers_ec <- as.data.frame(table(genomic_EC$Province))
colnames(seqnumbers_ec) <- c("Province", "Sequences")
seqs_cases_ec <- data.frame(Province = as.factor(casenumbers_dec_ec$Province),
                            Cases = as.numeric(casenumbers_dec_ec$Cases),
                            Deaths = as.numeric(casenumbers_dec_ec$Deaths),
                            Sequences = as.numeric(seqnumbers_ec$Sequences))
seqs_cases_ec$Seqs_per_1k_Deaths <- round((seqs_cases_ec$Sequences * 1000) / seqs_cases_ec$Deaths)

# TempEst tip to root divergence (heuristically rooted IQtree)
tempest <- read.csv("Data/genetic_data/TempEst_heuristic_GISAID_hCoV-19_20210107_syst_subset.tsv", sep = "\t")
tempest$loc <- ifelse(grepl("Ecuador", tempest$tip)==TRUE, "EC", "nonEC")

# Phylogenetically defined transmission lineages size, duration and spread
trans_lins_stats <- read.csv("Data/genetic_data/trans_lin_stats.csv", header = TRUE) %>% filter(category != "singleton")
trans_lins_stats$earliest <- as.Date(trans_lins_stats$earliest)
trans_lins_stats$latest <- as.Date(trans_lins_stats$latest)
trans_lins_stats$tmrca_median <- as.Date(trans_lins_stats$tmrca_median)
trans_lins_stats$tmrca_upper <- as.Date(trans_lins_stats$tmrca_upper)
trans_lins_stats$tmrca_lower <- as.Date(trans_lins_stats$tmrca_lower)
trans_lins_stats$duration <- as.numeric(trans_lins_stats$latest - trans_lins_stats$earliest)
trans_lins_stats$lag <- as.numeric(trans_lins_stats$earliest - trans_lins_stats$tmrca_median)

# Phylogenetically defined transmission lineages size, duration and spread
singletons <- genomic_EC %>% filter(Cluster == "singletons")
singletons$Province <- factor(singletons$Province, levels = c("El Oro", "Guayas", "Los Rios", "Manabi", "Santa Elena", "Santo Domingo",
                                                              "Imbabura", "Loja", "Pichincha", "Tungurahua",
                                                              "Morona Santiago", "Pastaza", "Galapagos"))

## Create colour schemes
# Colour by province palette
genomic_EC$Province <- factor(genomic_EC$Province, levels = c("El Oro", "Esmeraldas", "Guayas", "Los Rios", "Manabi",
                                                              "Santa Elena", "Santo Domingo",
                                                              "Azuay", "Bolivar", "Canar", "Carchi", "Chimborazo", "Cotopaxi",
                                                              "Imbabura", "Loja", "Pichincha", "Tungurahua",
                                                              "Morona Santiago", "Napo", "Orellana", "Pastaza", "Sucumbios",
                                                              "Zamora Chinchipe",
                                                              "Galapagos")) # <-- re-order provinces per region
genomic_EC$Region <- factor(genomic_EC$Region, levels = c("Coast", "Highlands", "Amazon", "Galapagos")) # <-- re-order regions

ECColors <- c("#095D53", "#0B6F64", "#0D8274", "#0F9585", "#11A796", "#12BAA6", "#14CCB7",
              "#513A01", "#795701", "#A27402", "#CA9102", "#F2AE02", "#FDBE21", "#FDCA49", "#FED672", "#FEE29A", "#FEEDC2",
              "#91540C", "#BD690F", "#E27E12", "#EE922F", "#F1A655", "#F4BA7B",
              "#0F5077")
names(ECColors) <- levels(genomic_EC$Province)

# Colours for provinces in map
map_colors_ec <- data.frame(prov = as.vector(levels(genomic_EC$Province)[1:length(levels(genomic_EC$Province))-1]),
                            map = c("ECU_22", "ECU_23", "ECU_2", "ECU_5", "ECU_6", "ECU_13", "ECU_14", "ECU_1", "ECU_12", "ECU_18", "ECU_19",
                            "ECU_20", "ECU_21", "ECU_3", "ECU_4", "ECU_11", "ECU_16", "ECU_7", "ECU_8", "ECU_9", "ECU_10", "ECU_15",
                            "ECU_17"))

map_ec$id <- mapvalues(map_ec$id, from = map_colors_ec$map, to = map_colors_ec$prov)

# Colour by region palette
RegionColors <- c("#0F9585", "#F2AE02", "#E27E12", "#0F5077")
names(RegionColors) <- levels(genomic_EC$Region)

# Colour by global lineage palette
lin.cols <- length(levels(genomic_EC$Lineage))
LineageColors <- colorRampPalette(brewer.pal(8, "Accent"))(lin.cols)

# Colour by transmission lineage palette
clust.cols <- length(levels(genomic_EC$Cluster))-1
ClusterColors <- colorRampPalette(brewer.pal(8, "Blues"))(clust.cols)
ClusterColors <- c(ClusterColors, "#C6AF88")


######### FIGURE 1 - Sampling and sequencing ##########################################################################

## Number of sequences versus cumulative number of deaths
fig1a <- ggplot(seqs_cases_ec, aes(x = Deaths, y = Sequences, size = Cases, color = Province)) + geom_point(alpha = 0.7) +
  theme_light() + scale_color_manual(values = ECColors) + scale_size_area(max_size = 30) +
  xlab("Cumulative number of confirmed deaths") + ylab("Number of sequences") # <-- confirm this plot through a Spearman correlation


## Sampling across time and provinces
sampling_start <- min(genomic_EC$Collection.date)
sampling_end <- max(genomic_EC$Collection.date)

fig1b <- ggplot(genomic_EC, aes(x = Collection.date, fill = Province)) +
  geom_dotplot(method = "histodot", stackgroups = TRUE, binpositions = "all",
               binwidth = length(unique(genomic_EC$Collection.date))/round((as.numeric(sampling_end - sampling_start)/14))) +
  theme_light() + scale_fill_manual(values = ECColors) + xlab("Collection date") + ylab("No. of samples (per 2-week periods periods)") +
  scale_y_continuous(breaks = c(0:1))

ecuacovid2020 <- ecuacovid[ecuacovid$created_at<sampling_end,]
ecuacovid2020$defunciones_promedio <- as.vector((ecuacovid2020$defunciones_2015_nuevas + ecuacovid2020$defunciones_2016_nuevas +
                                                       ecuacovid2020$defunciones_2017_nuevas + ecuacovid2020$defunciones_2018_nuevas +
                                                       ecuacovid2020$defunciones_2019_nuevas)/5)

fig1b2 <- ggplot(ecuacovid2020, aes(x = created_at, y = defunciones_2020_nuevas - defunciones_promedio)) +
  geom_col(fill = "grey")+ geom_line(data = ecuacovid2020, aes(x = created_at, y = positivas_pcr/17643060*50000)) +
  theme_light() + xlab("Collection date") + ylab("No. of daily excess deaths/Incidence per 50 000 people")


## Global lineages identified and geographical representation
global_lineages_ec <- genomic_EC %>% group_by(Lineage, Region, Province) %>% dplyr::summarise(n = n())
global_lineages_ec$Region <- factor(global_lineages_ec$Region, levels = c("Galapagos", "Amazon", "Highlands", "Coast"))

fig1c <- ggplot(global_lineages_ec, aes(x = Lineage, y = n, fill = Province, group = Region)) + geom_bar(stat = "identity") +
  theme_light() + theme(axis.text.x=element_text(angle=30, hjust=1), legend.position = "bottom") + scale_fill_manual(values = ECColors) +
  xlab("Global lineage") + ylab("No. of sequences")


## Map of Ecuador provinces and number of sequences
fig1d <- ggplot() + geom_polygon(data = map_ec, aes(x = long, y = lat, group = group, fill = id), color = '#FFFFFF', size = 0.4) +
  ggforce::geom_circle(data = trans_lineages_ec, aes(x0 = Long1, y0 = Lat1, r = Total/110), fill = '#999999', size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = ECColors) + theme_minimal() + theme(line = element_blank(), axis.text = element_blank(),
                                                                 axis.title = element_blank(), legend.position = "none")

fig1d_gal <- ggplot() + geom_polygon(data = map_gal, aes(x = long, y = lat, group = group), color = '#FFFFFF', fill = ECColors[names(ECColors)=="Galapagos"], size = 0.4) +
  ggforce::geom_circle(data = trans_lineages_gal, aes(x0 = Long1, y0 = Lat1, r = Total/110), fill = '#999999', size = 0.4, alpha = 0.9) +
  theme_minimal() + theme(line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "none")


######### FIGURE 2 - Transmission lineages correlations and trees ##############################################################

## Plot tip-to-root divergence
fig2b <- ggplot(tempest, aes(x = date, y = distance)) + geom_point(aes(color = loc)) + geom_smooth(method = lm, color = '#14748F')  +
  scale_color_manual(values = c('#8b4d4d', '#8e8d8d')) +
  theme_light() + theme(legend.position = "none") + xlab("Sample collection date") + ylab("Tip-to-root divergence")


## Plot pairwise comparison of transmission lineage sampling times, persistence and geographic spread
# Timeline for transmission lineage TMRCAs and first detections (earliest sequences)
trans_lins_lag <- trans_lins_stats %>% select(-category, -size, -latest, -tmrca_upper, -tmrca_lower, -num_provs, -duration, -lag) %>%
  gather("event", "date", -trans_lin)
trans_lins_lag$event <- as.factor(trans_lins_lag$event)
levels(trans_lins_lag$event) <- c("Transmission lineage first detection", "Transmission lineage TMRCA")

fig2c <- ggplot(trans_lins_lag) + geom_line(aes(x = date, y = event, group = trans_lin, color = '#d3d5d6')) +
  geom_point(aes(x = date, y = event, color = event), size = 3.0) + scale_color_manual(values = c('#d3d5d6', '#8b4d4d','#14748f')) +
  theme_light() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none",
                        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  scale_x_date(date_breaks = "months", date_labels = "%b-%y") + coord_flip()


######### FIGURE 3 - Overview of transmission lineages ##########################################################################

## Earliest collection date and TMRCA vs size, persistence and spread
# Earliest collection date vs geographic spread
fig3a_1 <- ggplot(trans_lins_stats, aes(x = earliest, y = num_provs)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') +
  theme_light() + xlab("Collection date of earliest sequence in transmission lineage") + ylab("Number of provinces where found")

# TMRCA vs geographic spread
fig3a_2 <- ggplot(trans_lins_stats, aes(x = tmrca_median, y = num_provs)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') +
  theme_light() + xlab("Transmission lineage median TMRCA") + ylab("Number of provinces where found")

fig3a <- (fig3a_1 / fig3a_2)

## Plot maps with sampled proportion of transmission lineages found in each province
# Fixed pie chart sizes
fig3b <- ggplot() + geom_polygon(data = map_ec, aes(x = long, y = lat, group = group), fill = '#999999', color = '#d3d5d6', size = 0.7) +
  geom_scatterpie(aes(x = Long1, y = Lat1, r = 0.2), data = trans_lineages_ec, cols = colnames(trans_lineages_ec[4:ncol(trans_lineages_ec)]),
                  alpha = 0.7) + scale_fill_manual(values = ClusterColors) +
  coord_fixed() + theme_minimal() + theme(line = element_blank(), axis.text = element_blank(), axis.title = element_blank())

fig3b_gal <- ggplot() + geom_polygon(data = map_gal, aes(x = long, y = lat, group = group), fill = '#999999', color = '#d3d5d6', size = 0.7) +
  geom_scatterpie(aes(x = Long1, y = Lat1, r = 0.2), data = trans_lineages_gal, cols = colnames(trans_lineages_gal[4:ncol(trans_lineages_gal)]),
                  alpha = 0.7) + scale_fill_manual(values = ClusterColors) +
  coord_fixed() + theme_minimal() + theme(line = element_blank(), axis.text = element_blank(), axis.title = element_blank())

## Geographic representation per transmission lineage
clusters_ec <- genomic_EC %>% group_by(Cluster, Region, Province) %>% dplyr::summarise(n = n())
clusters_ec$Region <- factor(clusters_ec$Region, levels = c("Galapagos", "Amazon", "Highlands", "Coast"))

fig3c <- ggplot(clusters_ec, aes(x = Cluster, y = n, fill = Province, group = Region)) + geom_bar(stat = "identity") +
  theme_light() + coord_flip() + scale_fill_manual(values = ECColors) +
  xlab("Transmission lineage") + ylab("No. of sequences")


######### SUPPLEMENTARY FIGURES ################################################################################################

## Fig S1 - Overview of LATAM sequencing
latamseq1 <- latamseq %>% filter(country!="Uruguay" & country!="Suriname" & country!="Guyana")
latamseq2 <- latamseq %>% filter(country=="Uruguay" | country=="Suriname" | country=="Guyana")

col.range <- c(0, max(latamseq$log_seqs_per_1k_deaths, na.rm = TRUE))

figS1a <- ggplot(latamseq1, aes(x = country, y = cases, fill = log_seqs_per_1k_deaths)) + geom_col() + theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
  scale_fill_gradient(limits = col.range) + coord_flip()

figS1b <- ggplot(latamseq1, aes(x = country, y = deaths, fill = log(seqs_per_1k_deaths))) + geom_col() +
  theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.title.x = element_blank(), legend.position = "none") + scale_fill_gradient(limits = col.range) + coord_flip()

figS1c <- ggplot(latamseq1, aes(x = factor(country), y = round(seqs_per_1k_deaths), fill = log(seqs_per_1k_deaths))) +
  geom_linerange(ymin = 0, ymax = round(latamseq1$seqs_per_1k_deaths)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.8) +
  theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.title.x = element_blank()) + scale_fill_gradient(limits = col.range, name = "") +
  coord_flip()

figS1d <- ggplot(latamseq2, aes(x = country, y = cases, fill = log(seqs_per_1k_deaths))) + geom_col() + theme_classic() +
  theme(axis.title.y = element_blank(), legend.position = "none") + scale_color_gradient(limits = col.range) +
  scale_fill_gradient(limits = col.range) + coord_flip() + xlab("Country") + ylab("No. of cases")

figS1e <- ggplot(latamseq2, aes(x = country, y = deaths, fill = log(seqs_per_1k_deaths))) + geom_col() +
  theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  scale_fill_gradient(limits = col.range) + coord_flip() + ylab("No. of deaths")

figS1f <- ggplot(latamseq2, aes(x = factor(country), y = round(seqs_per_1k_deaths), fill = log(seqs_per_1k_deaths))) +
  geom_linerange(ymin = 0, ymax = round(latamseq2$seqs_per_1k_deaths)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.8) +
  theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  scale_fill_gradient(limits = col.range) + ylab("No. of sequences per 1000 deaths") + coord_flip()
  
figS1 <- ((figS1a | figS1b | figS1c) / (figS1d | figS1e | figS1f)) + plot_layout(heights = c(11, 3), guides = "collect")

## Fig S2 - Overview of singleton distribution across time
figS2_1 <- ggplot(singletons) + geom_histogram(aes(Collection.date), binwidth = 14, fill = '#14748F', color = 'white') +
  scale_x_date(date_breaks = "months", date_labels = "%b-%y") + theme_light() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

figS2_2 <- ggplot(singletons, aes(x = Collection.date, y = Lineage, color = Province)) + geom_point(size = 2.0) +
  scale_color_manual(values = ECColors) + scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
  theme_light() + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),
                        axis.ticks.x = element_blank()) + xlab("Singleton sample collection date")

figS2_3 <- ggplot(singletons, aes(x = Collection.date, y = Province, color = Province)) + geom_point(size = 2.0) +
  scale_color_manual(values = ECColors) + scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
  theme_light() + theme(legend.position = "none") + xlab("Singleton sample collection date")

figS2 <- (figS2_1 / figS2_2 / figS2_3)

## Fig S10 - Timelines for transmission lineage lag times
figS10_1 <- ggplot(trans_lins_stats) + geom_histogram(aes(lag), binwidth = 5, fill = '#14748F', color = "grey", size = 0.2) +
  theme_light() + theme(axis.title.y = element_blank()) + xlab("Detection lag (days)")

figS10_2 <- ggplot(trans_lins_stats, aes(x = earliest, y = lag)) + geom_smooth(method = glm, color = '#14748F') +
  geom_point(size = 5.0, color = '#8b4d4d') + ylim(0,max(trans_lins_stats$lag)) +
  theme_light() + xlab("Collection date of earliest sequence in transmission lineage") +
  ylab("Detection lag (days)")

figS10 <- figS10_2 + inset_element(figS10_1, left = 0.08, bottom = 0.6, right = 0.5, top = 0.97, align_to = 'full')

## Fig S11 - Transmission lineage persistence vs size
figS11_1 <- ggplot(trans_lins_stats, aes(x = duration, y = size)) + geom_point(size = 5.0, color = '#8b4d4d') +
  geom_smooth(method = lm, color = '#14748F') +
  theme_light() + xlab("Transmission lineage \n sampling period (no. of days)") + ylab("Transmission lineage size (no. of sequences)")

figS11_2 <- ggplot(trans_lins_stats, aes(x = duration + lag, y = size)) + geom_point(size = 5.0, color = '#8b4d4d') +
  geom_smooth(method = lm, color = '#14748F') + theme_light() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                                                                      axis.ticks.y = element_blank()) +
  xlab("Transmission lineage sampling period + \n detection lag (no. of days)")

figS11 <- (figS11_1 | figS11_2)

## Fig S12 - Earliest collection date and TMRCA vs transmission lineage size
figS12_1 <- ggplot(trans_lins_stats, aes(x = earliest, y = size)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') +
  theme_light() + xlab("Collection date of earliest sequence in transmission lineage") + ylab("Transmission lineage size\n(no. of sequences)")

figS12_2 <- ggplot(trans_lins_stats, aes(x = tmrca_median, y = size)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') +
  theme_light() + xlab("Transmission lineage median TMRCA") +
  ylab("Transmission lineage size\n(no. of sequences)")

figS12 <- (figS12_1 / figS12_2)

## Fig S13 - Earliest collection date and TMRCA vs transmission lineage persistance
figS13_2 <- ggplot(trans_lins_stats, aes(x = tmrca_median, y = duration)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') + 
  theme_light() + xlab("Transmission lineage median TMRCA") +
  ylab("Transmission lineage\npersistence (no. of days)")

figS13_1 <- ggplot(trans_lins_stats, aes(x = earliest, y = duration)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') +
  theme_light() + xlab("Collection date of earliest sequence in transmission lineage") +
  ylab("Transmission lineage\npersistence(no. of days)")

figS13 <- (figS13_1 / figS13_2)

## Fig S14 - Earliest collection date and TMRCA vs transmission lineage persistance
figS14_2 <- ggplot(trans_lins_stats, aes(x = tmrca_median, y = num_provs)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') + 
  theme_light() + xlab("Transmission lineage median TMRCA") +
  ylab("Number of provinces\nwhere found")

figS14_1 <- ggplot(trans_lins_stats, aes(x = earliest, y = num_provs)) + geom_smooth(method = lm, color = '#14748F', se = FALSE) +
  geom_smooth(method = loess, color = '#0A3A47', se = FALSE) + geom_point(size = 5.0, color = '#8b4d4d') +
  theme_light() + xlab("Collection date of earliest sequence in transmission lineage") +
  ylab("Number of provinces\nwhere found")

figS14 <- (figS14_1 / figS14_2)

## Fig S15 - Pie chart sizes proportional to number of sequences
figS15 <- ggplot() + geom_polygon(data = map_ec, aes(x = long, y = lat, group = group), fill = '#999999', color = '#d3d5d6', size = 0.7) +
  geom_scatterpie(aes(x = Long1, y = Lat1, r = Total/72), data = trans_lineages_ec, cols = colnames(trans_lineages_ec[4:ncol(trans_lineages_ec)]),
                  alpha = 0.7) + scale_fill_manual(values = ClusterColors) +
  coord_fixed() + theme_minimal() + theme(line = element_blank(), axis.text = element_blank(), axis.title = element_blank())

## Fig S16 - B.1.1.74 sampling and frequency timeline
figS16_1 <- ggplot(genomic_EC[genomic_EC$Lineage=="B.1.1.74",]) +
  geom_histogram(aes(Collection.date), binwidth = 7, fill = '#14748F', color = 'white') +
  scale_x_date(date_breaks = "months", date_labels = "%b-%y") + scale_y_continuous(breaks = c(2,4,6,8)) +
  theme_light() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())

figS16_2 <- ggplot(genomic_EC[genomic_EC$Lineage=="B.1.1.74",], aes(x = Collection.date, y = Province, color = Province)) +
              geom_point(size = 2.0) + scale_color_manual(values = ECColors) + scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
              theme_light() + theme(legend.position = "none") + xlab("Sample collection date")

figS16 <- (figS16_1 / figS16_2) + plot_annotation(title = "B.1.1.74")
  
######### Export figures #######################################################################################################

ggsave("Plots/fig1a.pdf", plot = fig1a, width = 11.0, height = 10.38, units = 'in', dpi = 600)
ggsave("Plots/fig1b.pdf", plot = fig1b, dpi = 600)
ggsave("Plots/fig1b2.pdf", plot = fig1b2, dpi = 600)
ggsave("Plots/fig1c.pdf", plot = fig1c, dpi = 600)
ggsave("Plots/fig1d.pdf", plot = fig1d, width = 9.65, height = 10.38, units = 'in', scale = 2, dpi = 600)
ggsave("Plots/fig1d_gal.pdf", plot = fig1d_gal, width = 9.65, height = 10.38, units = 'in', scale = 2, dpi = 600)

ggsave("Plots/fig2b.pdf", plot = fig2b, width = 11.0, height = 10.38, units = 'in', dpi = 600)
ggsave("Plots/fig2c.pdf", plot = fig2c, width = 5.5, height = 10.38, units = 'in', dpi = 600)

ggsave("Plots/fig3a.pdf", plot = fig3a, width = 4.00, height = 7.38, units = 'in', scale = 2, dpi = 600)
ggsave("Plots/fig3b.pdf", plot = fig3b, width = 9.65, height = 10.38, units = 'in', scale = 2, dpi = 600)
ggsave("Plots/fig3b_gal.pdf", plot = fig3b_gal, width = 9.65, height = 10.38, units = 'in', scale = 2, dpi = 600)
ggsave("Plots/fig3c.pdf", plot = fig3c, dpi = 600)

ggsave("Plots/figS1.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS1, dpi = 600)
ggsave("Plots/figS2.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS2, dpi = 600)
ggsave("Plots/figS10.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS10, dpi = 600)
ggsave("Plots/figS11.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS11, dpi = 600)
ggsave("Plots/figS12.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS12, dpi = 600)
ggsave("Plots/figS13.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS13, dpi = 600)
ggsave("Plots/figS14.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS14, dpi = 600)
ggsave("Plots/figS15.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS15, dpi = 600)
ggsave("Plots/figS16.pdf", width = 280.00, height = 150.00, units = 'mm', scale = 0.72, plot = figS16, dpi = 600)

write_csv(latamseq, "Analyses/File_S1.csv")
write_csv(seqs_cases_ec, "Analyses/File_S2.csv")
#-------------------------------------------------------- END -----------------------------------------------------------------