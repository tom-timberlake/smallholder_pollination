###############################################################################################################################################
##### Paper title: Agricultural specialisation increases the vulnerability of pollination services for smallholder farmers
##### Journal: Journal of Applied ecology   
##### Date: 2024
###############################################################################################################################################

#The following script provides the code for all figures and data analyses included in the paper

#####################################################################
## Setting up environment and loading files and packages
#####################################################################

##Load required packages
library(readxl)
library(data.table)
library(skimr)
library(openxlsx)
library(ggalluvial)
library(ggplot2)
library(gridExtra)
library(cowplot)
library("ggpubr")
library(reshape2)
library(Rmisc) 
library(vegan)
library(ggrepel)
library(tidyverse)
library(ggthemes)
library(scales)
library(bipartite)
library(svglite)


##Set a working directory
setwd("XXXX")  #Edit this before running

# Paths (Edit these before running)
input.path <- "XXXX" # Change this to the file path where source data is stored
output.path <- "XXXX" # Change this to the file path where results and figures will be saved


#Load in data files

##Pollinator visitation data - this is the main plant-pollinator interaction dataset from the visitation survey data 
#Download this file from Dryad: https://doi.org/10.5061/dryad.0rxwdbs91    
poll_visitation_data <- data.table(read_excel(file.path(input.path, "Plant_pollinator_visitation_data.xlsx"),
                                              sheet = "Visitation data", 
                                              na = c("", "---", NA)))  

##Import pollen capacity data - this is a database listing the mean number of pollen grains found on each insect taxon (termed pollen carrying capacity) 
#File available from GitHub: https://github.com/tom-timberlake/smallholder_pollination 
pollen_data <- read.csv(file.path(input.path, "pollen_capacity_OTU.csv")) 

##Import crop importance data - this is a dataset of each pollinator-dependent crop, listing their nutritional and economic importance values  
#File available from GitHub: https://github.com/tom-timberlake/smallholder_pollination 
crop_importance_data <- data.table(read_excel(file.path(input.path, "Jumla_crop_importance_values.xlsx"),
                                              sheet = "Sheet1", 
                                              na = c("", "---", NA)))  


######################################################################################################################
##################           Step 1: Calculating importance of each pollinator by crop      ##########################
######################################################################################################################

# Merge visitation and pollen carrying capacity data so that each insect-plant interaction has associated pollen data for that insect
poll_visit_pollen_data <- merge(x=poll_visitation_data, y=pollen_data,by="insect_OTU",all.x=TRUE)

#Create new variable called visits and assign each interaction as 1 visit
poll_visit_pollen_data$visits <- 1

##Create new column with only plant genus name
poll_visit_pollen_data$plant_genus <- gsub( " .*$", "", poll_visitation_data$plant_sci_name) 

###Calculate the number of visits made by each pollinator to each plant species, as a proportion of all visits the plant has received.
# This is done separately using visitation-only data, but also using pollen-weighted visitation data (i.e visitation * pollen-carrying-capacity)

# Summarising at the Operational Taxonomic Unit (OTU) level
total_plant_visits <- poll_visit_pollen_data %>%   group_by(plant_sci_name) %>%  dplyr::summarise(total_plant_visits= sum(visits, na.rm=TRUE), total_plant_pollen= sum(mean_pollen_load, na.rm=TRUE))
plant_poll_visits_OTU_village <- poll_visit_pollen_data %>%   group_by(insect_OTU, plant_sci_name) %>%  summarise(plant_poll_visits= sum(visits, na.rm=TRUE), plant_poll_pollen= sum(mean_pollen_load, na.rm=TRUE))
poll_plant_prop_visits_OTU <- merge(x=total_plant_visits, y=plant_poll_visits_OTU_village,by="plant_sci_name",all.x=TRUE)
poll_plant_prop_visits_OTU$prop_visits <- poll_plant_prop_visits_OTU$plant_poll_visits / poll_plant_prop_visits_OTU$total_plant_visits
poll_plant_prop_visits_OTU$prop_pollen <- poll_plant_prop_visits_OTU$plant_poll_pollen / poll_plant_prop_visits_OTU$total_plant_pollen

# Summarising at the Genus level
total_plant_visits <- poll_visit_pollen_data %>%   group_by(plant_sci_name) %>%  summarise(total_plant_visits= sum(visits, na.rm=TRUE), total_plant_pollen= sum(mean_pollen_load, na.rm=TRUE))
plant_poll_visits_genus <- poll_visit_pollen_data %>%   group_by(insect_genus, plant_sci_name) %>%  summarise(plant_poll_visits= sum(visits, na.rm=TRUE), plant_poll_pollen= sum(mean_pollen_load, na.rm=TRUE))
poll_plant_prop_visits_genus <- merge(x=total_plant_visits, y=plant_poll_visits_genus,by="plant_sci_name",all.x=TRUE)
poll_plant_prop_visits_genus$prop_visits <- poll_plant_prop_visits_genus$plant_poll_visits / poll_plant_prop_visits_genus$total_plant_visits
poll_plant_prop_visits_genus$prop_pollen <- poll_plant_prop_visits_genus$plant_poll_pollen / poll_plant_prop_visits_genus$total_plant_pollen

#Summarising at the broad pollinator guild level
total_plant_visits <- poll_visit_pollen_data %>%   group_by(plant_sci_name) %>%  summarise(total_plant_visits= sum(visits, na.rm=TRUE), total_plant_pollen= sum(mean_pollen_load, na.rm=TRUE))
plant_poll_visits_taxa <- poll_visit_pollen_data %>%   group_by(pollinator_taxa, plant_sci_name) %>%  summarise(plant_poll_visits= sum(visits, na.rm=TRUE), plant_poll_pollen= sum(mean_pollen_load, na.rm=TRUE))
poll_plant_prop_visits_taxa <- merge(x=total_plant_visits, y=plant_poll_visits_taxa,by="plant_sci_name",all.x=TRUE)
poll_plant_prop_visits_taxa$prop_visits <- poll_plant_prop_visits_taxa$plant_poll_visits / poll_plant_prop_visits_taxa$total_plant_visits
poll_plant_prop_visits_taxa$prop_pollen <- poll_plant_prop_visits_taxa$plant_poll_pollen / poll_plant_prop_visits_taxa$total_plant_pollen


########################################################################################################
##################                Step 2: Calculating importance of each crop         ##################
########################################################################################################

#Rename the plant name column
colnames(crop_importance_data)[2]  <- "plant_sci_name"

##Subset to relevant data - remove unnecessary columns
crop_importance_subset <- select(crop_importance_data, plant_sci_name, plant_genus,	plant_barcode, eng_name, poll_dependence,	final_HH_income, 	mean_rank_nutr_importance)

##Normalise data for each metric (divide all values by highest value)
crop_importance_subset$econ_importance <- crop_importance_subset$final_HH_income / (max(crop_importance_subset$final_HH_income))

crop_importance_subset$nutr_importance <- crop_importance_subset$mean_rank_nutr_importance / (max(crop_importance_subset$mean_rank_nutr_importance))


##Weight each crop's importance metric by its pollinator dependence so that more pollinator-dependent crops are given a higher weighting than less PD ones
crop_importance_subset$econ_importance_PD <- crop_importance_subset$econ_importance * crop_importance_subset$poll_dependence

crop_importance_subset$nutr_importance_PD <- crop_importance_subset$nutr_importance * crop_importance_subset$poll_dependence


##############################################################################################################################
###############      Step 3: Calculating multi-crop importance of each pollinator (data from all villages merged)     ########
##############################################################################################################################

##Merge pollinator importance data with crop importance data at each of the three taxonomic levels
pollinator_importance_merged_OTU <-merge(x=crop_importance_subset, y=poll_plant_prop_visits_OTU,by="plant_sci_name",all.x=TRUE)
pollinator_importance_merged_genus <-merge(x=crop_importance_subset, y=poll_plant_prop_visits_genus,by="plant_sci_name",all.x=TRUE)
pollinator_importance_merged_taxa <-merge(x=crop_importance_subset, y=poll_plant_prop_visits_taxa,by="plant_sci_name",all.x=TRUE)

##Multiply importance value by proportion visits and proportion pollen
# For pollinator_importance_merged_OTU
pollinator_importance_merged_OTU <- pollinator_importance_merged_OTU %>%
  mutate(poll_econ_importance = econ_importance_PD * prop_visits,
    poll_nutr_importance = nutr_importance_PD * prop_visits,
    poll_econ_importance_pollen = econ_importance_PD * prop_pollen,
    poll_nutr_importance_pollen = nutr_importance_PD * prop_pollen)

# For pollinator_importance_merged_genus
pollinator_importance_merged_genus <- pollinator_importance_merged_genus %>%
  mutate(poll_econ_importance = econ_importance_PD * prop_visits,
    poll_nutr_importance = nutr_importance_PD * prop_visits,
    poll_econ_importance_pollen = econ_importance_PD * prop_pollen,
    poll_nutr_importance_pollen = nutr_importance_PD * prop_pollen)

# For pollinator_importance_merged_taxa
pollinator_importance_merged_taxa <- pollinator_importance_merged_taxa %>%
  mutate(poll_econ_importance = econ_importance_PD * prop_visits,
    poll_nutr_importance = nutr_importance_PD * prop_visits,
    poll_econ_importance_pollen = econ_importance_PD * prop_pollen,
    poll_nutr_importance_pollen = nutr_importance_PD * prop_pollen)

#Remove all entries with no visits
pollinator_importance_merged_OTU <- pollinator_importance_merged_OTU[!(pollinator_importance_merged_OTU$prop_visits == ""), ]
pollinator_importance_merged_genus <- pollinator_importance_merged_genus[!(pollinator_importance_merged_genus$prop_visits == ""), ]
pollinator_importance_merged_taxa <- pollinator_importance_merged_taxa[!(pollinator_importance_merged_taxa$prop_visits == ""), ]

#Add a new column which sums the total importance value for each pollinator
importance_sum_OTU <- pollinator_importance_merged_OTU %>%   group_by(insect_OTU) %>%  summarise(sum_econ_value_pollen= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                       sum_nutr_value_pollen= sum(poll_nutr_importance_pollen, na.rm=TRUE)) %>%  ungroup()
#Merge this in to the original OTU dataset
pollinator_importance_merged_OTU <- merge(pollinator_importance_merged_OTU, importance_sum_OTU, by = "insect_OTU", all.x = TRUE)

#----Summarise data at the pollinator level#
pollinator_importance_sp <- pollinator_importance_merged_OTU %>%   group_by(insect_OTU) %>%  summarise(sum_econ_value= sum(poll_econ_importance, na.rm=TRUE),
                                                                                                sum_nutr_value= sum(poll_nutr_importance, na.rm=TRUE),
                                                                                                sum_econ_value_pollen= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                sum_nutr_value_pollen= sum(poll_nutr_importance_pollen, na.rm=TRUE))

pollinator_importance_genus <- pollinator_importance_merged_genus %>%   group_by(insect_genus) %>%  summarise(sum_econ_value= sum(poll_econ_importance, na.rm=TRUE),
                                                                                                              sum_nutr_value= sum(poll_nutr_importance, na.rm=TRUE),
                                                                                                              sum_econ_value_pollen= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                              sum_nutr_value_pollen= sum(poll_nutr_importance_pollen, na.rm=TRUE))

pollinator_importance_taxa <- pollinator_importance_merged_taxa %>%   group_by(pollinator_taxa) %>%  summarise(sum_econ_value= sum(poll_econ_importance, na.rm=TRUE),
                                                                                                               sum_nutr_value= sum(poll_nutr_importance, na.rm=TRUE),
                                                                                                               sum_econ_value_pollen= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                               sum_nutr_value_pollen= sum(poll_nutr_importance_pollen, na.rm=TRUE))


#Remove all entries with no data
pollinator_importance_sp <- pollinator_importance_sp[!(pollinator_importance_sp$insect_OTU == ""), ]
pollinator_importance_genus <- pollinator_importance_genus[!(pollinator_importance_genus$insect_genus == ""), ]
pollinator_importance_taxa <- pollinator_importance_taxa[!(pollinator_importance_taxa$pollinator_taxa == ""), ]
pollinator_importance_genus <- pollinator_importance_genus %>% drop_na(insect_genus)


##Summarise data for each metric (divide all values by total value so that they are expressed as a proportion)
#Visit data - economic
pollinator_importance_sp$econ_importance_prop <- pollinator_importance_sp$sum_econ_value / (sum(pollinator_importance_sp$sum_econ_value))
pollinator_importance_genus$econ_importance_prop <- pollinator_importance_genus$sum_econ_value / (sum(pollinator_importance_genus$sum_econ_value))
pollinator_importance_taxa$econ_importance_prop <- pollinator_importance_taxa$sum_econ_value / (sum(pollinator_importance_taxa$sum_econ_value))
#Pollen data - economic
pollinator_importance_sp$econ_importance_pollen_prop <- pollinator_importance_sp$sum_econ_value_pollen / (sum(pollinator_importance_sp$sum_econ_value_pollen))
pollinator_importance_genus$econ_importance_pollen_prop <- pollinator_importance_genus$sum_econ_value_pollen / (sum(pollinator_importance_genus$sum_econ_value_pollen))
pollinator_importance_taxa$econ_importance_pollen_prop <- pollinator_importance_taxa$sum_econ_value_pollen / (sum(pollinator_importance_taxa$sum_econ_value_pollen))
#Visit data - nutrition
pollinator_importance_sp$nutr_importance_prop <- pollinator_importance_sp$sum_nutr_value / (sum(pollinator_importance_sp$sum_nutr_value))
pollinator_importance_genus$nutr_importance_prop <- pollinator_importance_genus$sum_nutr_value / (sum(pollinator_importance_genus$sum_nutr_value))
pollinator_importance_taxa$nutr_importance_prop <- pollinator_importance_taxa$sum_nutr_value / (sum(pollinator_importance_taxa$sum_nutr_value))
#Pollen data - nutrition
pollinator_importance_sp$nutr_importance_pollen_prop <- pollinator_importance_sp$sum_nutr_value_pollen / (sum(pollinator_importance_sp$sum_nutr_value_pollen))
pollinator_importance_genus$nutr_importance_pollen_prop <- pollinator_importance_genus$sum_nutr_value_pollen / (sum(pollinator_importance_genus$sum_nutr_value_pollen))
pollinator_importance_taxa$nutr_importance_pollen_prop <- pollinator_importance_taxa$sum_nutr_value_pollen / (sum(pollinator_importance_taxa$sum_nutr_value_pollen))

##Calculate overall importance value - econ importance + nutr importance
#Visit-based
pollinator_importance_sp$overall_importance <- pollinator_importance_sp$nutr_importance_prop + pollinator_importance_sp$econ_importance_prop
pollinator_importance_genus$overall_importance <- pollinator_importance_genus$nutr_importance_prop + pollinator_importance_genus$econ_importance_prop
pollinator_importance_taxa$overall_importance <- pollinator_importance_taxa$nutr_importance_prop + pollinator_importance_taxa$econ_importance_prop
#Pollen-based
pollinator_importance_sp$overall_importance_pollen <- pollinator_importance_sp$nutr_importance_pollen_prop + pollinator_importance_sp$econ_importance_pollen_prop
pollinator_importance_genus$overall_importance_pollen <- pollinator_importance_genus$nutr_importance_pollen_prop + pollinator_importance_genus$econ_importance_pollen_prop
pollinator_importance_taxa$overall_importance_pollen <- pollinator_importance_taxa$nutr_importance_pollen_prop + pollinator_importance_taxa$econ_importance_pollen_prop


###############################################################################
#######      Step 4:  - Identifying best plants under each scenario    ########
###############################################################################

#Subset visitation data to leave only plant info and insect OTU
poll_visitation_plant_subset <- select(poll_visit_pollen_data, specimen_code,	pollinator_taxa,	insect_order,	insect_family,	insect_genus,	insect_species,	insect_OTU,	plant_family,	plant_sci_name,	plant_eng_name,	plant_category,	plant_barcode, visits)

##Merge visitation data with pollinator importance data based on insect OTU
plant_importance_merged<-merge(x=pollinator_importance_sp, y=poll_visitation_plant_subset, by="insect_OTU",all.x=TRUE)

#Remove all entries with no plant species names
plant_importance_merged <- plant_importance_merged[!(plant_importance_merged$plant_sci_name == ""), ]

#----Summarise data at the plant level#
plant_importance <- plant_importance_merged %>%   group_by(plant_family,plant_sci_name,plant_eng_name,plant_category) %>%  summarise(sum_visits= sum(visits, na.rm=TRUE),
                                                                                                                                     sum_econ_value= sum(econ_importance_prop, na.rm=TRUE),
                                                                                                                                     sum_nutr_value= sum(nutr_importance_prop, na.rm=TRUE),
                                                                                                                                     sum_econ_value_pollen= sum(econ_importance_pollen_prop, na.rm=TRUE),
                                                                                                                                   sum_nutr_value_pollen= sum(nutr_importance_pollen_prop, na.rm=TRUE))

##Summarise data for each metric by dividing each importance value by the sum of all importance values so that it is expressed as a proportion of total importance
#Visit data - economic
plant_importance$econ_importance_prop <- plant_importance$sum_econ_value / (sum(plant_importance$sum_econ_value))
#Pollen data - economic
plant_importance$econ_importance_pollen_prop <- plant_importance$sum_econ_value_pollen / (sum(plant_importance$sum_econ_value_pollen))
#Visit data - nutrition
plant_importance$nutr_importance_prop <- plant_importance$sum_nutr_value / (sum(plant_importance$sum_nutr_value))
#Pollen data - nutrition
plant_importance$nutr_importance_pollen_prop <- plant_importance$sum_nutr_value_pollen / (sum(plant_importance$sum_nutr_value_pollen))

#----Remove crop plants from plant importance scores
plant_importance_wild <- plant_importance[!(plant_importance$plant_category == "crop"), ]


#################################################################
#######      Printing out results of importance analysis    ########
#################################################################

#----Export importance summaries#
wb <- createWorkbook()
addWorksheet(wb, sheet = "Info", gridLines = TRUE)
addWorksheet(wb, sheet = "Raw", gridLines = TRUE)
addWorksheet(wb, sheet = "pollinator_importance_sp", gridLines = TRUE)
addWorksheet(wb, sheet = "pollinator_importance_genus", gridLines = TRUE)
addWorksheet(wb, sheet = "pollinator_importance_taxa", gridLines = TRUE)
addWorksheet(wb, sheet = "plant_importance", gridLines = TRUE)
writeData(wb, sheet = "Raw", x = pollinator_importance_merged_OTU, rowNames = FALSE)
writeData(wb, sheet = "pollinator_importance_sp", x = pollinator_importance_sp, rowNames = FALSE)
writeData(wb, sheet = "pollinator_importance_genus", x = pollinator_importance_genus, rowNames = FALSE)
writeData(wb, sheet = "pollinator_importance_taxa", x = pollinator_importance_taxa, rowNames = FALSE)
writeData(wb, sheet = "plant_importance", x = plant_importance_wild, rowNames = FALSE)
freezePane(wb, sheet = "Info", firstRow = TRUE)
freezePane(wb, sheet = "Raw", firstRow = TRUE)
freezePane(wb, sheet = "pollinator_importance_sp", firstRow = TRUE)
freezePane(wb, sheet = "pollinator_importance_genus", firstRow = TRUE)
freezePane(wb, sheet = "pollinator_importance_taxa", firstRow = TRUE)
freezePane(wb, sheet = "plant_importance", firstRow = TRUE)
saveWorkbook(wb, file.path(output.path, paste0( "MP_poll-importance_analysis_", Sys.Date(), ".xlsx")), overwrite = TRUE)


#################################################################
#######      Figure 1a: Crop importance  heat map       ########
#################################################################

#Restructure data to create two convert it to long format, with all importance values shown in a single column and the variable (econ/nutr) specified in a separate column
crop_econ_importance <- select(crop_importance_subset, eng_name, econ_importance)
crop_nutr_importance <- select(crop_importance_subset, eng_name, nutr_importance)
crop_econ_importance$priority <- "Economic"
crop_nutr_importance$priority <- "Nutrition"
colnames(crop_econ_importance)[2]  <- "importance_value"
colnames(crop_nutr_importance)[2]  <- "importance_value"
crop_importance <- rbind(crop_econ_importance, crop_nutr_importance)

#Reorder values from most to least important
crop_importance <- crop_importance %>%
  mutate(eng_name = fct_reorder(eng_name, importance_value))

crop_importance$priority <- factor(crop_importance$priority, levels = c("Nutrition", "Economic"))

# Plot heatmap
crop_heatmap<- ggplot(crop_importance, aes(priority, eng_name, fill=importance_value)) +
  geom_tile(color="black", lwd=0.5, linetype=1) + scale_fill_gradientn(colors=hcl.colors(10, "YlOrRd"), trans='reverse') + theme_bw() +
  guides(fill=guide_colourbar(title="", barwidth=1, barheight=10)) +
  theme(axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11, face="bold"),
        axis.title=element_text(size=11, face="bold"),
        plot.title = element_text(color="black", size=14, face="bold", hjust=0.5)) +
  theme(strip.text.x = element_text(size = 11, face = "bold")) + theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("")+
  scale_x_discrete(labels=c("Nutritional\nValue", "Economic\nValue"))

ggsave(plot=crop_heatmap, filename= file.path(output.path,"crop_importance_heatmap.svg"), width=3, height=6, dpi=500)
ggsave(plot=crop_heatmap, filename= file.path(output.path,"crop_importance_heatmap.png"), width=3, height=6, dpi=500)

#################################################################
#######      Figure 1b - Crop importance  correlation    ########
#################################################################

# Testing for a correlation between the economic and nutritional importance scores for crops, pollinators and plants
lm_crop_importance_correlation <-  lm(log10(nutr_importance) ~ log10(econ_importance), data = crop_importance_subset)
lm_poll_importance_correlation <-  lm(nutr_importance_pollen_prop ~ econ_importance_pollen_prop, data = pollinator_importance_sp)
lm_plant_importance_correlation <-  lm(sum_nutr_value_pollen ~ sum_econ_value_pollen, data = plant_importance_wild)

summary(lm_crop_importance_correlation)
summary(lm_poll_importance_correlation)
summary(lm_plant_importance_correlation)

## Plotting the relationship between the economic and nutritional value of each crop (Figure 1b)
crop_importance_correlation <- ggplot(crop_importance_subset, aes(x = nutr_importance, y = econ_importance)) + 
  geom_point(color='black', size = 3, shape = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +  # Add 1:1 line to show equal economic & nutritional value
  geom_text(x = 0.5, y = 0.55, label = expression(italic("Economic value = Nutritional value")), 
            color = "black", size = 4, angle = 45) +  # Add text label
  theme_classic(base_size = 22) +
  labs(x = bquote("Nutritional Value"),
       y = bquote("Economic Value")) +
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +
    geom_text_repel(aes(label = eng_name),
                  box.padding = 0.3, 
                  point.padding = 0.5,
                  segment.color = 'black',
                  segment.size = 0.2)+
  ylim(0,1)+   xlim(0,1)

ggsave(plot=crop_importance_correlation, filename= file.path(output.path, "Crop_importance_correlation.png"), width=6, height=6, dpi=500)
ggsave(plot=crop_importance_correlation, filename= file.path(output.path, "Crop_importance_correlation.svg"), width=6, height=6, dpi=500)

#########################################################################################
#######      Figure 2a - Plotting importance of each pollinator guild by crop    ########
#########################################################################################

crop_poll_visit_prop <- pollinator_importance_merged_taxa

#Remove irrelevant taxa from pollinator importance data (bugs,moths & unknown) and crops with <5 visits #
crop_poll_visit_prop <- crop_poll_visit_prop[!(crop_poll_visit_prop$pollinator_taxa == "Bug" | crop_poll_visit_prop$pollinator_taxa == "Unknown" | crop_poll_visit_prop$pollinator_taxa == "Moth" | crop_poll_visit_prop$eng_name == "Tree tomato"  | crop_poll_visit_prop$eng_name == "Aubergine" | crop_poll_visit_prop$eng_name == "Green bean" | crop_poll_visit_prop$eng_name == "Plum"), ]

#Specify order of variables
crop_poll_visit_prop$pollinator_taxa <- factor(crop_poll_visit_prop$pollinator_taxa, levels = c("Fly","Wasp",	"Beetle",	"Sawfly","Butterfly", "Solitary bee",	"Bumblebee","Apis laboriosa", "Apis cerana"))
crop_poll_visit_prop$eng_name <- factor(crop_poll_visit_prop$eng_name, levels = c("Jumli bean","Sunflower","Buckwheat", "Chilli", "Pumpkin","Scarlet bean", "Mustard seed","Tomato","Peach","Slipper gourd","Cucumber",	"Apple"))

col <- c("#8098A2", "#CDD4DC","#F5DC83", "#BD338F","#EB8252","#8FA33F",  "#5F7929", "#616020", "#014820")

#Plotting based on visitation data alone (note that this plot is not shown in the MS)
crop_taxon_visits_plot <- ggplot(crop_poll_visit_prop, aes(fill=pollinator_taxa, y=eng_name, x=prop_visits)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = col, name= "Pollinator Taxa") +
  xlab("Proportion visits") + ylab("Crop")+
  theme_bw() +
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14)) +
  theme(plot.title = element_text(color="black", size=15, face="bold", hjust=0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("Visitation Data")

#Plotting based on pollen-weighted visitation data (i.e. visitation freq*pollen carrying capacity) - this is Figure 2a
crop_taxon_pollen_plot <- ggplot(crop_poll_visit_prop, aes(fill=pollinator_taxa, y=eng_name, x=prop_pollen)) +
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values = col, name= "Pollinator Taxa") +
  xlab("Proportion pollen transport") + ylab("")+
  theme_bw() +
  theme(legend.position = "none")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14)) +
  theme(plot.title = element_text(color="black", size=15, face="bold", hjust=0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

crop_taxon_importance_grouped <- ggarrange(crop_taxon_pollen_plot,  nrow = 1, labels= "", common.legend = TRUE, legend="right")


ggsave(plot=crop_taxon_importance_grouped, filename= file.path(output.path, "Crop_poll_taxa_plot.svg"), width=6, height=5, dpi=500)
ggsave(plot=crop_taxon_importance_grouped, filename= file.path(output.path, "Crop_poll_taxa_plot.png"), width=6, height=5, dpi=500)

#######################################################################################################
#######   Figure 2b: Calculating & plotting Shannon diversity of pollinators to each crop       ########
#######################################################################################################

#Removing NA values in dataset
pollinator_importance_clean <- pollinator_importance_merged_OTU %>% replace(is.na(.), 0)

#Remove crop plants with only one visitation event
pollinator_importance_clean <- pollinator_importance_clean[!(pollinator_importance_clean$plant_sci_name == "Cyphomandra betacea"), ]

crop_poll_visits <- dcast(pollinator_importance_clean, eng_name ~ insect_OTU, value.var="prop_visits")
crop_poll_pollen <- dcast(pollinator_importance_clean, eng_name ~ insect_OTU, value.var="prop_pollen")

crop_poll_visits <- crop_poll_visits %>%   column_to_rownames("eng_name")
crop_poll_pollen <- crop_poll_pollen %>%   column_to_rownames("eng_name")

crop_poll_visits <- crop_poll_visits %>% replace(is.na(.), 0)
crop_poll_pollen <- crop_poll_pollen %>% replace(is.na(.), 0)

#Calculating Shannon diversity index for each crop and each metric (visitation and pollen)
crop_poll_div_visits <- diversity(crop_poll_visits)
crop_poll_div_pollen <- diversity(crop_poll_pollen)

#Restructuring results to allow plotting and analysis
crop_poll_div_visits <- data.frame(crop_poll_div_visits)
crop_poll_div_visits$crop_species <- rownames(crop_poll_div_visits)
crop_poll_div_visits$data_source <- "visitation"
rownames(crop_poll_div_visits) <- NULL
colnames(crop_poll_div_visits)[1] = "diversity_score"
crop_poll_div_pollen <- data.frame(crop_poll_div_pollen)
crop_poll_div_pollen$crop_species <- rownames(crop_poll_div_pollen)
crop_poll_div_pollen$data_source <- "pollen"
rownames(crop_poll_div_pollen) <- NULL
colnames(crop_poll_div_pollen)[1] = "diversity_score"

## Re-order crops based on their pollinator community diversity
visit_diversity_df <- crop_poll_div_visits %>%  mutate(crop_species = fct_reorder(crop_species, diversity_score))
pollen_diversity_df <- crop_poll_div_pollen %>%  mutate(crop_species = fct_reorder(crop_species, diversity_score))

# Plot diversity of pollinators of each crop based on visitation data alone (note that this plot is not used in the MS but is almost identical to figure 2b)
visit_diversity <- ggplot(visit_diversity_df) +  geom_col(aes(diversity_score, crop_species), fill = "#8fa33fff", width = 0.6)+
  theme_classic() + ylab("Crop Plant") + xlab("Diversity of Insect Visitors") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18),
        plot.title = element_text(color="black", size=18, face="bold", hjust=0.5))+
  ggtitle("Visitation Data")


#Plotting diversity based on pollen-weighted visitation data (i.e. visitation freq*pollen carrying capacity) - this is Figure 2b
pollen_diversity <- ggplot(pollen_diversity_df) +  geom_col(aes(diversity_score, crop_species), fill = "#8fa33fff", width = 0.6)+
  theme_classic() + ylab("") + xlab("Diversity of Pollen Transport") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18),
        plot.title = element_text(color="black", size=18, face="bold", hjust=0.5))

crop_poll_diversity_plot <- plot_grid(pollen_diversity,  nrow = 1, labels="")

ggsave(plot=crop_poll_diversity_plot, filename= file.path(output.path, "Crop_pollinator_diversity.svg"), width=6, height=7, dpi=500)
ggsave(plot=crop_poll_diversity_plot, filename= file.path(output.path, "Crop_pollinator_diversity.png"), width=6, height=7, dpi=500)


#################################################################################################
#######        Figure 3:   Plotting crop-pollinator network as alluvial diagrams       ########
#################################################################################################

#Subset data to only show the pollinator taxa scoring above 0.01 on each scale (this removes the taxa with very few visits)
alluvial_data <- pollinator_importance_merged_OTU
alluvial_data$insect_OTU <- gsub("_", " ", alluvial_data$insect_OTU)
alluvial_data_econ <- subset(alluvial_data, sum_econ_value_pollen >= 0.01)
alluvial_data_nutr <- subset(alluvial_data, sum_nutr_value_pollen >= 0.01)

# This plots the crop-pollinator interactions with crops weighted by their NUTRITIONAL VALUE
nutr_web <- ggplot(data = alluvial_data_nutr,
                   aes(axis1 = reorder(eng_name, -nutr_importance_PD), 
                       axis2 = reorder(insect_OTU, -sum_nutr_value_pollen), 
                       y = poll_nutr_importance_pollen)) +
  geom_alluvium(aes(fill = insect_OTU), curve_type = "cubic") +
  geom_stratum(fill = "white")+
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_text(x = 1, y = 4.3, label = "Crop Importance", size = 4,fontface = "bold", color = "black")+  
  geom_text(x = 2, y = 4.3, label = "Pollinator Importance", size = 4, fontface = "bold", color = "black")+  
  scale_x_discrete(limits = c("Crop", "Pollinator"), expand = c(0.2, 0.1)) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  scale_fill_viridis_d()+
  theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10)) +  # Increase top margin for the title
  labs(title = "Food production")

# This plots the crop-pollinator interactions with crops weighted by their ECONOMIC VALUE
econ_web <- ggplot(data = alluvial_data_econ,
                   aes(axis1 = reorder(eng_name, -econ_importance_PD), 
                       axis2 = reorder(insect_OTU, -sum_econ_value_pollen), 
                       y = poll_econ_importance_pollen)) +
  geom_alluvium(aes(fill = insect_OTU), curve_type = "cubic") +
  geom_stratum(fill = "white")+
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_text(x = 1, y = 1.64, label = "Crop Importance", size = 4,fontface = "bold", color = "black")+  
  geom_text(x = 2, y = 1.64, label = "Pollinator Importance", size = 4, fontface = "bold", color = "black")+  
  scale_x_discrete(limits = c("Crop", "Pollinator"), expand = c(0.2, 0.1)) +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  scale_fill_viridis_d()+
  theme(plot.margin = margin(t = 15, r = 10, b = 10, l = 10)) +  # Increase top margin for the title
  labs(title = "Income generation")

econ_nutr_webs_grouped <- ggarrange(nutr_web,econ_web,  nrow = 1, labels= "", common.legend = FALSE)

ggsave(plot=econ_nutr_webs_grouped, filename= file.path(output.path, "Crop_poll_web_plot.svg"), width=13, height=7, dpi=500, bg="white")
ggsave(plot=econ_nutr_webs_grouped, filename= file.path(output.path, "Crop_poll_web_plot.png"), width=13, height=7, dpi=500, bg="white")


##########################################################################################
#######     Figures 4 & 5: Pollinator importance and diversity values by village    #########
##########################################################################################

##Note these analyses are done on a village-by-village basis - see bottom of script for this village-level analyis


########################################################################################################
#######    Figure 6: Simulating agricultural specialisation and analysing network structure ########
########################################################################################################

### Preparing the raw dataset
full_web <- pollinator_importance_merged_OTU
full_web$insect_OTU <- gsub("_", " ", full_web$insect_OTU)

# Summarise crop-pollinator interactions
web_summary <- full_web %>%
  group_by(eng_name, insect_OTU) %>%
  summarise(pollen_transport = sum(poll_econ_importance_pollen)) %>%
  spread(key = insect_OTU, value = pollen_transport, fill = 0)

#Make crop name the row name
web_summary <- column_to_rownames(web_summary, var = "eng_name")

######################################################################
### Remove different proportions of crops to show specialisation
####################################################################

#Calculate row sums 
row_sums <- rowSums(web_summary)
#Identify the indices of the rows with the lowest sums - i.e the crops with the lowest economic values
lowest_01 <- order(row_sums)[1:1]
lowest_02 <- order(row_sums)[1:2]
lowest_03 <- order(row_sums)[1:3]
lowest_04 <- order(row_sums)[1:4]
lowest_05 <- order(row_sums)[1:5]
lowest_06 <- order(row_sums)[1:6]
lowest_07 <- order(row_sums)[1:7]
lowest_08 <- order(row_sums)[1:8]
lowest_09 <- order(row_sums)[1:9]
lowest_10 <- order(row_sums)[1:10]
lowest_11 <- order(row_sums)[1:11]
lowest_12 <- order(row_sums)[1:12]
lowest_13 <- order(row_sums)[1:13]
lowest_14 <- order(row_sums)[1:14]

# Remove those rows from the table  - i.e remove crops in order of low to high economic value
specialisation_01 <- web_summary[-lowest_01, ]
specialisation_02 <- web_summary[-lowest_02, ]
specialisation_03 <- web_summary[-lowest_03, ]
specialisation_04 <- web_summary[-lowest_04, ]
specialisation_05 <- web_summary[-lowest_05, ]
specialisation_06 <- web_summary[-lowest_06, ]
specialisation_07 <- web_summary[-lowest_07, ]
specialisation_08 <- web_summary[-lowest_08, ]
specialisation_09 <- web_summary[-lowest_09, ]
specialisation_10 <- web_summary[-lowest_10, ]
specialisation_11 <- web_summary[-lowest_11, ]
specialisation_12 <- web_summary[-lowest_12, ]
specialisation_13 <- web_summary[-lowest_13, ]
specialisation_14 <- web_summary[-lowest_14, ]

##Convert to all dataframes to matrices
matrix_full <- as.matrix(web_summary)
matrix_01 <- as.matrix(specialisation_01)
matrix_02 <- as.matrix(specialisation_02)
matrix_03 <- as.matrix(specialisation_03)
matrix_04 <- as.matrix(specialisation_04)
matrix_05 <- as.matrix(specialisation_05)
matrix_06 <- as.matrix(specialisation_06)
matrix_07 <- as.matrix(specialisation_07)
matrix_08 <- as.matrix(specialisation_08)
matrix_09 <- as.matrix(specialisation_09)
matrix_10 <- as.matrix(specialisation_10)
matrix_11 <- as.matrix(specialisation_11)
matrix_12 <- as.matrix(specialisation_12)
matrix_13 <- as.matrix(specialisation_13)
matrix_14 <- as.matrix(specialisation_14)

# Caclulate Network stats on each matrix
diversity_datasets <- list(diversity_full <- networklevel(matrix_full, index = "Shannon diversity",level = "higher"),
networklevel(matrix_01, index = "Shannon diversity",level = "higher"),
networklevel(matrix_02, index = "Shannon diversity",level = "higher"),
networklevel(matrix_03, index = "Shannon diversity",level = "higher"),
networklevel(matrix_04, index = "Shannon diversity",level = "higher"),
networklevel(matrix_05, index = "Shannon diversity",level = "higher"),
networklevel(matrix_06, index = "Shannon diversity",level = "higher"),
networklevel(matrix_07, index = "Shannon diversity",level = "higher"),
networklevel(matrix_08, index = "Shannon diversity",level = "higher"),
networklevel(matrix_09, index = "Shannon diversity",level = "higher"),
networklevel(matrix_10, index = "Shannon diversity",level = "higher"),
networklevel(matrix_11, index = "Shannon diversity",level = "higher"),
networklevel(matrix_12, index = "Shannon diversity",level = "higher"),
networklevel(matrix_13, index = "Shannon diversity",level = "higher"),
networklevel(matrix_14, index = "Shannon diversity",level = "higher"))

for (i in 1:length(diversity_datasets)) {
  diversity_datasets[[i]]$value <- i
}

# Combine datasets into a single data frame
combined_data <- do.call(rbind, diversity_datasets)
combined_data <- as.data.frame(combined_data)
combined_data$value <- as.numeric(combined_data$value)
combined_data$`Shannon diversity` <- as.numeric(combined_data$`Shannon diversity`)

#Calculate percentage specialisation by dividing the number of crops removed by the total number of crops
combined_data$specialisation <- ((combined_data$value)/16) * 100 

## Plot diversity against percentage specialisation
diversity_plot <- ggplot(combined_data, aes(x = specialisation, y = `Shannon diversity`)) +
  geom_line() +
  geom_point() +
  labs(x = "Percentage specialisation", y = "Shannon diversity of pollinator community") +
  ggtitle("")+
  ylim(0,4)+
  theme_minimal() +  # You can change this to another theme if needed
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

ggsave(plot=diversity_plot, filename= file.path(output.path, "diversity_plot.svg"), width=8, height=4, dpi=500)
ggsave(plot=diversity_plot, filename= file.path(output.path, "diversity_plot.png"), width=8, height=4, dpi=500)


###################################################################################################################
##### Plot the plant-pollinator networks for different stages of agricultural specialisation (Part of Fig. 6) #####
#################################################################################################################

#Full web
sort_matrix_full <- sortweb(matrix_full, sort.order="dec", sequence=NULL)

svg(file.path(output.path,"webplot_01.svg"), width=8,height=6)
webplot_01 <- plotweb(as.data.frame(sort_matrix_full),
                 method = "normal", empty = TRUE, labsize = 1, ybig = 1.2, y.width.low = 0.05, 
                 y.width.high = 0.05, col.interaction="grey80", col.high = "orange",text.rot=90, 
                 col.low="darkolivegreen", low.lab.dis = NULL, abuns.type="additional", high.lablength = 0)
dev.off()

#Three crops removed
sort_matrix_03 <- sortweb(matrix_03, sort.order="dec", sequence=NULL)
svg(file.path(output.path,"webplot_03.svg"), width=8,height=6)
webplot_03 <- plotweb(as.data.frame(sort_matrix_03),
                      method = "normal", empty = TRUE, labsize = 1, ybig = 1.2, y.width.low = 0.05, 
                      y.width.high = 0.05, col.interaction="grey80", col.high = "orange",text.rot=90, 
                      col.low="darkolivegreen", low.lab.dis = NULL, abuns.type="additional", high.lablength = 0)
dev.off()

#Six crops removed
sort_matrix_06 <- sortweb(matrix_06, sort.order="dec", sequence=NULL)
svg(file.path(output.path,"webplot_06.svg"), width=8,height=6)
webplot_06 <- plotweb(as.data.frame(sort_matrix_06),
                      method = "normal", empty = TRUE, labsize = 1, ybig = 1.2, y.width.low = 0.05, 
                      y.width.high = 0.05, col.interaction="grey80", col.high = "orange",text.rot=90, 
                      col.low="darkolivegreen", low.lab.dis = NULL, abuns.type="additional", high.lablength = 0)
dev.off()

#Nine crops removed
sort_matrix_09 <- sortweb(matrix_09, sort.order="dec", sequence=NULL)
svg(file.path(output.path,"webplot_09.svg"), width=8,height=6)
webplot_09 <- plotweb(as.data.frame(sort_matrix_09),
                      method = "normal", empty = TRUE, labsize = 1, ybig = 1.2, y.width.low = 0.05, 
                      y.width.high = 0.05, col.interaction="grey80", col.high = "orange",text.rot=90, 
                      col.low="darkolivegreen", low.lab.dis = NULL, abuns.type="additional", high.lablength = 0)
dev.off()

#Twelve crops removed
sort_matrix_12 <- sortweb(matrix_12, sort.order="dec", sequence=NULL)
svg(file.path(output.path,"webplot_12.svg"), width=8,height=6)
webplot_12 <- plotweb(as.data.frame(sort_matrix_12),
                      method = "normal", empty = TRUE, labsize = 1, ybig = 1.2, y.width.low = 0.05, 
                      y.width.high = 0.05, col.interaction="grey80", col.high = "orange",text.rot=90, 
                      col.low="darkolivegreen", low.lab.dis = NULL, abuns.type="additional", high.lablength = 0)
dev.off()

#Fourteen crops removed
sort_matrix_14 <- sortweb(matrix_14, sort.order="dec", sequence=NULL)
svg(file.path(output.path,"webplot_14.svg"), width=8,height=6)
webplot_14 <- plotweb(as.data.frame(sort_matrix_14),
                      method = "normal", empty = TRUE, labsize = 1, ybig = 1.2, y.width.low = 0.05, 
                      y.width.high = 0.05, col.interaction="grey80", col.high = "orange",text.rot=90, 
                      col.low="darkolivegreen", low.lab.dis = NULL, abuns.type="additional", high.lablength = 0)
dev.off()


### The diversity plot and these different webplots are combined into a single figure in Inkscape to show how the structure of the network changes as farmers specialise


#############################################################################################################################################
#################################################     VILLAGE-LEVEL ANALYSES  ###############################################################
#############################################################################################################################################

#The following script repeats the calculations of pollinator and plant importance scores and diversity values in steps 3 & 4 above
# It makes the calculations at a village-level however

#Create copy of datasets to conduct village-level analysis on
poll_visit_pollen_data_village <- poll_visit_pollen_data


#################################################################################################################################
###########     Step 3B: Calculating village-level multi-crop importance of each pollinator (each village separately)     ########
#################################################################################################################################

##Create new column with pollinator domestication status
poll_visit_pollen_data_village <- poll_visit_pollen_data_village %>%
  mutate(wild_managed = ifelse(insect_OTU == 'Apis_cerana', 'Apis cerana', 'Wild insects'))

###Calculate proportion of visits made by each pollinator to each plant species - by village#
#OTU level
total_plant_visits_village <- poll_visit_pollen_data_village %>%   group_by(plant_sci_name, village_code) %>%  dplyr::summarise(total_plant_visits_village= sum(visits, na.rm=TRUE), total_plant_pollen= sum(mean_pollen_load, na.rm=TRUE))
plant_poll_visits_OTU_village <- poll_visit_pollen_data_village %>%   group_by(insect_OTU, wild_managed, plant_sci_name, village_code) %>%  dplyr::summarise(plant_poll_visits= sum(visits, na.rm=TRUE), plant_poll_pollen= sum(mean_pollen_load, na.rm=TRUE))
total_plant_visits_village$village_plant <- paste(total_plant_visits_village$village_code,"_",total_plant_visits_village$plant_sci_name)
plant_poll_visits_OTU_village$village_plant <- paste(plant_poll_visits_OTU_village$village_code,"_",plant_poll_visits_OTU_village$plant_sci_name)
plant_poll_visits_OTU_village <- subset(plant_poll_visits_OTU_village, select = -c(plant_sci_name, village_code))

poll_plant_prop_visits_OTU_village <- merge(x=total_plant_visits_village, y=plant_poll_visits_OTU_village,by="village_plant",all.x=TRUE)
poll_plant_prop_visits_OTU_village$prop_visits <- poll_plant_prop_visits_OTU_village$plant_poll_visits / poll_plant_prop_visits_OTU_village$total_plant_visits_village
poll_plant_prop_visits_OTU_village$prop_pollen <- poll_plant_prop_visits_OTU_village$plant_poll_pollen / poll_plant_prop_visits_OTU_village$total_plant_pollen

#Genus level
total_plant_visits_village <- poll_visit_pollen_data_village %>%   group_by(plant_sci_name, village_code) %>%  dplyr::summarise(total_plant_visits_village= sum(visits, na.rm=TRUE), total_plant_pollen= sum(mean_pollen_load, na.rm=TRUE))
plant_poll_visits_genus_village <- poll_visit_pollen_data_village %>%   group_by(insect_genus, wild_managed, plant_sci_name, village_code) %>%  dplyr::summarise(plant_poll_visits= sum(visits, na.rm=TRUE), plant_poll_pollen= sum(mean_pollen_load, na.rm=TRUE))
total_plant_visits_village$village_plant <- paste(total_plant_visits_village$village_code,"_",total_plant_visits_village$plant_sci_name)
plant_poll_visits_genus_village$village_plant <- paste(plant_poll_visits_genus_village$village_code,"_",plant_poll_visits_genus_village$plant_sci_name)
plant_poll_visits_genus_village <- subset(plant_poll_visits_genus_village, select = -c(plant_sci_name, village_code))

poll_plant_prop_visits_genus_village <- merge(x=total_plant_visits_village, y=plant_poll_visits_genus_village,by="village_plant",all.x=TRUE)
poll_plant_prop_visits_genus_village$prop_visits <- poll_plant_prop_visits_genus_village$plant_poll_visits / poll_plant_prop_visits_genus_village$total_plant_visits_village
poll_plant_prop_visits_genus_village$prop_pollen <- poll_plant_prop_visits_genus_village$plant_poll_pollen / poll_plant_prop_visits_genus_village$total_plant_pollen

#Guild level
total_plant_visits_village <- poll_visit_pollen_data_village %>%   group_by(plant_sci_name, village_code) %>%  dplyr::summarise(total_plant_visits_village= sum(visits, na.rm=TRUE), total_plant_pollen= sum(mean_pollen_load, na.rm=TRUE))
plant_poll_visits_taxa_village <- poll_visit_pollen_data_village %>%   group_by(pollinator_taxa, wild_managed, plant_sci_name, village_code) %>%  dplyr::summarise(plant_poll_visits= sum(visits, na.rm=TRUE), plant_poll_pollen= sum(mean_pollen_load, na.rm=TRUE))
total_plant_visits_village$village_plant <- paste(total_plant_visits_village$village_code,"_",total_plant_visits_village$plant_sci_name)
plant_poll_visits_taxa_village$village_plant <- paste(plant_poll_visits_taxa_village$village_code,"_",plant_poll_visits_taxa_village$plant_sci_name)
plant_poll_visits_taxa_village <- subset(plant_poll_visits_taxa_village, select = -c(plant_sci_name, village_code))

poll_plant_prop_visits_taxa_village <- merge(x=total_plant_visits_village, y=plant_poll_visits_taxa_village,by="village_plant",all.x=TRUE)
poll_plant_prop_visits_taxa_village$prop_visits <- poll_plant_prop_visits_taxa_village$plant_poll_visits / poll_plant_prop_visits_taxa_village$total_plant_visits_village
poll_plant_prop_visits_taxa_village$prop_pollen <- poll_plant_prop_visits_taxa_village$plant_poll_pollen / poll_plant_prop_visits_taxa_village$total_plant_pollen

##Merge pollinator data with crop importance data
visitation_importance_merged_OTU_village <-merge(x=crop_importance_subset, y=poll_plant_prop_visits_OTU_village,by="plant_sci_name",all.x=TRUE)
visitation_importance_merged_genus_village <-merge(x=crop_importance_subset, y=poll_plant_prop_visits_genus_village,by="plant_sci_name",all.x=TRUE)
visitation_importance_merged_taxa_village <-merge(x=crop_importance_subset, y=poll_plant_prop_visits_taxa_village,by="plant_sci_name",all.x=TRUE)

##Multiply importance value by proportion visits and proportion pollen
visitation_importance_merged_OTU_village$poll_econ_importance <- visitation_importance_merged_OTU_village$econ_importance_PD * visitation_importance_merged_OTU_village$prop_visits
visitation_importance_merged_OTU_village$poll_nutr_importance <- visitation_importance_merged_OTU_village$nutr_importance_PD * visitation_importance_merged_OTU_village$prop_visits
visitation_importance_merged_OTU_village$poll_econ_importance_pollen <- visitation_importance_merged_OTU_village$econ_importance_PD * visitation_importance_merged_OTU_village$prop_pollen
visitation_importance_merged_OTU_village$poll_nutr_importance_pollen <- visitation_importance_merged_OTU_village$nutr_importance_PD * visitation_importance_merged_OTU_village$prop_pollen

visitation_importance_merged_genus_village$poll_econ_importance <- visitation_importance_merged_genus_village$econ_importance_PD * visitation_importance_merged_genus_village$prop_visits
visitation_importance_merged_genus_village$poll_nutr_importance <- visitation_importance_merged_genus_village$nutr_importance_PD * visitation_importance_merged_genus_village$prop_visits
visitation_importance_merged_genus_village$poll_econ_importance_pollen <- visitation_importance_merged_genus_village$econ_importance_PD * visitation_importance_merged_genus_village$prop_pollen
visitation_importance_merged_genus_village$poll_nutr_importance_pollen <- visitation_importance_merged_genus_village$nutr_importance_PD * visitation_importance_merged_genus_village$prop_pollen

visitation_importance_merged_taxa_village$poll_econ_importance <- visitation_importance_merged_taxa_village$econ_importance_PD * visitation_importance_merged_taxa_village$prop_visits
visitation_importance_merged_taxa_village$poll_nutr_importance <- visitation_importance_merged_taxa_village$nutr_importance_PD * visitation_importance_merged_taxa_village$prop_visits
visitation_importance_merged_taxa_village$poll_econ_importance_pollen <- visitation_importance_merged_taxa_village$econ_importance_PD * visitation_importance_merged_taxa_village$prop_pollen
visitation_importance_merged_taxa_village$poll_nutr_importance_pollen <- visitation_importance_merged_taxa_village$nutr_importance_PD * visitation_importance_merged_taxa_village$prop_pollen

#Remove all entries with no visits
visitation_importance_merged_OTU_village <- visitation_importance_merged_OTU_village[!(visitation_importance_merged_OTU_village$prop_visits == ""), ]
visitation_importance_merged_genus_village <- visitation_importance_merged_genus_village[!(visitation_importance_merged_genus_village$prop_visits == ""), ]
visitation_importance_merged_taxa_village <- visitation_importance_merged_taxa_village[!(visitation_importance_merged_taxa_village$prop_visits == ""), ]
visitation_importance_merged_genus_village <- visitation_importance_merged_genus_village %>% drop_na(insect_genus)


#----Summarise data at the village level to work our proportion of visits/pollen contributed by each taxon
village_level_sp <- visitation_importance_merged_OTU_village %>%   group_by(village_code) %>%  dplyr::summarise(sum_econ_value_village= sum(poll_econ_importance, na.rm=TRUE),
                                                                                                 sum_nutr_value_village= sum(poll_nutr_importance, na.rm=TRUE),
                                                                                                 sum_econ_value_pollen_village= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                 sum_nutr_value_pollen_village= sum(poll_nutr_importance_pollen, na.rm=TRUE))

village_level_genus <- visitation_importance_merged_genus_village %>%   group_by(village_code) %>%  dplyr::summarise(sum_econ_value_village= sum(poll_econ_importance, na.rm=TRUE),
                                                                                                      sum_nutr_value_village= sum(poll_nutr_importance, na.rm=TRUE),
                                                                                                      sum_econ_value_pollen_village= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                      sum_nutr_value_pollen_village= sum(poll_nutr_importance_pollen, na.rm=TRUE))

village_level_taxa <- visitation_importance_merged_taxa_village %>%   group_by(village_code) %>%  dplyr::summarise(sum_econ_value_village= sum(poll_econ_importance, na.rm=TRUE),
                                                                                                    sum_nutr_value_village= sum(poll_nutr_importance, na.rm=TRUE),
                                                                                                    sum_econ_value_pollen_village= sum(poll_econ_importance_pollen, na.rm=TRUE),
                                                                                                    sum_nutr_value_pollen_village= sum(poll_nutr_importance_pollen, na.rm=TRUE))


#Merge village-level info into species-level, genus-level and taxa-level files
village_poll_sp <- merge(visitation_importance_merged_OTU_village, village_level_sp, by = "village_code")
village_poll_genus <- merge(visitation_importance_merged_genus_village, village_level_genus, by = "village_code")
village_poll_taxa <- merge(visitation_importance_merged_taxa_village, village_level_taxa, by = "village_code")


#----Summarise data at different taxonomic levels
pollinator_importance_sp_village <- village_poll_sp %>%   group_by(insect_OTU, wild_managed, village_code) %>%  reframe(sum_econ_value= (sum(poll_econ_importance, na.rm=TRUE)/sum_econ_value_village),
                                                                                                                sum_nutr_value= (sum(poll_nutr_importance, na.rm=TRUE)/sum_nutr_value_village),
                                                                                                                sum_econ_value_pollen= (sum(poll_econ_importance_pollen, na.rm=TRUE)/sum_econ_value_pollen_village),
                                                                                                                sum_nutr_value_pollen= (sum(poll_nutr_importance_pollen, na.rm=TRUE)/sum_nutr_value_pollen_village))

pollinator_importance_genus_village <- village_poll_genus %>%   group_by(insect_genus, wild_managed, village_code) %>%  reframe(sum_econ_value= (sum(poll_econ_importance, na.rm=TRUE)/sum_econ_value_village),
                                                                                                                        sum_nutr_value= (sum(poll_nutr_importance, na.rm=TRUE)/sum_nutr_value_village),
                                                                                                                        sum_econ_value_pollen= (sum(poll_econ_importance_pollen, na.rm=TRUE)/sum_econ_value_pollen_village),
                                                                                                                        sum_nutr_value_pollen= (sum(poll_nutr_importance_pollen, na.rm=TRUE)/sum_nutr_value_pollen_village))

pollinator_importance_taxa_village <- village_poll_taxa %>%   group_by(pollinator_taxa, wild_managed, village_code) %>%  reframe(sum_econ_value= (sum(poll_econ_importance, na.rm=TRUE)/sum_econ_value_village),
                                                                                                                         sum_nutr_value= (sum(poll_nutr_importance, na.rm=TRUE)/sum_nutr_value_village),
                                                                                                                         sum_econ_value_pollen= (sum(poll_econ_importance_pollen, na.rm=TRUE)/sum_econ_value_pollen_village),
                                                                                                                         sum_nutr_value_pollen= (sum(poll_nutr_importance_pollen, na.rm=TRUE)/sum_nutr_value_pollen_village))

pollinator_importance_wild_managed_village <- village_poll_taxa %>%   group_by(wild_managed, village_code) %>%  reframe(sum_econ_value= (sum(poll_econ_importance, na.rm=TRUE)/sum_econ_value_village),
                                                                                                                sum_nutr_value= (sum(poll_nutr_importance, na.rm=TRUE)/sum_nutr_value_village),
                                                                                                                sum_econ_value_pollen= (sum(poll_econ_importance_pollen, na.rm=TRUE)/sum_econ_value_pollen_village),
                                                                                                                sum_nutr_value_pollen= (sum(poll_nutr_importance_pollen, na.rm=TRUE)/sum_nutr_value_pollen_village))


#Remove duplicated rows
pollinator_importance_sp_village <- distinct(pollinator_importance_sp_village)
pollinator_importance_genus_village <- distinct(pollinator_importance_genus_village)
pollinator_importance_taxa_village <- distinct(pollinator_importance_taxa_village)
pollinator_importance_wild_managed_village <- distinct(pollinator_importance_wild_managed_village)

#Remove all entries with no data
pollinator_importance_sp_village <- pollinator_importance_sp_village[!(pollinator_importance_sp_village$insect_OTU == ""), ]
pollinator_importance_genus_village <- pollinator_importance_genus_village[!(pollinator_importance_genus_village$insect_genus == ""), ]
pollinator_importance_taxa_village <- pollinator_importance_taxa_village[!(pollinator_importance_taxa_village$pollinator_taxa == ""), ]


#################################################################
#######      Identifying best plants under each scenario    ########
#################################################################

#Subset visitation data to leave only plant info and insect OTU
poll_visitation_plant_subset_village <- select(poll_visitation_data, village_code, specimen_code,	pollinator_taxa,	insect_order,	insect_family,	insect_genus,	insect_species,	insect_OTU,	plant_family,	plant_sci_name,	plant_eng_name,	plant_category,	plant_barcode)

#Create new column with concatenated insect OTU and village code for both datasets
pollinator_importance_sp_village$village_insect <- paste(pollinator_importance_sp_village$village_code,"_",pollinator_importance_sp_village$insect_OTU)
poll_visitation_plant_subset_village$village_insect <- paste(poll_visitation_plant_subset_village$village_code,"_",poll_visitation_plant_subset_village$insect_OTU)
poll_visitation_plant_subset_village <- subset(poll_visitation_plant_subset_village, select = -c(insect_OTU, village_code))

##Merge visitation data with pollinator importance data using appropriate merge level 
plant_importance_merged_village<-merge(x=pollinator_importance_sp_village, y=poll_visitation_plant_subset_village,by="village_insect",all.x=TRUE)

#Remove all entries with no plant species names or no village name
plant_importance_merged_village <- plant_importance_merged_village[!(plant_importance_merged_village$plant_sci_name == ""), ]
plant_importance_merged_village <- plant_importance_merged_village %>% drop_na(village_code)

#----dplyr::summarise data at the village level to work our proportion of visits/pollen contributed by each taxon
village_level_plants <- plant_importance_merged_village %>%   group_by(village_code) %>%  dplyr::summarise(sum_econ_value_village= sum(sum_econ_value, na.rm=TRUE),
                                                                                            sum_nutr_value_village= sum(sum_nutr_value, na.rm=TRUE),
                                                                                            sum_econ_value_pollen_village= sum(sum_econ_value_pollen, na.rm=TRUE),
                                                                                            sum_nutr_value_pollen_village= sum(sum_nutr_value_pollen, na.rm=TRUE))

#Merge village-level info into plant importance file
village_plants <- merge(plant_importance_merged_village, village_level_plants, by = "village_code")

#----dplyr::summarise data at the plant level#
plant_importance_village <- village_plants %>%   group_by(village_code, plant_family,plant_sci_name,plant_eng_name,plant_category) %>%  reframe(sum_econ_value= (sum(sum_econ_value, na.rm=TRUE)/sum_econ_value_village),
                                                                                                                                        sum_nutr_value= (sum(sum_nutr_value, na.rm=TRUE)/sum_nutr_value_village),
                                                                                                                                        sum_econ_value_pollen= (sum(sum_econ_value_pollen, na.rm=TRUE)/sum_econ_value_pollen_village),
                                                                                                                                        sum_nutr_value_pollen= (sum(sum_nutr_value_pollen, na.rm=TRUE)/sum_nutr_value_pollen_village))
#Remove duplicated rows
plant_importance_village <- distinct(plant_importance_village)


#----Remove crop plants from plant importance scores#
plant_importance_village_wild <- plant_importance_village[!(plant_importance_village$plant_category == "crop"), ]

#################################################################
#######      Printing out results    ########
#################################################################

#----Export importance summaries#
wb <- createWorkbook()
addWorksheet(wb, sheet = "Info", gridLines = TRUE)
addWorksheet(wb, sheet = "Raw", gridLines = TRUE)
addWorksheet(wb, sheet = "pollinator_importance_sp", gridLines = TRUE)
addWorksheet(wb, sheet = "pollinator_importance_genus", gridLines = TRUE)
addWorksheet(wb, sheet = "pollinator_importance_taxa", gridLines = TRUE)
addWorksheet(wb, sheet = "plant_importance", gridLines = TRUE)
writeData(wb, sheet = "Raw", x = visitation_importance_merged_OTU_village, rowNames = FALSE)
writeData(wb, sheet = "pollinator_importance_sp", x = pollinator_importance_sp_village, rowNames = FALSE)
writeData(wb, sheet = "pollinator_importance_genus", x = pollinator_importance_genus_village, rowNames = FALSE)
writeData(wb, sheet = "pollinator_importance_taxa", x = pollinator_importance_taxa_village, rowNames = FALSE)
writeData(wb, sheet = "plant_importance", x = plant_importance_village, rowNames = FALSE)
freezePane(wb, sheet = "Info", firstRow = TRUE)
freezePane(wb, sheet = "Raw", firstRow = TRUE)
freezePane(wb, sheet = "pollinator_importance_sp", firstRow = TRUE)
freezePane(wb, sheet = "pollinator_importance_genus", firstRow = TRUE)
freezePane(wb, sheet = "pollinator_importance_taxa", firstRow = TRUE)
freezePane(wb, sheet = "plant_importance", firstRow = TRUE)
saveWorkbook(wb, file.path(output.path, paste0("MP_poll-importance_analysis_village_level_", Sys.Date(), ".xlsx")), overwrite = TRUE)


###############################################################################################
#######     Figure 4: Plotting village-level importance  of each pollinator guild        ########
###############################################################################################

# Calculate mean values for each category
mean_econ_values_taxa <- pollinator_importance_taxa_village %>%
  group_by(pollinator_taxa) %>%
  dplyr::summarize(mean_value = mean(sum_econ_value_pollen))

mean_nutr_values_taxa <- pollinator_importance_taxa_village %>%
  group_by(pollinator_taxa) %>%
  summarize(mean_value = mean(sum_nutr_value_pollen))

# Arrange categorical levels based on mean economic values
pollinator_importance_taxa_village$pollinator_taxa <- factor(pollinator_importance_taxa_village$pollinator_taxa, levels = mean_econ_values_taxa$pollinator_taxa[order(mean_econ_values_taxa$mean_value, decreasing = TRUE)])

#Plot economic importance of pollinator taxa
taxa_importance_village_econ <- ggplot(pollinator_importance_taxa_village, aes(x = pollinator_taxa, y = sum_econ_value_pollen)) + 
  geom_boxplot(width = 0.4, fill = "lightgray", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 1, color = "#ee7621ff", size = 1) +
  labs(x = "Pollinator taxon", y = "") +
  theme_minimal() +  ylim(0, 0.8) +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange categorical levels based on mean nutritional values
pollinator_importance_taxa_village$pollinator_taxa <- factor(pollinator_importance_taxa_village$pollinator_taxa, levels = mean_nutr_values_taxa$pollinator_taxa[order(mean_nutr_values_taxa$mean_value, decreasing = TRUE)])

#Plot nutritional importance of pollinator taxa
taxa_importance_village_nutr <- ggplot(pollinator_importance_taxa_village, aes(x = pollinator_taxa, y = sum_nutr_value_pollen)) +
  geom_boxplot(width = 0.4, fill = "lightgray", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 1, color = "#8fa33fff", size = 1) +
  labs(x = "", y = "") +
  theme_minimal() + ylim(0, 0.8) +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot Apis cerena versus wild pollinators for each priority
domestic_importance_village_econ <- ggplot(pollinator_importance_wild_managed_village, aes(x = wild_managed, y = sum_econ_value_pollen)) + 
  geom_boxplot(width = 0.4, fill = "lightgray", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 1, color = "#ee7621ff", size = 1) +
  labs(x = "", y = "Economic importance score") +
  theme_minimal() + ylim(0, 0.8) +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, face="bold"))

domestic_importance_village_nutr <- ggplot(pollinator_importance_wild_managed_village, aes(x = wild_managed, y = sum_nutr_value_pollen)) +
  geom_boxplot(width = 0.4, fill = "lightgray", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 1, color = "#8fa33fff", size = 1) +
  labs(x = "", y = "Nutritional importance score") +
  theme_minimal() + ylim(0, 0.8) +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, face="bold"))

#Join all plots
poll_importance_village_plot <- plot_grid(  domestic_importance_village_nutr, taxa_importance_village_nutr, domestic_importance_village_econ,  taxa_importance_village_econ,   
                                           nrow = 2, labels="", rel_widths = c(0.3, 1, 0.3, 1))


ggsave(plot=poll_importance_village_plot, filename= file.path(output.path, "Poll_importance_village.svg"), width=6, height=8, dpi=500, bg="white")
ggsave(plot=poll_importance_village_plot, filename= file.path(output.path, "Poll_importance_village.png"), width=6, height=8, dpi=500, bg="white")

#######################################################################################
#######    Figure S4: Plotting village-level insect genus importance          ########
#######################################################################################

##Define order of genera to plot (based on data from all villages combined)
genus_importance_order <- pollinator_importance_genus

# Calculate mean values for each category
mean_econ_values_genus <- pollinator_importance_genus_village %>%
  group_by(insect_genus) %>%
  summarize(mean_value = mean(sum_econ_value_pollen))

mean_nutr_values_genus <- pollinator_importance_genus_village %>%
  group_by(insect_genus) %>%
  summarize(mean_value = mean(sum_nutr_value_pollen))

##Subset to the top 40 most important pollinators
important_genera_top40 <- mean_econ_values_genus %>%  top_n(40, mean_value)
pollinator_importance_genus_village <- subset(pollinator_importance_genus_village, insect_genus %in% important_genera_top40$insect_genus)

# Arrange categorical levels based on mean economic values
pollinator_importance_genus_village$insect_genus <- factor(pollinator_importance_genus_village$insect_genus, levels = genus_importance_order$insect_genus[order(genus_importance_order$overall_importance_pollen, decreasing = TRUE)])

# Calculate the range for each category
econ_range_data <- pollinator_importance_genus_village %>%
  group_by(insect_genus) %>%
  summarize(min = min(sum_econ_value_pollen), max = max(sum_econ_value_pollen))

nutr_range_data <- pollinator_importance_genus_village %>%
  group_by(insect_genus) %>%
  summarize(min = min(sum_nutr_value_pollen), max = max(sum_nutr_value_pollen))


#Plot economic importance of pollinator taxa
genus_importance_village_econ <- ggplot(pollinator_importance_genus_village, aes(x = insect_genus, y = sum_econ_value_pollen)) + 
  geom_jitter(width = 0.1, alpha = 1, color = "#ee7621ff", size = 1) +
  labs(x = "", y = "Economic importance score") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  geom_segment(data = econ_range_data, aes(x = insect_genus, xend = insect_genus, y = min, yend = max),
               color = "black")

# Arrange categorical levels based on mean economic values
pollinator_importance_genus_village$insect_genus <- factor(pollinator_importance_genus_village$insect_genus, levels = genus_importance_order$insect_genus[order(genus_importance_order$overall_importance_pollen, decreasing = TRUE)])


#Plot nutritional importance of pollinator taxa
genus_importance_village_nutr <- ggplot(pollinator_importance_genus_village, aes(x = insect_genus, y = sum_nutr_value_pollen)) +
  geom_jitter(width = 0.1, alpha = 1, color = "#8fa33fff", size = 1) +
  labs(x = "Pollinator genus", y = "Nutritional importance score") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6))+
  geom_segment(data = nutr_range_data, aes(x = insect_genus, xend = insect_genus, y = min, yend = max),
               color = "black")


poll_genus_importance_village_plot <- plot_grid(genus_importance_village_econ, genus_importance_village_nutr,  nrow = 2, labels="auto")


ggsave(plot=poll_genus_importance_village_plot, filename= file.path(output.path, "Poll_genus_importance_village_top40.svg"), width=6, height=8, dpi=500, bg="white")
ggsave(plot=poll_genus_importance_village_plot, filename= file.path(output.path, "Poll_genus_importance_village_top40.png"), width=6, height=8, dpi=500, bg="white")


#######################################################################################
#######            Plotting village-level plant importance          ########
#######################################################################################

# Calculate mean values for each category
mean_econ_values_plant <- plant_importance_wild %>%
  group_by(plant_sci_name) %>%
  summarize(mean_value = mean(sum_econ_value_pollen))

mean_nutr_values_plant <- plant_importance_wild %>%
  group_by(plant_sci_name) %>%
  summarize(mean_value = mean(sum_nutr_value_pollen))

important_plants_top40 <- mean_econ_values_plant %>%  top_n(40, mean_value)
plant_importance_subset <- subset(plant_importance_village_wild, plant_sci_name %in% important_plants_top40$plant_sci_name)


# Arrange categorical levels based on mean economic values
plant_importance_subset$plant_sci_name <- factor(plant_importance_subset$plant_sci_name, levels = mean_econ_values_plant$plant_sci_name[order(mean_econ_values_plant$mean_value, decreasing = TRUE)])

# Calculate the range for each category
econ_range_data <- plant_importance_subset %>%
  group_by(plant_sci_name) %>%
  summarize(min = min(sum_econ_value_pollen), max = max(sum_econ_value_pollen))

nutr_range_data <- plant_importance_subset %>%
  group_by(plant_sci_name) %>%
  summarize(min = min(sum_nutr_value_pollen), max = max(sum_nutr_value_pollen))


#Plot economic importance of plant taxa
plant_importance_village_econ <- ggplot(plant_importance_subset, aes(x = plant_sci_name, y = sum_econ_value_pollen)) + 
  geom_jitter(width = 0.1, alpha = 1, color = "#ee7621ff", size = 1) +
  labs(x = "", y = "Economic importance score") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1))+
  geom_segment(data = econ_range_data, aes(x = plant_sci_name, xend = plant_sci_name, y = min, yend = max),
               color = "black")


#Plot nutritional importance of plant taxa
plant_importance_village_nutr <- ggplot(plant_importance_subset, aes(x = plant_sci_name, y = sum_nutr_value_pollen)) +
  geom_jitter(width = 0.1, alpha = 1, color = "#8fa33fff", size = 1) +
  labs(x = "Plant species", y = "Nutritional importance score") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1))+
  geom_segment(data = nutr_range_data, aes(x = plant_sci_name, xend = plant_sci_name, y = min, yend = max),
               color = "black")


plant_importance_village_plot <- plot_grid(plant_importance_village_econ, plant_importance_village_nutr,  nrow = 2, labels="auto")

ggsave(plot=plant_importance_village_plot, filename= file.path(output.path, "Plant_importance_village.svg"), width=7, height=8, dpi=500, bg="white")
ggsave(plot=plant_importance_village_plot, filename= file.path(output.path, "Plant_importance_village.png"), width=7, height=8, dpi=500, bg="white")


###############################################################################################################################################
####### Figure 5:  Calculating and plotting Shannon diversity scores for the pollinator community underpinning food production and income ########
###############################################################################################################################################

#Replace NA values with 0 to avoid errors in diversity calculation
pollinator_importance_sp_village <- pollinator_importance_sp_village %>% replace(is.na(.), 0)

poll_econ_visits <- dcast(pollinator_importance_sp_village, village_code~insect_OTU, value.var="sum_econ_value")
poll_nutr_visits <- dcast(pollinator_importance_sp_village, village_code~insect_OTU, value.var="sum_nutr_value")
poll_econ_pollen <- dcast(pollinator_importance_sp_village, village_code~insect_OTU, value.var="sum_econ_value_pollen")
poll_nutr_pollen <- dcast(pollinator_importance_sp_village, village_code~insect_OTU, value.var="sum_nutr_value_pollen")

poll_econ_visits <- poll_econ_visits %>%   column_to_rownames("village_code")
poll_nutr_visits <- poll_nutr_visits %>%   column_to_rownames("village_code")
poll_econ_pollen <- poll_econ_pollen %>%   column_to_rownames("village_code")
poll_nutr_pollen <- poll_nutr_pollen %>%   column_to_rownames("village_code")

poll_econ_visits <- poll_econ_visits %>% replace(is.na(.), 0)
poll_nutr_visits <- poll_nutr_visits %>% replace(is.na(.), 0)
poll_econ_pollen <- poll_econ_pollen %>% replace(is.na(.), 0)
poll_nutr_pollen <- poll_nutr_pollen %>% replace(is.na(.), 0)

#Calculating Shannon diversity index for each metric
poll_econ_div_visits <- diversity(poll_econ_visits)
poll_nutr_div_visits <- diversity(poll_nutr_visits)
poll_econ_div_pollen <- diversity(poll_econ_pollen)
poll_nutr_div_pollen <- diversity(poll_nutr_pollen)


#Restructuring results to allow analysis
poll_econ_div_visits <- data.frame(poll_econ_div_visits)
poll_econ_div_visits$village_code <- rownames(poll_econ_div_visits)
poll_econ_div_visits$priority <- "economic"
poll_econ_div_visits$data_source <- "visitation"
rownames(poll_econ_div_visits) <- NULL
colnames(poll_econ_div_visits)[1] = "diversity_score"

poll_nutr_div_visits <- data.frame(poll_nutr_div_visits)
poll_nutr_div_visits$village_code <- rownames(poll_nutr_div_visits)
poll_nutr_div_visits$priority <- "nutrition"
poll_nutr_div_visits$data_source <- "visitation"
rownames(poll_nutr_div_visits) <- NULL
colnames(poll_nutr_div_visits)[1] = "diversity_score"

poll_econ_div_pollen <- data.frame(poll_econ_div_pollen)
poll_econ_div_pollen$village_code <- rownames(poll_econ_div_pollen)
poll_econ_div_pollen$priority <- "economic"
poll_econ_div_pollen$data_source <- "pollen"
rownames(poll_econ_div_pollen) <- NULL
colnames(poll_econ_div_pollen)[1] = "diversity_score"

poll_nutr_div_pollen <- data.frame(poll_nutr_div_pollen)
poll_nutr_div_pollen$village_code <- rownames(poll_nutr_div_pollen)
poll_nutr_div_pollen$priority <- "nutrition"
poll_nutr_div_pollen$data_source <- "pollen"
rownames(poll_nutr_div_pollen) <- NULL
colnames(poll_nutr_div_pollen)[1] = "diversity_score"

visitation_results <- rbind(poll_econ_div_visits, poll_nutr_div_visits)
pollen_results <- rbind(poll_econ_div_pollen, poll_nutr_div_pollen)

## Test for significant difference between nutrition and economic-based importance
#Visitation data only
visitation.aov <- aov(diversity_score ~ priority, data = visitation_results)
summary(visitation.aov)

#Pollen-weighted visitation data
pollen.aov <- aov(diversity_score ~ priority, data = pollen_results)
summary(pollen.aov)

## Plot differences between nutrition and economic-based importance using pollen-weighted data
sum_pollen_data = summarySE(pollen_results, measurevar = "diversity_score", groupvars = c("priority"))
pollen_results$priority <- factor(pollen_results$priority, levels = c( "nutrition","economic"), labels = c("Food production", "Income generation"))
sum_pollen_data$priority <- factor(sum_pollen_data$priority, levels = c( "nutrition","economic"), labels = c("Food production", "Income generation"))

## Figure 5 - Shannon diversity of 
pollen_plot <- sum_pollen_data %>%  
  ggplot(aes(x=priority, y=diversity_score, color=priority)) + geom_point(size=5, shape=9, stroke=0.7, color="black") + #This bit plots the mean value for each habitat
  geom_jitter(data=pollen_results, mapping = aes(x=priority, y=diversity_score), size=4, width = 0.15) + 
  geom_point(data=sum_pollen_data, mapping = aes(x=priority, y=diversity_score), size=5, shape=9, stroke=0.7, color="black") +
  #This adds the actual real data points with a bit of sideways 'jitter' so that they aren't all overlapping
  theme_classic() + ylab("Shannon Diversity Score") + xlab("Agricultural Priority") + #This sets the plot theme (black and white) and labels the x and y axis manually so you can change the labels if you want
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),
        plot.title = element_text(color="black", size=14, face="bold", hjust=0.5))+
  ylim(1.3, 3.5)

poll_importance_diversity_plot <- pollen_plot + scale_colour_manual(values = c("#8fa33fff","#ee7621ff"))+
  theme(legend.position = "none")

ggsave(plot=poll_importance_diversity_plot, filename= file.path(output.path, "Pollinator_importance_diversity.svg"), width=5, height=4, dpi=500)
ggsave(plot=poll_importance_diversity_plot, filename= file.path(output.path, "Pollinator_importance_diversity.png"), width=5, height=4, dpi=500)


#############################################################################################################################################
###########################################################    END OF SCRIPT  ###############################################################
#############################################################################################################################################
