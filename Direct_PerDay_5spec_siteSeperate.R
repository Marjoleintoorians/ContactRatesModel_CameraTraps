# Marjolein Toorians
# Calculate interspecific contact rates
# Dec 15, 2023
# DIRECT PER DAY PER SITE

library(dplyr)
library(ggplot2)
library(hms)
library(lubridate)
library(tidyverse)


##################### LOADING DATA ############################

# Load the CT_Select dataset
setwd("/Users/marjoleintoorians/Library/CloudStorage/OneDrive-Personal/UBC_workonOneDrive/Modelling/Ch4 -SI model/CT data analysis and contact rates/Contact rates/")
# A: 5 species only
#CT_direct <- read.csv("CT_direct_persite.csv")
# B: All 23 species
CT_direct <- read.csv("CT_direct_persite_allspecies.csv")

### Pick species! 
# 1: Load species tree
#setwd("/Users/marjoleintoorians/Library/CloudStorage/OneDrive-Personal/UBC_workonOneDrive/Modelling/Ch4 -SI model/Phylogenetic tree/")
#species_table <- read.csv("5_species.csv", header = TRUE)
#species_binom <- species_table$Binomial_name
#species_vec <- species_table$Common_name

# 2: all species
species_vec <- unique(CT_direct$species)
nspec <- length(species_vec)


################### set up matrices ###################

contact_matrix_temp <- matrix(0,length(species_vec),length(species_vec))
rownames(contact_matrix_temp) <- species_vec
colnames(contact_matrix_temp) <- species_vec
contact_matrix_direct_site <- matrix(0,length(species_vec),length(species_vec))
rownames(contact_matrix_direct_site) <- species_vec
colnames(contact_matrix_direct_site) <- species_vec
Mean_direct_day_site <- matrix(0,length(species_vec),length(species_vec))
rownames(Mean_direct_day_site) <- species_vec
colnames(Mean_direct_day_site) <- species_vec
Total_direct_day_site <- matrix(0,length(species_vec),length(species_vec))
rownames(Total_direct_day_site) <- species_vec
colnames(Total_direct_day_site) <- species_vec

# Pick the number of sites
# all 9 sites:
unique_sites <- unique(CT_direct$site)
# for the 5 species, delte DUK and SUN
# unique_sites <- c("De Laporte","Gomondwane","Kwaggas","Nwas","Hoyo Hoyo","Ngotso North","Nyamahri")
nsites <- length(unique_sites) 

################## Manually pick site #####################

site_temp <- unique_sites[1]                                         ### !!! ###
# Create a dataframe for this specific site
sitedf <- subset(CT_direct, site == site_temp)
date_vec <- unique(sitedf$date)

# Set up the matrix for calc the mean for this site
calcDimDF <- sitedf %>% group_by(date) %>% 
  summarize(n_present = n())
dim <- nrow(calcDimDF)
conMatrix_3D <- array(rep(0, nspec*nspec*dim), dim=c(nspec, nspec, dim))

################### Loop #############################

# loop over dates
for(d in 1:length(date_vec)){
  # Subset dates
  datedf <- filter(sitedf, date == date_vec[d])
  #try
  #datedf <- filter(sitedf, date == date_vec[2])
  
  # Initialize empty contact_matrix_temp for this day
  contact_matrix_temp <- matrix(0,length(species_vec),length(species_vec))
  rownames(contact_matrix_temp) <- species_vec
  colnames(contact_matrix_temp) <- species_vec
  
  # find unique photos
  unique_photos <- unique(datedf$timestamp)
  
  # 3 Make a subset df for each photo:
  for(i in 1:length(unique_photos)){
    # subset
    photo_df <- subset(datedf,timestamp == unique_photos[i])
    # try
    #photo_df <- subset(datedf,timestamp == "2015-06-18 15:49")
    
    # count nr of species
    unique_species <- unique(photo_df$species)
    
    # INTRA-specific for each, diagonal. Add n_present-1 to the index
    for(s in 1:length(unique_species)){
      # safe the species name and add to contact matrix
      spec <- unique_species[s]
      photo_df_spec <-  subset(photo_df,species == spec) # df with only 1 row
      contact_matrix_temp[spec,spec] <- contact_matrix_temp[spec,spec] + (photo_df_spec$n_present-1)
    }
    
    # 2 species interspecific direct rates
    if(length(unique_species)==2){
      # names
      spec_1 <- unique_species[1]
      spec_2 <- unique_species[2]
      # collect n_present per species
      x12 <- subset(photo_df,species == spec_1)$n_present # indvs of spec 1 in frame, so per capita contact rate spec 2 has with 1
      x21 <- subset(photo_df,species == spec_2)$n_present # indvs of spec 2 in frame, so per capita contact rate spec 1 has with 2
      # Add the observed to the contact matrix
      contact_matrix_temp[spec_1,spec_2] <- contact_matrix_temp[spec_1,spec_2] + x12
      contact_matrix_temp[spec_2,spec_1] <- contact_matrix_temp[spec_2,spec_1] + x21        
    }
    
    # 3 species interspecific direct rates
    if(length(unique_species)==3){
      spec_1 <- unique_species[1]
      spec_2 <- unique_species[2]
      spec_3 <- unique_species[3]
      # n_present
      x12 <- subset(photo_df,species == spec_1)$n_present
      x13 <- subset(photo_df,species == spec_1)$n_present
      x21 <- subset(photo_df,species == spec_2)$n_present
      x23 <- subset(photo_df,species == spec_2)$n_present
      x31 <- subset(photo_df,species == spec_3)$n_present
      x32 <- subset(photo_df,species == spec_3)$n_present
      # Add the observed to the contact matrix
      contact_matrix_temp[spec_1,spec_2] <- contact_matrix_temp[spec_1,spec_2] + x12
      contact_matrix_temp[spec_1,spec_3] <- contact_matrix_temp[spec_1,spec_3] + x13
      contact_matrix_temp[spec_2,spec_1] <- contact_matrix_temp[spec_2,spec_1] + x21  
      contact_matrix_temp[spec_2,spec_3] <- contact_matrix_temp[spec_2,spec_3] + x23
      contact_matrix_temp[spec_3,spec_1] <- contact_matrix_temp[spec_3,spec_1] + x31  
      contact_matrix_temp[spec_3,spec_2] <- contact_matrix_temp[spec_3,spec_2] + x32
    }
    
    # 4 species interspecific direct rates
    if(length(unique_species)==4){
      spec_1 <- unique_species[1]
      spec_2 <- unique_species[2]
      spec_3 <- unique_species[3]
      spec_4 <- unique_species[4]
      # n_present
      x12 <- subset(photo_df,species == spec_1)$n_present
      x13 <- subset(photo_df,species == spec_1)$n_present
      x14 <- subset(photo_df,species == spec_1)$n_present
      x21 <- subset(photo_df,species == spec_2)$n_present
      x23 <- subset(photo_df,species == spec_2)$n_present
      x24 <- subset(photo_df,species == spec_2)$n_present
      x31 <- subset(photo_df,species == spec_3)$n_present
      x32 <- subset(photo_df,species == spec_3)$n_present
      x34 <- subset(photo_df,species == spec_3)$n_present
      x41 <- subset(photo_df,species == spec_4)$n_present
      x42 <- subset(photo_df,species == spec_4)$n_present
      x43 <- subset(photo_df,species == spec_4)$n_present
      
      # Add the observed to the contact matrix
      contact_matrix_temp[spec_1,spec_2] <- contact_matrix_temp[spec_1,spec_2] + x12
      contact_matrix_temp[spec_1,spec_3] <- contact_matrix_temp[spec_1,spec_3] + x13
      contact_matrix_temp[spec_1,spec_4] <- contact_matrix_temp[spec_1,spec_4] + x14
      contact_matrix_temp[spec_2,spec_1] <- contact_matrix_temp[spec_2,spec_1] + x21  
      contact_matrix_temp[spec_2,spec_3] <- contact_matrix_temp[spec_2,spec_3] + x23
      contact_matrix_temp[spec_2,spec_4] <- contact_matrix_temp[spec_2,spec_4] + x24
      contact_matrix_temp[spec_3,spec_1] <- contact_matrix_temp[spec_3,spec_1] + x31  
      contact_matrix_temp[spec_3,spec_2] <- contact_matrix_temp[spec_3,spec_2] + x32
      contact_matrix_temp[spec_3,spec_4] <- contact_matrix_temp[spec_3,spec_4] + x34
      contact_matrix_temp[spec_4,spec_1] <- contact_matrix_temp[spec_4,spec_1] + x41  
      contact_matrix_temp[spec_4,spec_2] <- contact_matrix_temp[spec_4,spec_2] + x42
      contact_matrix_temp[spec_4,spec_3] <- contact_matrix_temp[spec_4,spec_3] + x43
    }
    
  } # end timestamp (i)
  
  # save the contact matrix at end of this day
  conMatrix_3D[,,d] <- contact_matrix_temp # Matrix for each day
  
} # end date (d)

###################### save per site ###########################

# Calculate the average contact rate over of all days, of the days # matrices.
for(x in 1:nspec){
  for(y in 1:nspec){
    contact_matrix_direct_site[x,y] <- mean(conMatrix_3D[x,y,])
  }
}

print(contact_matrix_direct_site)
print(site_temp)


# save each site matrix                                             ### !!! ###
GOM <- contact_matrix_direct_site
write.csv(GOM,"GOM.csv")


################# save all separate sites #################
setwd("/Users/marjoleintoorians/Library/CloudStorage/OneDrive-Personal/UBC_workonOneDrive/Modelling/Ch4 -SI model/CT data analysis and contact rates/Contact rates/Direct/Contact matrices per site 16dec/")

# list sites:
list_sites <- array(rep(0, nspec*nspec*nsites), dim=c(nspec, nspec, nsites))
list_sites[,,1] <- KWA
list_sites[,,2] <- NGO
list_sites[,,3] <- NYA
list_sites[,,4] <- NWA
list_sites[,,5] <- HOY
list_sites[,,6] <- DLP
list_sites[,,7] <- GOM
# optional
list_sites[,,8] <- SUN
list_sites[,,9] <- DUK


# Merge sites: sum over all sites
for(l in 1:nspec){
  for(m in 1:nspec){
    Total_direct_day_site[l,m] <- sum(list_sites[l,m,])  }}
print(Total_direct_day_site)


# Save total of all sites
write.csv(Total_direct_day_site, "23species_9sites.csv")


