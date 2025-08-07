#### Analysis of food web robustness based on node properties #### 
#### Code developed by Fredric Windsor #### 
#### Please contact fmwindsor@gmail.com for any questions ####

# Apologies in advance, there is a lot of repetition, etc. that I am sure 
# someone more proficient could have avoided in an more elegant way!! 

#### Setup #### 

## Get everything clean and tidy for this script

# Clear the working environment
rm(list=ls())

# Set the working directory
setwd("/Users/c1513054/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Realistic robustness (Ecology Letters)")

# Load the libraries required
library(dplyr); library(magrittr); library(igraph); library(reshape2)
library(ggplot2); library(bipartite); library(gridExtra)
library(rmangal); library(rtry) # Data sources


#### Data wrangling #### 

### Food webs

## Read in data from Brose et al. 2005 

# Get the data
foodwebs <- read.csv("Data/Foodwebs/bodysizes_2008.csv")

# Create two blank columns with clean names for organisms 
foodwebs$resource.name <- NULL; foodwebs$consumer.name <- NULL

# Some networks have taxonomic names, others common names so we need to merge
# the two columns to make life easier later on! 
foodwebs[1:191, "resource.name"] <- foodwebs[1:191, "Common.name.s..resource"]
foodwebs[192:16866, "resource.name"] <- foodwebs[192:16866, "Taxonomy.resource"]
foodwebs[1:191, "consumer.name"] <- foodwebs[1:191,"Common.name.s..consumer"]
foodwebs[192:16866, "consumer.name"] <- foodwebs[192:16866, "Taxonomy.consumer"]

# Some nodes have mass or length so we need to merge the size units to make it
# easier to run through the later script without having to tailor everything to 
# the individual networks
foodwebs[c(1:347, 864:898, 1518:2145, 2480:14343, 14452:14871), "resource.size"] <- foodwebs[c(1:347, 864:898, 1518:2145, 2480:14343, 14452:14871), "Mean.mass..g..resource"]
foodwebs[c(348:863, 899:1517, 2146:2479, 14344:14451, 14872:16866), "resource.size"] <- foodwebs[c(348:863, 899:1517, 2146:2479, 14344:14451, 14872:16866), "Mean.length..m..resource"]
foodwebs[c(1:347, 864:898, 1518:2145, 2480:14343, 14452:14871), "consumer.size"] <- foodwebs[c(1:347, 864:898, 1518:2145, 2480:14343, 14452:14871), "Mean.mass..g..consumer"]
foodwebs[c(348:863, 899:1517, 2146:2479, 14344:14451, 14872:16866), "consumer.size"] <- foodwebs[c(348:863, 899:1517, 2146:2479, 14344:14451, 14872:16866), "Mean.length..m..consumer"]

# There is a zero entry as a resource in one dataset 
foodwebs_clean <- foodwebs[foodwebs$Link.ID != 3010,]

# For the following tests, we are going to focus on freshwater food webs as they
# are strongly size-structured and there is clear theoretical and empirical 
# evidence that large-bodied organisms are lost first under elevated temperature

# Subset to the freshwater networks collected in the field from 
freshwater_foodwebs <- subset(foodwebs_clean, General.habitat == "freshwater" &
                                Specific.habitat != "laboratory study" &
                                Specific.habitat != "freshwater communities")

# Split into the different networks
foodweb_edgelists <- split(freshwater_foodwebs, 
                           f = as.factor(freshwater_foodwebs$Link.reference))

# Change the order of the columns
foodweb_edgelists_correct <- lapply(foodweb_edgelists, 
                                    function(x){dplyr::relocate(x, 
                                                                resource.name,
                                                                consumer.name, 
                                                        .before = Link.ID)})

# Remove one food web with low taxonomic resolution
foodweb_edgelists_sub <- foodweb_edgelists_correct[-c(2, 5)]

# Create networks from the edgelists
foodweb_networks <- lapply(foodweb_edgelists_sub, graph_from_data_frame)

# Create incidence matrices for the networks
foodweb_matrices <- lapply(foodweb_networks, function(x){bipartite::empty(as.matrix(as_adjacency_matrix(x)))})

# Body size scenario
bodysizes <- lapply(foodweb_edgelists_sub, function(x){aggregate(resource.size ~ resource.name, data = x, FUN = mean)})
bodysizes_ext_list <- lapply(bodysizes, function(x){x[base::order(x$resource.size, decreasing = T),]})


#### Robustness (from Bane et al. 2018) ####

# For each network extract the different extinction scenarios
# 1. Based on large to small body size
# 2. Based on most to least connected (out degree to create bottom-up losses)
# 3. Based on least to most connected (out degree to create bottom-up losses)
# 4. Mean from 1000 random extinction scenarios

## Random extinction sequence

# Set up the empty lists
fw_rand_robustnessQ_curves<-NULL; fw_rand_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(foodweb_matrices)){
  myobsmatrix <- t(foodweb_matrices[[k]])
  nperm<-1000    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix          # make copy of original matrix to work on
    resource<-colSums(mymat)    # save resource degrees
    consumer<-rowSums(mymat)    # save consumer degrees
    survivors<-NULL             # create survivors to save consumer counts
    triggerhat<-colnames(mymat) # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-sample(triggerhat,1) # pick a trigger from the hat
      mymat[,trigger]<-0            # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fw_rand_robustnessQ_Rvalues[[k]]<-sRvalues
  fw_rand_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fw_rand_robustnessQ_plot_wide<-NULL
for (p in 1:length(fw_rand_robustnessQ_curves)){
  trial <- fw_rand_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fw_rand_robustnessQ_plot_wide[[p]]<-melt(trial)
}

# Rename the list to be the names of the channels 
names(fw_rand_robustnessQ_plot_wide) <- names(foodweb_networks)

# Bind the results together and label based on the channels
fw_rand_robustnessQ_plot_long <- bind_rows(fw_rand_robustnessQ_plot_wide, .id = "Foodweb")

## Worst case (high to low degree) extinction sequence

# Set up the empty lists
fw_worst_robustnessQ_curves<-NULL; fw_worst_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(foodweb_matrices)){
  myobsmatrix <- t(foodweb_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=TRUE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger 
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fw_worst_robustnessQ_Rvalues[[k]]<-sRvalues
  fw_worst_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fw_worst_robustnessQ_plot_wide<-NULL
for (p in 1:length(fw_worst_robustnessQ_curves)){
  trial <- fw_worst_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fw_worst_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fw_worst_robustnessQ_plot_wide) <- names(foodweb_networks)

# Bind the results together and label based on the channels
fw_worst_robustnessQ_plot_long <- bind_rows(fw_worst_robustnessQ_plot_wide, .id = "Foodweb")

## Best case (high to low degree) extinction sequence

# Set up the empty lists
fw_best_robustnessQ_curves<-NULL; fw_best_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(foodweb_matrices)){
  myobsmatrix <- t(foodweb_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=FALSE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1]   # pick a trigger from the hat
      mymat[,trigger]<-0      # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fw_best_robustnessQ_Rvalues[[k]]<-sRvalues
  fw_best_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fw_best_robustnessQ_plot_wide<-NULL
for (p in 1:length(fw_best_robustnessQ_curves)){
  trial <- fw_best_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fw_best_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fw_best_robustnessQ_plot_wide) <- names(foodweb_networks)

# Bind the results together and label based on the channels
fw_best_robustnessQ_plot_long <- bind_rows(fw_best_robustnessQ_plot_wide, .id = "Foodweb")

## Size dependent extinction sequence

# Set up the empty lists
fw_size_robustnessQ_curves<-NULL; fw_size_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(foodweb_matrices)){
  myobsmatrix <- t(foodweb_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                  # make copy of original matrix to work on
    resource<-colSums(mymat)                            # save resource degrees
    consumer<-rowSums(mymat)                            # save consumer degrees
    res.names<- colnames(mymat)                         # names of the resources
    survivors<-NULL                                     # create survivors to save consumer counts
    chosen.order<-bodysizes_ext_list[[k]]$resource.name # body size based extinction lists                 
    triggerhat<-chosen.order                            # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger  
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fw_size_robustnessQ_Rvalues[[k]]<-sRvalues
  fw_size_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fw_size_robustnessQ_plot_wide<-NULL
for (p in 1:length(fw_size_robustnessQ_curves)){
  trial <- fw_size_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fw_size_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fw_size_robustnessQ_plot_wide) <- names(foodweb_networks)

# Bind the results together and label based on the channels
fw_size_robustnessQ_plot_long <- bind_rows(fw_size_robustnessQ_plot_wide, .id = "Foodweb")

## Wrangle the data to analyse later

# Sort the robustness R value results 
fw_Rvalues_dframe <- data.frame(foodweb = names(fw_rand_robustnessQ_plot_wide),
                                R_rand_mean = unlist(lapply(fw_rand_robustnessQ_Rvalues, mean)),
                                R_rand_sd = unlist(lapply(fw_rand_robustnessQ_Rvalues, sd)),
                                R_best_mean = unlist(lapply(fw_best_robustnessQ_Rvalues, mean)),
                                R_best_sd = unlist(lapply(fw_best_robustnessQ_Rvalues, sd)),
                                R_worst_mean = unlist(lapply(fw_worst_robustnessQ_Rvalues, mean)),
                                R_sd_mean = unlist(lapply(fw_worst_robustnessQ_Rvalues, sd)),
                                R_size_mean = unlist(lapply(fw_size_robustnessQ_Rvalues, mean)),
                                R_size_sd = unlist(lapply(fw_size_robustnessQ_Rvalues, sd)))

# Labels 
foodweb.labs <- c("FW_HS", "FW_LE", "FW_RA")
names(foodweb.labs) <- fw_Rvalues_dframe$foodweb

# ggplot rendering 
plot1 <- ggplot() +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "black", linewidth = 1, inherit.aes = FALSE,
               data = fw_rand_robustnessQ_plot_long) + 
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "darkred", linewidth = 1, inherit.aes = FALSE,
               data = fw_worst_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "seagreen", linewidth = 1, inherit.aes = FALSE,
               data = fw_best_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "orange", linewidth = 1, inherit.aes = FALSE,
               data = fw_size_robustnessQ_plot_long) +
  ylab("Proportion of consumers remaining") + 
  xlab("Proportion of resources removed") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 8, colour = "black"), 
        strip.text = element_text(face = "bold", hjust = 0, size = 10),
        strip.background = element_blank()) + 
  facet_wrap(~ Foodweb, ncol = 3, labeller = labeller(Foodweb = foodweb.labs))
plot1

## Clean directory
rm(mymat, myobsmatrix, resourcesurvivors, trial, chosen.order, consumer, i, j, 
   k, nperm, p, resource, sRvalues, survivors, survivorsx, threshold, trigger,
   triggerhat, v.ext, extlists, exttimes, bodysizes, bodysizes_ext_list, 
   foodweb_edgelists, foodweb_edgelists_correct, foodweb_matrices, res.names,
   foodweb_networks, foodwebs, foodwebs_clean, freshwater_foodwebs,
   fw_best_robustnessQ_plot_wide, fw_rand_robustnessQ_plot_wide,
   fw_worst_robustnessQ_plot_wide, fw_size_robustnessQ_plot_wide,
   fw_rand_robustnessQ_Rvalues, fw_best_robustnessQ_Rvalues, 
   fw_worst_robustnessQ_Rvalues, fw_size_robustnessQ_Rvalues, 
   foodweb_edgelists_sub, fw_best_robustnessQ_curves, 
   fw_rand_robustnessQ_curves, fw_worst_robustnessQ_curves,
   fw_size_robustnessQ_curves, foodweb.labs)


### Plant-pollinator networks

## Trait data 

# UK plant trait data
plantat_traits <- read.csv("Data/Flower visitor networks/PLANTATT_19_Nov_08.csv")
core_traits <- read.csv("Data/Flower visitor networks/CoRRE_continuousTraitData_Nov2023.csv") %>% filter(trait == "SLA")

## Norwood farm plant-pollinator data (Pocock et al. 2012)

# Load in the data (from Windsor et al. 2022 but generated by Pocock et al. 2012)
pocock <- read.csv("Data/Flower visitor networks/pocock2012.csv")
pocock_fv <- subset(pocock, upper.guild == "02FV." | upper.guild == "12BFLY.")

# Replace the names in the norwood data to make sure they match the traits in the heights dataset
pocock_fv[pocock_fv$lower.species == "Conolvulus/Calystegia", "lower.species"] <- "Convolvulus arvensis"
pocock_fv[pocock_fv$lower.species == "Medicago sativa ssp. sativa", "lower.species"] <- "Medicago sativa subsp.sativa"
pocock_fv[pocock_fv$lower.species == "Papaver sp.", "lower.species"] <- "Papaver rhoeas"
pocock_fv[pocock_fv$lower.species == "Rosa sp.", "lower.species"] <- "Rosa arvensis"
pocock_fv[pocock_fv$lower.species == "Rubus fruticosus", "lower.species"] <- "Rubus fruticosus agg."
pocock_fv[pocock_fv$lower.species == "Senacio jacobaea", "lower.species"] <- "Senecio jacobaea"
pocock_fv[pocock_fv$lower.species == "Sonchus asper/oleraceus", "lower.species"] <- "Sonchus asper"
pocock_fv[pocock_fv$lower.species == "Trifolium pratense/repens", "lower.species"] <- "Trifolium arvense"
pocock_fv[pocock_fv$lower.species == "Viburnum lanatum", "lower.species"] <- "Viburnum lantana"
pocock_fv[pocock_fv$lower.species == "Vicia hirsute", "lower.species"] <- "Vicia hirsuta"
pocock_fv[pocock_fv$lower.species == "Epilobium sp.", "lower.species"] <- "Epilobium lanceolatum"
pocock_fv[pocock_fv$lower.species == "Hypochaeris/ Crepis sp.", "lower.species"] <- "Hypochaeris radicata"
pocock_fv[pocock_fv$lower.species == "Taraxacum officinale", "lower.species"] <- "Taraxacum"

# Clean the dataset to remove additional columns
pocock_fv_clean <- pocock_fv %>% 
  dplyr::select(-lower.guild, -upper.guild, -lower, -upper) %>%
  rename(resource.name = lower.species, consumer.name = upper.species, 
         weight = fortotals)

# Create a list of plants and their heights
pocock_heights <- left_join(pocock_fv_clean, plantat_traits, by = c("resource.name" = "Taxon.name")) %>%
  dplyr::select(resource.name, Hght) %>%
  filter(!duplicated(resource.name)) 

# Replace the names in the norwood data to make sure they match the traits in the SLA dataset
pocock_fv[pocock_fv$lower.species == "Brassica napus", "lower.species"] <- "Brassica nigra" # Closely related species
pocock_fv[pocock_fv$lower.species == "Chaerophyllum temulum", "lower.species"] <- "Anthriscus sylvestris" # Closely related species
pocock_fv[pocock_fv$lower.species == "Crataegus monogyna", "lower.species"] <- "Cotoneaster integerrimus" # Same subtribe
pocock_fv[pocock_fv$lower.species == "Galium odoratum", "lower.species"] <- "Galium album" # Morphologically similar
pocock_fv[pocock_fv$lower.species == "Geranium robertianum", "lower.species"] <- "Geranium pratense" # Closely related species
pocock_fv[pocock_fv$lower.species == "Lamiastrum galeobdolon", "lower.species"] <- "Lamium purpureum" # Closely related species
pocock_fv[pocock_fv$lower.species == "Lamium album", "lower.species"] <- "Lamium purpureum" # Closely related species
pocock_fv[pocock_fv$lower.species == "Leontodon autumnalis", "lower.species"] <- "Scorzoneroides autumnalis" # Correct name
pocock_fv[pocock_fv$lower.species == "Ligustrum vulgare", "lower.species"] <- "Jasminum odoratissimum" # Closest related species
pocock_fv[pocock_fv$lower.species == "Matricaria recutita", "lower.species"] <- "Matricaria breviradiata"
pocock_fv[pocock_fv$lower.species == "Medicago sativa subsp.sativa", "lower.species"] <- "Medicago sativa"
pocock_fv[pocock_fv$lower.species == "Mentha arvensis", "lower.species"] <- "Mentha pulegium"
pocock_fv[pocock_fv$lower.species == "Prunus spinosa", "lower.species"] <- "Prunus pumila"
pocock_fv[pocock_fv$lower.species == "Rosa arvensis", "lower.species"] <- "Rosa acicularis"
pocock_fv[pocock_fv$lower.species == "Rubus fruticosus agg.", "lower.species"] <- "Rubus vestitus"
pocock_fv[pocock_fv$lower.species == "Senecio jacobaea", "lower.species"] <- "Jacobaea vulgaris"
pocock_fv[pocock_fv$lower.species == "Silene dioica", "lower.species"] <- "Silene latifolia"
pocock_fv[pocock_fv$lower.species == "Stachys sylvatica", "lower.species"] <- "Stachys sylvatica"
pocock_fv[pocock_fv$lower.species == "Viburnum lantana", "lower.species"] <- "Viburnum nudum"
pocock_fv[pocock_fv$lower.species == "Epilobium lanceolatum", "lower.species"] <- "Epilobium ciliatum"
pocock_fv[pocock_fv$lower.species == "Hyacinthoides non-scripta", "lower.species"] <- "Dipcadi serotinum"
pocock_fv[pocock_fv$lower.species == "Taraxacum", "lower.species"] <- "Taraxacum campylodes"

# Create the database
pocock_sla <- left_join(pocock_fv_clean, core_traits, by = c("resource.name" = "species")) %>%
  dplyr::select(resource.name, trait_value) %>%
  filter(!duplicated(resource.name)) 

## Data from a grassland system (Memmott et al. 1999)

# Get data
data("memmott1999")

# Load in the original data
memmott <- data.frame(lower.species = rownames(memmott1999), memmott1999, 
                      row.names = seq(1, 25, 1))

# Convert the matrix to an edgelist
memmott_fv <- melt(memmott, value = "count") %>% filter(value > 0)

# Replace periods with spaces 
memmott_fv$lower.species <- gsub(".", " ", memmott_fv$lower.species, fixed = TRUE)

# Change the names of some of the plants to match the traits 
memmott_fv[memmott_fv$lower.species == "Agrimonium eupatorium", "lower.species"] <- "Agrimonia eupatoria"
memmott_fv[memmott_fv$lower.species == "Euphrasia officinalis", "lower.species"] <- "Euphrasia officinalis agg."
memmott_fv[memmott_fv$lower.species == "Rubus fruticosus", "lower.species"] <- "Rubus fruticosus agg."

# Clean the dataset to remove additional columns
memmott_fv_clean <- memmott_fv %>%
  rename(resource.name = lower.species, consumer.name = variable, 
         weight = value)

# Create a list of plants and their heights
memmott_heights <- left_join(memmott_fv_clean, plantat_traits, by = c("resource.name" = "Taxon.name")) %>%
  dplyr::select(resource.name, Hght) %>%
  filter(!duplicated(resource.name)) 

memmott_fv[memmott_fv$lower.species == "Leontodon saxatilis", "lower.species"] <- "Leontodon taraxacoides"
memmott_fv[memmott_fv$lower.species == "Leontodon autumnalis", "lower.species"] <- "Scorzoneroides autumnalis" # Correct name
memmott_fv[memmott_fv$lower.species == "Senecio jacobaea", "lower.species"] <- "Jacobaea vulgaris"
memmott_fv[memmott_fv$lower.species == "Aethusa cynapium", "lower.species"] <- "Anthriscus sylvestris"
memmott_fv[memmott_fv$lower.species == "Chamerion angustifolium", "lower.species"] <- "Epilobium ciliatum"
memmott_fv[memmott_fv$lower.species == "Agrimonia eupatoria", "lower.species"] <- "Agrimonia procera"
memmott_fv[memmott_fv$lower.species == "Euphrasia officinalis agg.", "lower.species"] <- "Euphrasia officinalis"
memmott_fv[memmott_fv$lower.species == "Clematis vitalba", "lower.species"] <- "Clematis fremontii"
memmott_fv[memmott_fv$lower.species == "Aethusa cynapium", "lower.species"] <- "Scorzoneroides autumnalis"
memmott_fv[memmott_fv$lower.species == "Rubus fruticosus agg.", "lower.species"] <- "Rubus vestitus"

# Create a list of plants and their heights
memmott_sla <- left_join(memmott_fv_clean, core_traits, by = c("resource.name" = "species")) %>%
  dplyr::select(resource.name, trait_value) %>%
  filter(!duplicated(resource.name)) 


## Agricultural grassland plant-pollinator interactions (Magrach et al. 2017)

# Load the original data
magrach <- read.csv("Data/Flower visitor networks/magrach2017.csv")

# Subset the data to only the 'non-crop' or semi-natural grasslands
magrach_sub <- subset(magrach, country == "UK" & landscape_type == "low")

# Merge the genus and species names for the plants and pollinators
magrach_sub$plants <- paste(magrach_sub$plant_genus, magrach_sub$plant_species, sep = " ")
magrach_sub$pollinators <- paste(magrach_sub$pollinator_genus, magrach_sub$pollinator_species, sep = " ")

# Aggregate interactions across sampling rounds and across sites to make a meta-network
magrach_fv <- aggregate(. ~ plants + pollinators, 
                      dplyr::select(magrach_sub, -pollinator_species.1, 
                                    -country, -year, -landscape_type, -transect,
                                    -period, -plant_species, -plant_genus,
                                    -pollinator_species, -pollinator_genus,
                                    -pollinator_group, -site_id), 
                      FUN = mean)

# Change the names of some of the plants to match the traits 
magrach_fv[magrach_fv$plants == "Hieracium agg.", "plants"] <- "Pilosella flagellaris"
magrach_fv[magrach_fv$plants == "Taraxacum agg.", "plants"] <- "Taraxacum"
magrach_fv[magrach_fv$plants == "Myosotis Unknown", "plants"] <- "Myosotis arvensis" # Middle ground of height for the genus
magrach_fv[magrach_fv$plants == "Hypochaeris Unknown", "plants"] <- "Hypochaeris maculata" # Middle ground of height for the genus
magrach_fv[magrach_fv$plants == "Polygala Unknown", "plants"] <- "Polygala serpyllifolia" # Middle ground of height for the genus
magrach_fv[magrach_fv$plants == "Cruciata verum", "plants"] <- "Cruciata laevipes"
magrach_fv[magrach_fv$plants == "Galium Unknown", "plants"] <- "Galium spurium" # Middle ground of height for the genus

# Remove a plant 'unknown unknown'
magrach_fv <- magrach_fv[magrach_fv$plants != "Unknown Unknown",]

# Clean the dataset to remove additional columns
magrach_fv_clean <- magrach_fv %>%
  dplyr::select(-round) %>%
  rename(resource.name = plants, consumer.name = pollinators)

# Add a weight column
magrach_fv_clean$weight <- rep(1, nrow(magrach_fv_clean))

# Create a list of plants and their heights
magrach_heights <- left_join(magrach_fv_clean, plantat_traits, by = c("resource.name" = "Taxon.name")) %>%
  dplyr::select(resource.name, Hght) %>%
  filter(!duplicated(resource.name))

magrach_fv[magrach_fv$plants == "Pilosella flagellaris", "plants"] <- "Pilosella caespitosa"
magrach_fv[magrach_fv$plants == "Thymus polytrichus", "plants"] <- "Thymus praecox"
magrach_fv[magrach_fv$plants == "Taraxacum", "lower.species"] <- "Taraxacum campylodes"
magrach_fv[magrach_fv$plants == "Clinopodium vulgare", "plants"] <- "Clinopodium acinos"
magrach_fv[magrach_fv$plants == "Senecio jacobaea", "lower.species"] <- "Jacobaea vulgaris"
magrach_fv[magrach_fv$plants == "Centaurea scabiosa", "plants"] <- "Centaurea nigra"
magrach_fv[magrach_fv$plants == "Cirsium acaule", "plants"] <- "Cirsium arvense"
magrach_fv[magrach_fv$plants == "Campanula glomerata", "plants"] <- "Campanula rotundifolia"
magrach_fv[magrach_fv$plants == "Carlina vulgaris", "plants"] <- "Carlina acaulis"
magrach_fv[magrach_fv$plants == "Hippocrepis comosa", "plants"] <- "Hippocrepis multisiliquosa"
magrach_fv[magrach_fv$plants == "Hypochaeris maculata", "plants"] <- "Hypochaeris radicata"
magrach_fv[magrach_fv$plants == "Linaria repens", "plants"] <- "Linaria canadensis"
magrach_fv[magrach_fv$plants == "Lamium album", "lower.species"] <- "Lamium purpureum" # Closely related species
magrach_fv[magrach_fv$plants == "Odontites vernus", "plants"] <- "Bartsia alpina"
magrach_fv[magrach_fv$plants == "Aquilegia vulgaris", "plants"] <- "Aquilegia canadensis"
magrach_fv[magrach_fv$plants == "Agrimonia eupatoria", "lower.species"] <- "Agrimonia procera"
magrach_fv[magrach_fv$plants == "Pastinaca sativa", "plants"] <- "Daucus carota"

# Create the database
magrach_sla <- left_join(magrach_fv_clean, core_traits, by = c("resource.name" = "species")) %>%
  dplyr::select(resource.name, trait_value) %>%
  filter(!duplicated(resource.name)) 

## Create lists of data

# Organise the data into a consistent format with consistent column names

# Create a list of edgelists
fv_edgelists <- list(pocock_fv_clean, memmott_fv_clean, magrach_fv_clean)
names(fv_edgelists) <- c("FV_PO", "FV_ME", "FV_MA")

# Create networks from the edgelists
fv_networks <- lapply(fv_edgelists, graph_from_data_frame)
names(fv_networks) <- c("FV_PO", "FV_ME", "FV_MA")

# Create incidence matrices for the networks
fv_matrices <- lapply(fv_networks, function(x){bipartite::empty(as.matrix(as_adjacency_matrix(x)))})
names(fv_matrices) <- c("FV_PO", "FV_ME", "FV_MA")

# Plant height scenario
plantheights <- list(pocock_heights, memmott_heights,  magrach_heights)
names(plantheights) <- c("FV_PO", "FV_ME", "FV_MA")

# Create the height ordered extinction lists (shorter plants are)
height_ext_list <- lapply(plantheights, function(x){x[base::order(x$Hght, decreasing = T),]})

# Plant sla scenario
plantsla <- list(pocock_sla, memmott_sla, magrach_sla)
names(plantsla) <- c("FV_PO", "FV_ME", "FV_MA")

# Create the height ordered extinction lists (shorter plants are)
sla_ext_list <- lapply(plantsla, function(x){x[base::order(x$trait_value, decreasing = T),]})


#### Robustness (from Bane et al. 2018) ####

# For each network extract the different extinction scenarios
# 1. Based on large to small body size
# 2. Based on most to least connected (out degree to create bottom-up losses)
# 3. Based on least to most connected (out degree to create bottom-up losses)
# 4. Mean from 1000 random extinction scenarios

## Random extinction sequence

# Set up the empty lists
fv_rand_robustnessQ_curves<-NULL; fv_rand_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(fv_matrices)){
  myobsmatrix <- t(fv_matrices[[k]])
  nperm<-1000    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix          # make copy of original matrix to work on
    resource<-colSums(mymat)    # save resource degrees
    consumer<-rowSums(mymat)    # save consumer degrees
    survivors<-NULL             # create survivors to save consumer counts
    triggerhat<-colnames(mymat) # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-sample(triggerhat,1) # pick a trigger from the hat
      mymat[,trigger]<-0            # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fv_rand_robustnessQ_Rvalues[[k]]<-sRvalues
  fv_rand_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fv_rand_robustnessQ_plot_wide<-NULL
for (p in 1:length(fv_rand_robustnessQ_curves)){
  trial <- fv_rand_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fv_rand_robustnessQ_plot_wide[[p]]<-melt(trial)
}

# Rename the list to be the names of the channels 
names(fv_rand_robustnessQ_plot_wide) <- names(fv_networks)

# Bind the results together and label based on the channels
fv_rand_robustnessQ_plot_long <- bind_rows(fv_rand_robustnessQ_plot_wide, .id = "network")


## Worst case (high to low degree) extinction sequence

# Set up the empty lists
fv_worst_robustnessQ_curves<-NULL; fv_worst_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(fv_matrices)){
  myobsmatrix <- t(fv_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=TRUE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger 
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fv_worst_robustnessQ_Rvalues[[k]]<-sRvalues
  fv_worst_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fv_worst_robustnessQ_plot_wide<-NULL
for (p in 1:length(fv_worst_robustnessQ_curves)){
  trial <- fv_worst_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fv_worst_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fv_worst_robustnessQ_plot_wide) <- names(fv_networks)

# Bind the results together and label based on the channels
fv_worst_robustnessQ_plot_long <- bind_rows(fv_worst_robustnessQ_plot_wide, .id = "network")


## Best case (high to low degree) extinction sequence

# Set up the empty lists
fv_best_robustnessQ_curves<-NULL; fv_best_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(fv_matrices)){
  myobsmatrix <- t(fv_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=FALSE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1]   # pick a trigger from the hat
      mymat[,trigger]<-0      # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fv_best_robustnessQ_Rvalues[[k]]<-sRvalues
  fv_best_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fv_best_robustnessQ_plot_wide<-NULL
for (p in 1:length(fv_best_robustnessQ_curves)){
  trial <- fv_best_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fv_best_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fv_best_robustnessQ_plot_wide) <- names(fv_networks)

# Bind the results together and label based on the channels
fv_best_robustnessQ_plot_long <- bind_rows(fv_best_robustnessQ_plot_wide, .id = "network")


## Height dependent extinction sequence

# Set up the empty lists
fv_height_robustnessQ_curves<-NULL; fv_height_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(fv_matrices)){
  myobsmatrix <- t(fv_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                               # make copy of original matrix to work on
    resource<-colSums(mymat)                         # save resource degrees
    consumer<-rowSums(mymat)                         # save consumer degrees
    res.names<- colnames(mymat)                      # names of the resources
    survivors<-NULL                                  # create survivors to save consumer counts
    chosen.order<-height_ext_list[[k]]$resource.name # height based extinction lists                 
    triggerhat<-chosen.order                         # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger  
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fv_height_robustnessQ_Rvalues[[k]]<-sRvalues
  fv_height_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fv_height_robustnessQ_plot_wide<-NULL
for (p in 1:length(fv_height_robustnessQ_curves)){
  trial <- fv_height_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fv_height_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fv_height_robustnessQ_plot_wide) <- names(fv_networks)

# Bind the results together and label based on the channels
fv_height_robustnessQ_plot_long <- bind_rows(fv_height_robustnessQ_plot_wide, .id = "network")


## SLA dependent extinction sequence

# Set up the empty lists
fv_sla_robustnessQ_curves<-NULL; fv_sla_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(fv_matrices)){
  myobsmatrix <- t(fv_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                            # make copy of original matrix to work on
    resource<-colSums(mymat)                      # save resource degrees
    consumer<-rowSums(mymat)                      # save consumer degrees
    res.names<- colnames(mymat)                   # names of the resources
    survivors<-NULL                               # create survivors to save consumer counts
    chosen.order<-sla_ext_list[[k]]$resource.name # sla based extinction lists                 
    triggerhat<-chosen.order                      # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger  
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  fv_sla_robustnessQ_Rvalues[[k]]<-sRvalues
  fv_sla_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
fv_sla_robustnessQ_plot_wide<-NULL
for (p in 1:length(fv_sla_robustnessQ_curves)){
  trial <- fv_sla_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  fv_sla_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(fv_sla_robustnessQ_plot_wide) <- names(fv_networks)

# Bind the results together and label based on the channels
fv_sla_robustnessQ_plot_long <- bind_rows(fv_sla_robustnessQ_plot_wide, .id = "network")

## Wrangle the data to analyse later

# Sort the robustness R value results 
fv_Rvalues_dframe <- data.frame(network = names(fv_matrices),
                                R_rand_mean = unlist(lapply(fv_rand_robustnessQ_Rvalues, mean)),
                                R_rand_sd = unlist(lapply(fv_rand_robustnessQ_Rvalues, sd)),
                                R_best_mean = unlist(lapply(fv_best_robustnessQ_Rvalues, mean)),
                                R_best_sd = unlist(lapply(fv_best_robustnessQ_Rvalues, sd)),
                                R_worst_mean = unlist(lapply(fv_worst_robustnessQ_Rvalues, mean)),
                                R_sd_mean = unlist(lapply(fv_worst_robustnessQ_Rvalues, sd)),
                                R_height_mean = unlist(lapply(fv_height_robustnessQ_Rvalues, mean)),
                                R_height_sd = unlist(lapply(fv_height_robustnessQ_Rvalues, sd)), 
                                R_sla_mean = unlist(lapply(fv_sla_robustnessQ_Rvalues, mean)),
                                R_sla_sd = unlist(lapply(fv_sla_robustnessQ_Rvalues, sd)))

# Labels
fv.labs <- c("FV_PO", "FV_ME", "FV_MA")
names(fv.labs) <- fv_Rvalues_dframe$network

# ggplot rendering 
plot2 <- ggplot() +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "black", linewidth = 1, inherit.aes = FALSE,
               data = fv_rand_robustnessQ_plot_long) + 
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "darkred", linewidth = 1, inherit.aes = FALSE,
               data = fv_worst_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "seagreen", linewidth = 1, inherit.aes = FALSE,
               data = fv_best_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "orange", linewidth = 1, inherit.aes = FALSE,
               data = fv_height_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "orangered3", linewidth = 1, inherit.aes = FALSE,
               data = fv_sla_robustnessQ_plot_long) +
  ylab("Proportion of consumers remaining") + 
  xlab("Proportion of resources removed") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 8, colour = "black"), 
        strip.text = element_text(face = "bold", hjust = 0, size = 10),
        strip.background = element_blank()) + 
  facet_wrap(~ network, nrow = 1)
plot2


## Clean directory (remove everything that we don't need later)
rm(mymat, myobsmatrix, resourcesurvivors, trial, chosen.order, consumer, i, j, 
   k, nperm, p, resource, sRvalues, survivors, survivorsx, threshold, trigger,
   triggerhat, v.ext, extlists, exttimes, res.names, sla_ext_list, pocock_sla,
   pocock_heights, pocock_fv_clean, pocock_fv, pocock, memmott1999, memmott,
   memmott_fv, plantat_traits, plantheights, plantsla, magrach_sub, 
   fv_edgelists, core_traits, memmott_heights, memmott_fv_clean, memmott_sla, 
   magrach_sla, magrach_heights, magrach_fv_clean, magrach_fv, magrach,
   height_ext_list, fv_networks, fv_matrices,
   fv_best_robustnessQ_plot_wide, fv_rand_robustnessQ_plot_wide,
   fv_worst_robustnessQ_plot_wide, fv_sla_robustnessQ_plot_wide, 
   fv_height_robustnessQ_plot_wide, fv_rand_robustnessQ_Rvalues, 
   fv_best_robustnessQ_Rvalues, fv_worst_robustnessQ_Rvalues, 
   fv_sla_robustnessQ_Rvalues, fv_height_robustnessQ_Rvalues, 
   fv_best_robustnessQ_curves, fv_rand_robustnessQ_curves, 
   fv_worst_robustnessQ_curves, fv_sla_robustnessQ_curves, 
   fv_height_robustnessQ_curves, fv.labs)


## Tetrapod meta-webs and longevity (age)

# Sort out the datasets
tetraweb <- read.csv("Data/Tetrapods/TetraEU_pairwise_interactions.csv", sep = ";")
age <- readr::read_delim("Data/Tetrapods/anage_data.txt")
age$name <- paste(age$Genus, age$Species, sep = " ") # create a merged species column

# Merge the two datasets together
TP_age <- left_join(tetraweb, age, by = c("targetTaxonName" = "name")) %>%
  dplyr::select(targetTaxonName, `Maximum longevity (yrs)`) %>%
  filter(!duplicated(targetTaxonName)) %>%
  filter(!is.na(`Maximum longevity (yrs)`))

# Subset the database by the resources that have longevity data (525 taxa)
tetraweb_age <- tetraweb[tetraweb$targetTaxonName %in% TP_age$targetTaxonName,]

# Randomly sample a number of links 
set.seed(123)

# Randomly select a a number of nodes between 30 and 100 (to keep it relatively)
# consistent with the size of the other empirical datasets (25 to 74 nodes)
rand <- round(runif(3, min = 20, max = 75)) # 40, 57 and 71

# Create a list of the resource nodes (those which will be the primary extinction)
nodes <- unique(tetraweb_age$targetTaxonName)

# Create the subset lists of nodes (i.e., sample N nodes from the network)
TP1_nodes <- sample(nodes, rand[1])
TP2_nodes <- sample(nodes, rand[2])
TP3_nodes <- sample(nodes, rand[3])

# Subset the networks to the number of links from above
TP_1 <- tetraweb_age[tetraweb_age$targetTaxonName %in% TP1_nodes,]
TP_2 <- tetraweb_age[tetraweb_age$targetTaxonName %in% TP2_nodes,]
TP_3 <- tetraweb_age[tetraweb_age$targetTaxonName %in% TP3_nodes,]

# Number of nodes 
length(unique(c(TP_1$sourceTaxonId, TP_1$targetTaxonId)))
length(unique(c(TP_2$sourceTaxonId, TP_2$targetTaxonId)))
length(unique(c(TP_3$sourceTaxonId, TP_3$targetTaxonId)))

# Clean the dataset to remove additional columns
TP_1_clean <- TP_1 %>%
  dplyr::select(targetTaxonName, sourceTaxonName) %>%
  rename(resource.name = targetTaxonName, consumer.name = sourceTaxonName)

# Add a weight column
TP_1_clean$weight <- rep(1, nrow(TP_1_clean))

# Create a list of resources and their longevity
TP_1_age <- left_join(TP_1_clean, age, by = c("resource.name" = "name")) %>%
  dplyr::select(resource.name, `Maximum longevity (yrs)`) %>%
  filter(!duplicated(resource.name))

# Clean the dataset to remove additional columns
TP_2_clean <- TP_2 %>%
  dplyr::select(targetTaxonName, sourceTaxonName) %>%
  rename(resource.name = targetTaxonName, consumer.name = sourceTaxonName)

# Add a weight column
TP_2_clean$weight <- rep(1, nrow(TP_2_clean))

# Create a list of resources and their longevity
TP_2_age <- left_join(TP_2_clean, age, by = c("resource.name" = "name")) %>%
  dplyr::select(resource.name, `Maximum longevity (yrs)`) %>%
  filter(!duplicated(resource.name))

# Clean the dataset to remove additional columns
TP_3_clean <- TP_3 %>%
  dplyr::select(targetTaxonName, sourceTaxonName) %>%
  rename(resource.name = targetTaxonName, consumer.name = sourceTaxonName)

# Add a weight column
TP_3_clean$weight <- rep(1, nrow(TP_3_clean))

# Create a list of resources and their longevity
TP_3_age <- left_join(TP_3_clean, age, by = c("resource.name" = "name")) %>%
  dplyr::select(resource.name, `Maximum longevity (yrs)`) %>%
  filter(!duplicated(resource.name))

## Create lists of data

# Organise the data into a consistent format with consistent column names

# Create a list of edgelists
tp_edgelists <- list(TP_1_clean, TP_2_clean, TP_3_clean)
names(tp_edgelists) <- c("TP_1", "TP_2", "TP_3")

# Create networks from the edgelists
tp_networks <- lapply(tp_edgelists, graph_from_data_frame)
names(tp_networks) <- c("TP_1", "TP_2", "TP_3")

# Create incidence matrices for the networks
tp_matrices <- lapply(tp_networks, function(x){bipartite::empty(as.matrix(as_adjacency_matrix(x)))})
names(tp_matrices) <- c("TP_1", "TP_2", "TP_3")

# Plant height scenario
ages <- list(TP_1_age, TP_2_age, TP_3_age)
names(ages) <- c("TP_1", "TP_2", "TP_3")

# Create the height ordered extinction lists (shorter plants are)
age_ext_list <- lapply(ages, function(x){x[base::order(x$`Maximum longevity (yrs)`, decreasing = T),]})


#### Robustness (from Bane et al. 2018) ####

# For each network extract the different extinction scenarios
# 1. Based on large to small body size
# 2. Based on most to least connected (out degree to create bottom-up losses)
# 3. Based on least to most connected (out degree to create bottom-up losses)
# 4. Mean from 1000 random extinction scenarios

## Random extinction sequence

# Set up the empty lists
tp_rand_robustnessQ_curves<-NULL; tp_rand_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(tp_matrices)){
  myobsmatrix <- t(tp_matrices[[k]])
  nperm<-1000    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix          # make copy of original matrix to work on
    resource<-colSums(mymat)    # save resource degrees
    consumer<-rowSums(mymat)    # save consumer degrees
    survivors<-NULL             # create survivors to save consumer counts
    triggerhat<-colnames(mymat) # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-sample(triggerhat,1) # pick a trigger from the hat
      mymat[,trigger]<-0            # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  tp_rand_robustnessQ_Rvalues[[k]]<-sRvalues
  tp_rand_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
tp_rand_robustnessQ_plot_wide<-NULL
for (p in 1:length(tp_rand_robustnessQ_curves)){
  trial <- tp_rand_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  tp_rand_robustnessQ_plot_wide[[p]]<-melt(trial)
}

# Rename the list to be the names of the channels 
names(tp_rand_robustnessQ_plot_wide) <- names(tp_networks)

# Bind the results together and label based on the channels
tp_rand_robustnessQ_plot_long <- bind_rows(tp_rand_robustnessQ_plot_wide, .id = "network")


## Worst case (high to low degree) extinction sequence

# Set up the empty lists
tp_worst_robustnessQ_curves<-NULL; tp_worst_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(tp_matrices)){
  myobsmatrix <- t(tp_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=TRUE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger 
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  tp_worst_robustnessQ_Rvalues[[k]]<-sRvalues
  tp_worst_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
tp_worst_robustnessQ_plot_wide<-NULL
for (p in 1:length(tp_worst_robustnessQ_curves)){
  trial <- tp_worst_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  tp_worst_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(tp_worst_robustnessQ_plot_wide) <- names(tp_networks)

# Bind the results together and label based on the channels
tp_worst_robustnessQ_plot_long <- bind_rows(tp_worst_robustnessQ_plot_wide, .id = "network")


## Best case (high to low degree) extinction sequence

# Set up the empty lists
tp_best_robustnessQ_curves<-NULL; tp_best_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(tp_matrices)){
  myobsmatrix <- t(tp_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=FALSE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1]   # pick a trigger from the hat
      mymat[,trigger]<-0      # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  tp_best_robustnessQ_Rvalues[[k]]<-sRvalues
  tp_best_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
tp_best_robustnessQ_plot_wide<-NULL
for (p in 1:length(tp_best_robustnessQ_curves)){
  trial <- tp_best_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  tp_best_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(tp_best_robustnessQ_plot_wide) <- names(tp_networks)

# Bind the results together and label based on the channels
tp_best_robustnessQ_plot_long <- bind_rows(tp_best_robustnessQ_plot_wide, .id = "network")


## Age dependent extinction sequence

# Set up the empty lists
tp_age_robustnessQ_curves<-NULL; tp_age_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(tp_matrices)){
  myobsmatrix <- t(tp_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                            # make copy of original matrix to work on
    resource<-colSums(mymat)                      # save resource degrees
    consumer<-rowSums(mymat)                      # save consumer degrees
    res.names<- colnames(mymat)                   # names of the resources
    survivors<-NULL                               # create survivors to save consumer counts
    chosen.order<-age_ext_list[[k]]$resource.name # sla based extinction lists                 
    triggerhat<-chosen.order                      # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger  
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  tp_age_robustnessQ_Rvalues[[k]]<-sRvalues
  tp_age_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
tp_age_robustnessQ_plot_wide<-NULL
for (p in 1:length(tp_age_robustnessQ_curves)){
  trial <- tp_age_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  tp_age_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(tp_age_robustnessQ_plot_wide) <- names(tp_networks)

# Bind the results together and label based on the channels
tp_age_robustnessQ_plot_long <- bind_rows(tp_age_robustnessQ_plot_wide, .id = "network")

## Wrangle the data to analyse later

# Sort the robustness R value results 
tp_Rvalues_dframe <- data.frame(network = names(tp_matrices),
                                R_rand_mean = unlist(lapply(tp_rand_robustnessQ_Rvalues, mean)),
                                R_rand_sd = unlist(lapply(tp_rand_robustnessQ_Rvalues, sd)),
                                R_best_mean = unlist(lapply(tp_best_robustnessQ_Rvalues, mean)),
                                R_best_sd = unlist(lapply(tp_best_robustnessQ_Rvalues, sd)),
                                R_worst_mean = unlist(lapply(tp_worst_robustnessQ_Rvalues, mean)),
                                R_sd_mean = unlist(lapply(tp_worst_robustnessQ_Rvalues, sd)),
                                R_age_mean = unlist(lapply(tp_age_robustnessQ_Rvalues, mean)),
                                R_age_sd = unlist(lapply(tp_age_robustnessQ_Rvalues, sd)))

# Labels
tp.labs <- c("TP_1", "TP_2", "TP_3")
names(tp.labs) <- tp_Rvalues_dframe$network

# ggplot rendering 
plot3 <- ggplot() +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "black", linewidth = 1, inherit.aes = FALSE,
               data = tp_rand_robustnessQ_plot_long) + 
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "darkred", linewidth = 1, inherit.aes = FALSE,
               data = tp_worst_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "seagreen", linewidth = 1, inherit.aes = FALSE,
               data = tp_best_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "orange", linewidth = 1, inherit.aes = FALSE,
               data = tp_age_robustnessQ_plot_long) +
  ylab("Proportion of consumers remaining") + 
  xlab("Proportion of resources removed") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 8, colour = "black"), 
        strip.text = element_text(face = "bold", hjust = 0, size = 10),
        strip.background = element_blank()) + 
  facet_wrap(~ network, nrow = 1)
plot3

## Clear environment
rm(age, age_ext_list, ages, chosen.order, consumer, extlists, exttimes, i, j, 
   k, mymat, myobsmatrix, nodes, nperm, p, rand, res.names, resource,
   resourcesurvivors, sRvalues, survivors, survivorsx, tetraweb,
   tetraweb_age, threshold, TP_1, TP_1_age, TP_1_clean ,TP_2,
   TP_2_age,TP_2_clean,TP_3, TP_3_age, TP_3_clean, TP_age,
   tp_age_robustnessQ_curves, tp_age_robustnessQ_plot_wide,
   tp_age_robustnessQ_Rvalues,  tp_best_robustnessQ_curves,
   tp_best_robustnessQ_plot_wide, tp_best_robustnessQ_Rvalues, tp_edgelists,
   tp_matrices, tp_networks,  tp_rand_robustnessQ_curves, 
   tp_rand_robustnessQ_plot_wide, tp_rand_robustnessQ_Rvalues, 
   tp_worst_robustnessQ_curves, tp_worst_robustnessQ_plot_wide, 
   tp_worst_robustnessQ_Rvalues, tp.labs, TP1_nodes, TP2_nodes, TP3_nodes,
   trial, trigger, triggerhat, v.ext)


## Pollution tolerance for stream macroinvertebrates

# Read in the edgelists
SW_T1 <- read.csv("Data/Streams/T1_edgelist.csv") %>% filter(!(resource == "CPOM" | resource == "FPOM" | consumer == "Cottus gobio"))
SW_T2 <- read.csv("Data/Streams/T2_edgelist.csv") %>% filter(!(resource == "CPOM" | resource == "FPOM" | consumer == "Cottus gobio"))
SW_W1 <- read.csv("Data/Streams/W1_edgelist.csv") %>% filter(!(resource == "CPOM" | resource == "FPOM" | consumer == "Cottus gobio"))

# Number of nodes
length(unique(c(SW_T1$resource, SW_T1$consumer))) # 42
length(unique(c(SW_T2$resource, SW_T2$consumer))) # 48
length(unique(c(SW_W1$resource, SW_W1$consumer))) # 52

# Abundances
SW_T1_abun <- read.csv("Data/Streams/T1_nodes.csv") %>% filter(!(node == "CPOM" | node == "FPOM" | node == "Cottus gobio"))
SW_T2_abun <- read.csv("Data/Streams/T2_nodes.csv") %>% filter(!(node == "CPOM" | node == "FPOM" | node == "Cottus gobio"))
SW_W1_abun <- read.csv("Data/Streams/W1_nodes.csv") %>% filter(!(node == "CPOM" | node == "FPOM" | node == "Cottus gobio"))

# Clean the dataset to remove additional columns
SW_T1_clean <- SW_T1 %>%
  dplyr::select(resource, consumer) %>%
  rename(resource.name = resource, consumer.name = consumer) 

# Add a weight column
SW_T1_clean$weight <- rep(1, nrow(SW_T1_clean))

# Create a list of resources and their longevity
SW_T1_pop <- left_join(SW_T1_clean, SW_T1_abun, by = c("resource.name" = "node")) %>%
  dplyr::select(resource.name, WHPT) %>%
  filter(!duplicated(resource.name))

# Clean the dataset to remove additional columns
SW_T2_clean <- SW_T2 %>%
  dplyr::select(resource, consumer) %>%
  rename(resource.name = resource, consumer.name = consumer)

# Add a weight column
SW_T2_clean$weight <- rep(1, nrow(SW_T2_clean))

# Create a list of resources and their longevity
SW_T2_pop <- left_join(SW_T2_clean, SW_T2_abun, by = c("resource.name" = "node")) %>%
  dplyr::select(resource.name, WHPT) %>%
  filter(!duplicated(resource.name))

# Clean the dataset to remove additional columns
SW_W1_clean <- SW_W1 %>%
  dplyr::select(resource, consumer) %>%
  rename(resource.name = resource, consumer.name = consumer)

# Add a weight column
SW_W1_clean$weight <- rep(1, nrow(SW_W1_clean))

# Create a list of resources and their longevity
SW_W1_pop <- left_join(SW_W1_clean, SW_W1_abun, by = c("resource.name" = "node")) %>%
  dplyr::select(resource.name, WHPT) %>%
  filter(!duplicated(resource.name))


## Create lists of data

# Organise the data into a consistent format with consistent column names

# Create a list of edgelists
sw_edgelists <- list(SW_T1_clean, SW_T2_clean, SW_W1_clean)
names(sw_edgelists) <- c("SW_T1", "SW_T2", "SW_W3")

# Create networks from the edgelists
sw_networks <- lapply(sw_edgelists, graph_from_data_frame)
names(sw_networks) <- c("SW_T1", "SW_T2", "SW_W3")

# Create incidence matrices for the networks
sw_matrices <- lapply(sw_networks, function(x){bipartite::empty(as.matrix(as_adjacency_matrix(x)))})
names(sw_matrices) <- c("SW_T1", "SW_T2", "SW_W3")

# Plant height scenario
pops <- list(SW_T1_pop, SW_T2_pop, SW_W1_pop)
names(pops) <- c("SW_T1", "SW_T2", "SW_W3")

# Create the height ordered extinction lists (shorter plants are)
pop_ext_list <- lapply(pops, function(x){x[base::order(x$WHPT, decreasing = T),]})


#### Robustness (from Bane et al. 2018) ####

# For each network extract the different extinction scenarios
# 1. Based on large to small body size
# 2. Based on most to least connected (out degree to create bottom-up losses)
# 3. Based on least to most connected (out degree to create bottom-up losses)
# 4. Mean from 1000 random extinction scenarios

## Random extinction sequence

# Set up the empty lists
sw_rand_robustnessQ_curves<-NULL; sw_rand_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(sw_matrices)){
  myobsmatrix <- t(sw_matrices[[k]])
  nperm<-1000    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix          # make copy of original matrix to work on
    resource<-colSums(mymat)    # save resource degrees
    consumer<-rowSums(mymat)    # save consumer degrees
    survivors<-NULL             # create survivors to save consumer counts
    triggerhat<-colnames(mymat) # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-sample(triggerhat,1) # pick a trigger from the hat
      mymat[,trigger]<-0            # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  sw_rand_robustnessQ_Rvalues[[k]]<-sRvalues
  sw_rand_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
sw_rand_robustnessQ_plot_wide<-NULL
for (p in 1:length(sw_rand_robustnessQ_curves)){
  trial <- sw_rand_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  sw_rand_robustnessQ_plot_wide[[p]]<-melt(trial)
}

# Rename the list to be the names of the channels 
names(sw_rand_robustnessQ_plot_wide) <- names(sw_networks)

# Bind the results together and label based on the channels
sw_rand_robustnessQ_plot_long <- bind_rows(sw_rand_robustnessQ_plot_wide, .id = "network")


## Worst case (high to low degree) extinction sequence

# Set up the empty lists
sw_worst_robustnessQ_curves<-NULL; sw_worst_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(sw_matrices)){
  myobsmatrix <- t(sw_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=TRUE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger 
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  sw_worst_robustnessQ_Rvalues[[k]]<-sRvalues
  sw_worst_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
sw_worst_robustnessQ_plot_wide<-NULL
for (p in 1:length(sw_worst_robustnessQ_curves)){
  trial <- sw_worst_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  sw_worst_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(sw_worst_robustnessQ_plot_wide) <- names(sw_networks)

# Bind the results together and label based on the channels
sw_worst_robustnessQ_plot_long <- bind_rows(sw_worst_robustnessQ_plot_wide, .id = "network")


## Best case (high to low degree) extinction sequence

# Set up the empty lists
sw_best_robustnessQ_curves<-NULL; sw_best_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(sw_matrices)){
  myobsmatrix <- t(sw_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                                        # make copy of original matrix to work on
    resource<-colSums(mymat)                                  # save resource degrees
    consumer<-rowSums(mymat)                                  # save consumer degrees
    survivors<-NULL                                           # create survivors to save consumer counts
    chosen.order<-names(sort(colSums(myobsmatrix), decreasing=FALSE))
    triggerhat<-chosen.order                                  # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1]   # pick a trigger from the hat
      mymat[,trigger]<-0      # make it extinct in the matrix 
      exttimes[j,i]<-trigger
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  sw_best_robustnessQ_Rvalues[[k]]<-sRvalues
  sw_best_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
sw_best_robustnessQ_plot_wide<-NULL
for (p in 1:length(sw_best_robustnessQ_curves)){
  trial <- sw_best_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  sw_best_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(sw_best_robustnessQ_plot_wide) <- names(sw_networks)

# Bind the results together and label based on the channels
sw_best_robustnessQ_plot_long <- bind_rows(sw_best_robustnessQ_plot_wide, .id = "network")


## Population size dependent extinction sequence

# Set up the empty lists
sw_pop_robustnessQ_curves<-NULL; sw_pop_robustnessQ_Rvalues<-NULL

# Run the robustness analysis
for (k in 1:length(sw_matrices)){
  myobsmatrix <- t(sw_matrices[[k]])
  nperm<-1    
  threshold<-1                                                                  
  sRvalues<-NULL                                                                
  resourcesurvivors<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  exttimes<-matrix(nrow=nperm, ncol=ncol(myobsmatrix))
  extlists<-matrix(ncol=nperm, nrow=(ncol(myobsmatrix)+1))
  
  for(j in 1:nperm){
    
    mymat<-myobsmatrix                            # make copy of original matrix to work on
    resource<-colSums(mymat)                      # save resource degrees
    consumer<-rowSums(mymat)                      # save consumer degrees
    res.names<- colnames(mymat)                   # names of the resources
    survivors<-NULL                               # create survivors to save consumer counts
    chosen.order<-pop_ext_list[[k]]$resource.name # population size based extinction lists                 
    triggerhat<-chosen.order                      # create hat with resource names to pick from
    
    for(i in 1:length(triggerhat)) {
      
      trigger<-triggerhat[1] # pick a trigger from the hat
      mymat[,trigger]<-0     # make it extinct in the matrix 
      exttimes[j,i]<-trigger  
      
      v.ext<-which(((consumer-rowSums(mymat))/consumer)>=threshold) # check consumers
      mymat[v.ext,]<-0                                              # make consumers extinct
      survivors[i]<-length(which(rowSums(mymat)>0))                 # save number of survivors
      triggerhat<-triggerhat[triggerhat!=trigger]                   # update trigger hat
      resourcesurvivors[j,i]<-length(which(colSums(mymat)==0))
    }
    
    survivorsx<-c(nrow(mymat),survivors)                   # add starting total to survivors list
    extlists[,j]<-survivorsx/nrow(mymat)                   # Save the proportion of remaining consumers
    sRvalues[j]<-sum(survivorsx)/(ncol(mymat)*nrow(mymat)) # Calculate the area under the curve
  }
  print(k)
  sw_pop_robustnessQ_Rvalues[[k]]<-sRvalues
  sw_pop_robustnessQ_curves[[k]]<-extlists
}

# Prepare the data for plotting 
sw_pop_robustnessQ_plot_wide<-NULL
for (p in 1:length(sw_pop_robustnessQ_curves)){
  trial <- sw_pop_robustnessQ_curves[[p]]
  rownames(trial) <- seq(0,(nrow(trial)-1),1)/(nrow(trial)-1)
  sw_pop_robustnessQ_plot_wide[[p]] <- melt(trial)
}

# Rename the list to be the names of the channels 
names(sw_pop_robustnessQ_plot_wide) <- names(sw_networks)

# Bind the results together and label based on the channels
sw_pop_robustnessQ_plot_long <- bind_rows(sw_pop_robustnessQ_plot_wide, .id = "network")

## Wrangle the data to analyse later

# Sort the robustness R value results 
sw_Rvalues_dframe <- data.frame(network = names(sw_matrices),
                                R_rand_mean = unlist(lapply(sw_rand_robustnessQ_Rvalues, mean)),
                                R_rand_sd = unlist(lapply(sw_rand_robustnessQ_Rvalues, sd)),
                                R_best_mean = unlist(lapply(sw_best_robustnessQ_Rvalues, mean)),
                                R_best_sd = unlist(lapply(sw_best_robustnessQ_Rvalues, sd)),
                                R_worst_mean = unlist(lapply(sw_worst_robustnessQ_Rvalues, mean)),
                                R_sd_mean = unlist(lapply(sw_worst_robustnessQ_Rvalues, sd)),
                                R_pop_mean = unlist(lapply(sw_pop_robustnessQ_Rvalues, mean)),
                                R_pop_sd = unlist(lapply(sw_pop_robustnessQ_Rvalues, sd)))

# Labels
sw.labs <- c("SW_T1", "SW_T2", "SW_W1")
names(sw.labs) <- sw_Rvalues_dframe$network

# ggplot rendering 
plot4 <- ggplot() +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "black", linewidth = 1, inherit.aes = FALSE,
               data = sw_rand_robustnessQ_plot_long) + 
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "darkred", linewidth = 1, inherit.aes = FALSE,
               data = sw_worst_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "seagreen", linewidth = 1, inherit.aes = FALSE,
               data = sw_best_robustnessQ_plot_long) +
  stat_summary(aes(x=Var1, y=value), geom = "line",
               fun = mean, colour = "orange", linewidth = 1, inherit.aes = FALSE,
               data = sw_pop_robustnessQ_plot_long) +
  ylab("Proportion of consumers remaining") + 
  xlab("Proportion of resources removed") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  theme_bw() + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 8, colour = "black"), 
        strip.text = element_text(hjust = 0, size = 10, face = "bold"),
        strip.background = element_blank()) + 
  facet_wrap(~ network, nrow = 1)
plot4


## Clear environment
rm(pops, pop_ext_list, chosen.order, consumer, extlists, exttimes, i, j, 
   k, mymat, myobsmatrix, nperm, p, res.names, resource,
   resourcesurvivors, sRvalues, survivors, survivorsx,
   threshold, SW_T1, SW_T1_pop, SW_T1_clean, SW_T2,
   SW_T2_pop, SW_T2_clean, SW_W1, SW_W1_pop, SW_W1_clean,
   sw_pop_robustnessQ_curves, sw_pop_robustnessQ_plot_wide,
   sw_pop_robustnessQ_Rvalues, sw_best_robustnessQ_curves,
   sw_best_robustnessQ_plot_wide, sw_best_robustnessQ_Rvalues, sw_edgelists,
   sw_matrices, sw_networks, sw_rand_robustnessQ_curves, 
   sw_rand_robustnessQ_plot_wide, sw_rand_robustnessQ_Rvalues, 
   sw_worst_robustnessQ_curves, sw_worst_robustnessQ_plot_wide, 
   sw_worst_robustnessQ_Rvalues, sw.labs, SW_T1_abun, SW_T2_abun, SW_W1_abun,
   trial, trigger, triggerhat, v.ext)


#### Figure 1 #### 

# Remove axes and create share ones
gridExtra::grid.arrange(plot1, plot3, plot4, plot2, nrow = 4, 
                        left = textGrob("Proportion of consumers remaining", gp=gpar(fontsize=10.5), rot = 90), 
                        bottom = textGrob("Proportion of resources removed", gp=gpar(fontsize=10.5)))


#### Figure 2 #### 

# An ugly bit of code to bind together dataframes
names(fw_Rvalues_dframe)[names(fw_Rvalues_dframe) == "foodweb"] <- "network"
test <- bind_rows(fw_Rvalues_dframe, tp_Rvalues_dframe, sw_Rvalues_dframe, fv_Rvalues_dframe)
test[1:3, "network"] <- c("FW_HS", "FW_LE", "FW_RA")

# Plot the Rvalues
plot2 <- ggplot() + 
  geom_linerange(aes(xmin=R_worst_mean, xmax = R_best_mean, y = network), data = test, colour = "grey60") + 
  geom_linerange(aes(xmin=R_rand_mean - R_rand_sd, xmax = R_rand_mean + R_rand_sd, y = network), data = test) + 
  geom_point(aes(x=R_rand_mean, y = network), data = test, size = 2) + 
  geom_point(aes(x=R_age_mean, y = network), data = test, colour = "orange", size = 2) + 
  geom_point(aes(x=R_pop_mean, y = network), data = test, colour = "orange", size = 2) + 
  geom_point(aes(x=R_sla_mean, y = network), data = test, colour = "orangered3", size = 2) + 
  geom_point(aes(x=R_height_mean, y = network), data = test, colour = "orange", size = 2) + 
  geom_point(aes(x=R_size_mean, y = network), data = test, colour = "orange", size = 2) + 
  theme_bw() + theme(axis.text = element_text(size = 8, colour = "black"), 
                     strip.text = element_text(hjust = 0, size = 10, face = "bold"),
                     strip.background = element_blank(), axis.title.y = element_blank()) + 
  coord_cartesian(xlim = c(0.2, 1.0)) + 
  xlab("R value") 
plot2

#### Statistical analysis #### 

# Reduce the random extinction scenarios - making them have the same steps as other extinction sequences
# X = Var and Y = Value

# Calculate the mean curve for all of the simulations for fv
fv_rand_mean_robustness_long <- aggregate(value ~ network + Var1, 
                                          FUN = mean, 
                                          data = fv_rand_robustnessQ_plot_long)

# Calculate the mean curve for all of the simulations for fw
fw_rand_mean_robustness_long <- aggregate(value ~ Foodweb + Var1, 
                                          FUN = mean, 
                                          data = fw_rand_robustnessQ_plot_long)

# Calculate the mean curve for all of the simulations for tp
tp_rand_mean_robustness_long <- aggregate(value ~ network + Var1, 
                                          FUN = mean, 
                                          data = tp_rand_robustnessQ_plot_long)

# Calculate the mean curve for all of the simulations for sw
sw_rand_mean_robustness_long <- aggregate(value ~ network + Var1, 
                                          FUN = mean, 
                                          data = sw_rand_robustnessQ_plot_long)

## Statistical tests ## 

# Comparisons for fw networks 

# Realistic versus best
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Harper-Smith et al., in press", "value"], 
        fw_best_robustnessQ_plot_long[fw_best_robustnessQ_plot_long$Foodweb == "Harper-Smith et al., in press", "value"])
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Ledger, Edwards & Woodward unpublished", "value"], 
        fw_best_robustnessQ_plot_long[fw_best_robustnessQ_plot_long$Foodweb == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Rayner, unpublished", "value"], 
        fw_best_robustnessQ_plot_long[fw_best_robustnessQ_plot_long$Foodweb == "Rayner, unpublished", "value"])

# Realistic versus worst
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Harper-Smith et al., in press", "value"], 
        fw_worst_robustnessQ_plot_long[fw_worst_robustnessQ_plot_long$Foodweb == "Harper-Smith et al., in press", "value"])
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Ledger, Edwards & Woodward unpublished", "value"], 
        fw_worst_robustnessQ_plot_long[fw_worst_robustnessQ_plot_long$Foodweb == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Rayner, unpublished", "value"], 
        fw_worst_robustnessQ_plot_long[fw_worst_robustnessQ_plot_long$Foodweb == "Rayner, unpublished", "value"])

# Realistic versus random
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Harper-Smith et al., in press", "value"], 
        fw_rand_mean_robustness_long[fw_rand_mean_robustness_long$Foodweb == "Harper-Smith et al., in press", "value"])
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Ledger, Edwards & Woodward unpublished", "value"], 
        fw_rand_mean_robustness_long[fw_rand_mean_robustness_long$Foodweb == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(fw_size_robustnessQ_plot_long[fw_size_robustnessQ_plot_long$Foodweb == "Rayner, unpublished", "value"], 
        fw_rand_mean_robustness_long[fw_rand_mean_robustness_long$Foodweb == "Rayner, unpublished", "value"])

# Comparisons for tp networks 

# Realistic versus best
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"], 
        tp_best_robustnessQ_plot_long[tp_best_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"])
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"], 
        tp_best_robustnessQ_plot_long[tp_best_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Rayner, unpublished", "value"], 
        tp_best_robustnessQ_plot_long[tp_best_robustnessQ_plot_long$network == "Rayner, unpublished", "value"])

# Realistic versus worst
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"], 
        fw_worst_robustnessQ_plot_long[tp_worst_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"])
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"], 
        tp_worst_robustnessQ_plot_long[tp_worst_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Rayner, unpublished", "value"], 
        tp_worst_robustnessQ_plot_long[tp_worst_robustnessQ_plot_long$network == "Rayner, unpublished", "value"])

# Realistic versus random
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"], 
        tp_rand_mean_robustness_long[tp_rand_mean_robustness_long$network == "Harper-Smith et al., in press", "value"])
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"], 
        tp_rand_mean_robustness_long[tp_rand_mean_robustness_long$network == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(tp_age_robustnessQ_plot_long[tp_age_robustnessQ_plot_long$network == "Rayner, unpublished", "value"], 
        tp_rand_mean_robustness_long[tp_rand_mean_robustness_long$network == "Rayner, unpublished", "value"])

# Comparisons for fv networks 

# Height

# Realistic versus best
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_PO", "value"], 
        fv_best_robustnessQ_plot_long[fv_best_robustnessQ_plot_long$network == "FV_PO", "value"])
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_ME", "value"], 
        fv_best_robustnessQ_plot_long[fv_best_robustnessQ_plot_long$network == "FV_ME", "value"])
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_MA", "value"], 
        fv_best_robustnessQ_plot_long[fv_best_robustnessQ_plot_long$network == "FV_MA", "value"])

# Realistic versus worst
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_PO", "value"], 
        fv_worst_robustnessQ_plot_long[fv_worst_robustnessQ_plot_long$network == "FV_PO", "value"])
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_ME", "value"], 
        fv_worst_robustnessQ_plot_long[fv_worst_robustnessQ_plot_long$network == "FV_ME", "value"])
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_MA", "value"], 
        fv_worst_robustnessQ_plot_long[fv_worst_robustnessQ_plot_long$network == "FV_MA", "value"])

# Realistic versus random
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_PO", "value"], 
        fv_rand_mean_robustness_long[fv_rand_mean_robustness_long$network == "FV_PO", "value"])
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_ME", "value"], 
        fv_rand_mean_robustness_long[fv_rand_mean_robustness_long$network == "FV_ME", "value"])
ks.test(fv_height_robustnessQ_plot_long[fv_height_robustnessQ_plot_long$network == "FV_MA", "value"], 
        fv_rand_mean_robustness_long[fv_rand_mean_robustness_long$network == "FV_MA", "value"])

# SLA

# Realistic versus best
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_PO", "value"], 
        fv_best_robustnessQ_plot_long[fv_best_robustnessQ_plot_long$network == "FV_PO", "value"])
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_ME", "value"], 
        fv_best_robustnessQ_plot_long[fv_best_robustnessQ_plot_long$network == "FV_ME", "value"])
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_MA", "value"], 
        fv_best_robustnessQ_plot_long[fv_best_robustnessQ_plot_long$network == "FV_MA", "value"])

# Realistic versus worst
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_PO", "value"], 
        fv_worst_robustnessQ_plot_long[fv_worst_robustnessQ_plot_long$network == "FV_PO", "value"])
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_ME", "value"], 
        fv_worst_robustnessQ_plot_long[fv_worst_robustnessQ_plot_long$network == "FV_ME", "value"])
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_MA", "value"], 
        fv_worst_robustnessQ_plot_long[fv_worst_robustnessQ_plot_long$network == "FV_MA", "value"])

# Realistic versus random
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_PO", "value"], 
        fv_rand_mean_robustness_long[fv_rand_mean_robustness_long$network == "FV_PO", "value"])
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_ME", "value"], 
        fv_rand_mean_robustness_long[fv_rand_mean_robustness_long$network == "FV_ME", "value"])
ks.test(fv_sla_robustnessQ_plot_long[fv_sla_robustnessQ_plot_long$network == "FV_MA", "value"], 
        fv_rand_mean_robustness_long[fv_rand_mean_robustness_long$network == "FV_MA", "value"])


# Comparisons for sw networks 

# Realistic versus best
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"], 
        sw_best_robustnessQ_plot_long[sw_best_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"])
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"], 
        sw_best_robustnessQ_plot_long[sw_best_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Rayner, unpublished", "value"], 
        sw_best_robustnessQ_plot_long[sw_best_robustnessQ_plot_long$network == "Rayner, unpublished", "value"])

# Realistic versus worst
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"], 
        fw_worst_robustnessQ_plot_long[sw_worst_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"])
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"], 
        sw_worst_robustnessQ_plot_long[sw_worst_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Rayner, unpublished", "value"], 
        sw_worst_robustnessQ_plot_long[sw_worst_robustnessQ_plot_long$network == "Rayner, unpublished", "value"])

# Realistic versus random
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Harper-Smith et al., in press", "value"], 
        sw_rand_mean_robustness_long[sw_rand_mean_robustness_long$network == "Harper-Smith et al., in press", "value"])
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Ledger, Edwards & Woodward unpublished", "value"], 
        sw_rand_mean_robustness_long[sw_rand_mean_robustness_long$network == "Ledger, Edwards & Woodward unpublished", "value"])
ks.test(sw_age_robustnessQ_plot_long[sw_age_robustnessQ_plot_long$network == "Rayner, unpublished", "value"], 
        sw_rand_mean_robustness_long[sw_rand_mean_robustness_long$network == "Rayner, unpublished", "value"])



#### Data sources #### 

# PLANTATT - https://www.brc.ac.uk/biblio/plantatt-attributes-british-and-irish-plants-spreadsheet
# Magrach - https://datadryad.org/stash/dataset/doi:10.5061/dryad.k0q1n
# Memmott - https://onlinelibrary.wiley.com/doi/full/10.1046/j.1461-0248.1999.00087.x
# Pocock - https://www.science.org/doi/full/10.1126/science.1214915 
# Brose - https://esajournals.onlinelibrary.wiley.com/doi/10.1890/05-0379
# Core - https://search.dataone.org/view/https%3A%2F%2Fpasta.lternet.edu%2Fpackage%2Fmetadata%2Feml%2Fedi%2F1533%2F2
# Tetra EU - https://onlinelibrary.wiley.com/doi/10.1111/geb.13138 
# AnAge - https://genomics.senescence.info/species/index.html 
# Windsor - https://www.sciencedirect.com/science/article/pii/S0043135419306244 
