#### Source code for the analysis of food web robustness based on node ... #### 
#### Code developed by Fredric Windsor #### 
#### Please contact fmwindsor@gmail.com for any questions ####

#### Functions #### 

## Data management

# Extracting an incidence matrix from an mgNetwork
extract_im <- function(mgCollect){ 
  
  # create an edgelist
  el <- mgCollect$interactions
  
  # get the labels for rows and columns
  rows <- unique(as.character(el$node_from))
  cols <- unique(as.character(el$node_to))
  
  # create a blank matrix to store the data
  network <- matrix(0, nrow = length(rows), ncol = length(cols), 
                    dimnames = list(as.character(rows), as.character(cols)))
  
  # extract the link weights from the edgelist and assign them to the matrix
  for (r in 1:nrow(el)){
    lowerTaxon <- as.character(el$node_from[r])
    upperTaxon <- as.character(el$node_to[r])
    network[lowerTaxon, upperTaxon] <- el$value[r]
 
  }
  
  return(network)
  
}

# Extracting the node names from dataframes (with id table for later conversion)
extract_node_names <- function(mgCollect){ 
  
  # prevent messages (annoying with the left_join function)
  #suppressMessages()
  
  # create objects for nodes and links
  nodes <- mgCollect$nodes
  interactions <- mgCollect$interactions
  
  # create an id dataframe for later
  node_id <- dplyr::select(nodes, node_id, original_name)
  
  # separate out the nodes into their bipartite levels
  node_lower <- data.frame(node_id = unique(interactions$node_to))
  node_upper <- data.frame(node_id = unique(interactions$node_from))
  
  # combine the nodes into a long dataframe
  nodes_all_d <- bind_rows(node_lower, node_upper, .id = "level")
  nodes_all <- distinct(nodes_all_d)
  
  # create a final dataset with ids, orginal names and levels
  node_info <- left_join(nodes_all, node_id)
  colnames(node_info) <- c("level", "id", "name")
  
  return(node_info)
  
}

# Merge trait or other data with interactions

merge_inter_traits <- function(mgCollect, upper_traits = NULL, lowertraits = NULL){
  
  edgelist <- mgCollect$interactions
  int_traits <- data.frame()
  
  ## Trait data need to be in long format
  
  if(is.null(lowertraits) == F & is.null(uppertraits) == T){
    lower_trait_data <- dplyr::select(lower_traits, AccSpeciesID, TraitNum)
    int_traits <- left_join(edgelist, trait_data, by = c("node_to" = "id"))
  return(lower_int_traits)
  }
  
  if(is.null(uppertraits) == F & is.null(lowertraits) == T){
    upper_trait_data <- dplyr::select(upper_traits, AccSpeciesID, TraitNum)
    int_traits <- left_join(edgelist, trait_data, by = c("node_from" = "id"))
  }
  
  if(is.null(uppertraits) == F & is.null(lowertraits) == F){
    lower_trait_data <- dplyr::select(lower_traits, AccSpeciesID, TraitNum)
    upper_trait_data <- dplyr::select(upper_traits, AccSpeciesID, TraitNum) # Need to change this to variables for another trait database
    
    lower_int_traits <- left_join(edgelist, trait_data, by = c("node_to" = "id"))
    upper_int_traits <- left_join(edgelist, trait_data, by = c("node_from" = "id"))
    int_traits <- bind_rows(lower_int_traits, upper_int_traits)
  }
}
