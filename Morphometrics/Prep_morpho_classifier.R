library(tidyverse)
library(tools)

### The following functions create classifier CSVs to be used for geometric morphometric
### analyses using the geomorph R package and Stereomorph R package landmark files

### SET WORKING DIRECTORY: Morphometrics folder

########################################
#                                      #
#           ## FUNCTIONS ##            #                              
#                                      #
########################################

create_classifier <- function(view = "Dorsal", specimen_file = "Philippine_muridae_fmnh.csv"){
  ### Create a classifier for each view of the skull
  # Make a table with the names of the shapes files 
  # (which specimen records do we need for the classifier csv?)
  directory <- paste(view, '_shapes', sep = "")
  df <- data.frame(matrix(file_path_sans_ext(list.files(directory)), nrow = length(dir(directory)), ncol = 1, 
                          dimnames = list(1:length(dir(directory)), 'ID')),
                   stringsAsFactors = FALSE)
  # Extract relevant columns and specimen records and make a csv that 
  # matches the order of the shape files (necessary for geomorph analysis)
  classifier <- read_csv(specimen_file) %>% 
    rename(ID = CatMammalsCatalogNo, 
           Species = IdeFiledAsQualifiedName, 
           Sex = DarSex,
           StateProvince = DarStateProvince,
           Elevation = DarMinimumElevationInMeters,
           County = DarCounty) %>%
    select(ID, Species, StateProvince, Elevation, Sex, County) 
  view_initial <- substr(view, start = 1, stop = 1)
  classifier$ID <- paste(classifier$ID, view_initial, sep = "")
  classifier <- inner_join(df, classifier)
  write.csv(classifier, paste(view, "_classifier.csv", sep = ""))
}

### Add species pairs for ANOVAS (see morphometrics_stats.R)
  # (!) ONLY USE this function for Chrotomys morphometrics project
  # (!) will overwrite existing classifier csv
add_Chrotomys_pair_columns <- function(view = "Dorsal"){
  classifier <- read.csv(paste(view, '_classifier.csv', sep = ""))
  CM = "Chrotomys mindorensis"
  S = "Chrotomys silaceus"
  W = "Chrotomys whiteheadi"
  L = "Luzon I"
  M = "Mindoro I"
  sp <- classifier$Species
  isl <- classifier$StateProvince
  # Full species
  classifier$CMS <- ifelse(sp %in% c(CM, S), "yes",
                           "no")
  classifier$CMW <- ifelse(sp %in% c(CM, W), "yes",
                           "no")
  classifier$SW <- ifelse(sp %in% c(S, W), "yes",
                           "no")
  # C. mindorensis populations
  classifier$MS <- ifelse(sp == S, "yes",
                          ifelse(sp == CM & isl == M, "yes",
                            "no"))
  classifier$MW <- ifelse(sp == W, "yes",
                          ifelse(sp == CM & isl == M, "yes",
                                 "no"))
  classifier$LW <- ifelse(sp == W, "yes",
                          ifelse(sp == CM & isl == L, "yes",
                                 "no"))
  classifier$LS <- ifelse(sp == S, "yes",
                          ifelse(sp == CM & isl == L, "yes",
                                 "no"))
  classifier$ML <- ifelse(sp == CM, "yes",
                        "no")
  # overwrite classifier csv
  write.csv(classifier, paste(view, '_classifier.csv', sep = ""))
}

########################################
#                                      #
#  ## APPLY FUNCTIONS TO ALL VIEWS ##  #                              
#                                      #
########################################

views <- c("Dorsal", "Lateral", "Ventral")
for (i in views){
  create_classifier(view = i)
}
### FOR CHROTOMYS ONLY
for (i in views){
  add_Chrotomys_pair_columns(view = i)
}

