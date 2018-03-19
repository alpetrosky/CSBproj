library(StereoMorph)
library(geomorph)
library(Morpho)
library(tidyverse)

### Run these statistical tests after landmarks have been digitized using Stereomorph 
  # make sure you've already set up your classifier csv's with Prep_morpho_classifier.R

### SET WORKING DIRECTORY: Morphometrics folder

### view refers to the aspect of the skull that was photographed:
  # Dorsal, Ventral, or Lateral 

# create figure and summary folders if they don't exist already
dir.create("Figures")
dir.create("Summaries")

########################################
#                                      #
#           ## FUNCTIONS ##            #                              
#                                      #
########################################

### MUST run set_landmarks_and_classifier before any of the statistical analyses

set_landmarks_and_classifier <- function(view = "Dorsal"){
  ### read in all your files
  shapes_file <- paste(view, '_shapes', sep = "")
  tps = paste(view, '.tps', sep = "")
  allfactors <- paste(view, '_classifier.csv', sep = "")
    # make sure Stereomorph coordinates will be accepted by geomorph
  shapes <- readShapes(shapes_file)
  array <<- shapes$landmarks.scaled
  writeland.tps(array, tps)
  landmarks <- readland.tps(tps, specID = "ID")
  gpa <<- gpagen(landmarks)
    # create factors
  classifier <<- read.csv(allfactors, header=T)
  species <<- factor(classifier$Species)
  island <<- factor(classifier$StateProvince)
  elevation <<- factor(classifier$MinimumElevation)
    # create figure and summary folders if they don't exist already
  #dir.create("Figures", "Summaries")
}

#########################################################################################

### Unfortunately most of the output classes from geomorph can't be coerced
  # into dataframes (hence sink rather than write.csv), so the CSV files are
  # ugly af but at least the values get exported

### Plot Generalized Procrustes Alignment
  # label = TRUE to number landmarks
run_GPA <- function(view = "Dorsal", label = FALSE){
  gpa_figure <- file.path("Figures", paste("Procrustes_alignment_", view, ".pdf", sep = ""))
  pdf(file = gpa_figure)
  plot(gpa, label = label)
  dev.off()
}

### PCA colored by species. Add: label = TRUE to find IDs of any outliers
  # re-run for PCs explaining >5% variance
  # pc.variance can't be coerced into a dataframe, write.csv won't work :'(
run_PCA <- function(view = "Dorsal", axis1 = 1, axis2= 2, label = FALSE){
  pca_figure <- file.path("Figures", 
                          paste(view, "_PC", axis1,"_PC", axis2, ".pdf", sep = ""))
  pdf(file = pca_figure)
  pca_csv <- file.path("Summaries", paste(view, "PCA_", view, ".csv", sep = ""))
  sink(pca_csv)
  pca <<- plotTangentSpace(gpa$coords, axis1 = axis1, axis2 = axis2,
                          groups = species, label = label)
  print(pca$pc.summary)
  dev.off()
  sink()
}


### Morphological disparity
run_morph_disparity <- function(view = "Dorsal"){
    # Chrotomys species
  disparity <- file.path("Summaries", paste("MorphDisparity_", view, ".csv", sep = ""))
  sink(disparity)
  morph <<- morphol.disparity(gpa$coords ~ species, print.progress = FALSE)
  
    # C. mindorensis populations
  subset.species <- coords.subset(A = gpa$coords, group = species)
  mind_only <- subset.species$`Chrotomys mindorensis`
  ymind <<- two.d.array(mind_only)
  MindVsLuzon <<- classifier$StateProvince[classifier$Species == 'Chrotomys mindorensis']
  summary(morphol.disparity(ymind ~ MindVsLuzon, print.progress = FALSE))
  sink()
}

### CVA (package: Morpho)
run_CVA <- function(view = "Dorsal"){
  cva <- file.path("Summaries", paste("CVA_", view, ".csv", sep = ""))
  sink(cva)
  print(CVA(array, species))
  print(CVA(ymind, MindVsLuzon))
  sink()
}

### Procrustes ANOVA, all species 
  # (package: geomorph)
run_allspeciesANOVA <- function(view = "Dorsal"){
  anova <- procD.lm(gpa$coords ~ species, print.progress = FALSE)
  anova1 <- as.data.frame(anova$aov.table)
  filename = paste("Summaries/ANOVA_", view, ".csv", sep = "")
  write.csv(anova1, file = filename)
}

### pairwise ANOVAs for full species and C. mindorensis populations
  # I tried to loop this, but I couldn't figure out how to get around 
  # the limitation to a single factor of the coords.subset function from geomorph 
run_Crotomys_ANOVAs <- function(view = "Dorsal"){
  CMSfactor <- factor(classifier$CMS) # C. mindorensis and C. silaceus
  CMScoords <<- coords.subset(A = gpa$coords, group = CMSfactor)
  CMSspecies <<- classifier$Species[classifier$CMS == 'yes']
  CMS <<- procD.lm(CMScoords$yes ~ CMSspecies, print.progress = FALSE)
  CMS <<- as.data.frame(CMS$aov.table)
  CMWfactor <- factor(classifier$CMW) # C. mindorensis and C. whiteheadi
  CMWcoords <<- coords.subset(A = gpa$coords, group = CMWfactor)
  CMWspecies <<- classifier$Species[classifier$CMW == 'yes']
  CMW <<- procD.lm(CMWcoords$yes ~ CMWspecies, print.progress = FALSE)
  CMW <<- as.data.frame(CMW$aov.table)
  SWfactor <- factor(classifier$SW) # C. silaceus and C. whiteheadi
  SWcoords <<- coords.subset(A = gpa$coords, group = SWfactor)
  SWspecies <<- classifier$Species[classifier$SW == 'yes']
  SW <<- procD.lm(SWcoords$yes ~ SWspecies, print.progress = FALSE)
  SW <<- as.data.frame(SW$aov.table)
  MSfactor <- factor(classifier$MS) # C. mind FROM MINDORO and C. silaceus
  MScoords <<- coords.subset(A = gpa$coords, group = MSfactor)
  MSspecies <<- classifier$Species[classifier$MS == 'yes']
  MS <<- procD.lm(MScoords$yes ~ MSspecies, print.progress = FALSE)
  MS <<- as.data.frame(MS$aov.table)
  MWfactor <<- factor(classifier$MW) # C. mind FROM MINDORO and C. whiteheadi
  MWcoords <<- coords.subset(A = gpa$coords, group = MWfactor)
  MWspecies <<- classifier$Species[classifier$MW == 'yes']
  MW <<- procD.lm(MWcoords$yes ~ MWspecies, print.progress = FALSE)
  MW <<- as.data.frame(MW$aov.table)
  LSfactor <- factor(classifier$LS) # C. mind FROM LUZON and C. silaceus
  LScoords <<- coords.subset(A = gpa$coords, group = LSfactor)
  LSspecies <<- classifier$Species[classifier$LS == 'yes']
  LS <<- procD.lm(LScoords$yes ~ LSspecies, print.progress = FALSE)
  LS <<- as.data.frame(LS$aov.table)
  LWfactor <- factor(classifier$LW) # C. mind FROM LUZON and C. whiteheadi
  LWcoords <<- coords.subset(A = gpa$coords, group = LWfactor)
  LWspecies <<- classifier$Species[classifier$LW == 'yes']
  LW <<- procD.lm(LWcoords$yes ~ LWspecies, print.progress = FALSE)
  LW <<- as.data.frame(LW$aov.table)
  MLfactor <- factor(classifier$ML) # C. mind FROM MINDORO and C. mind FROM LUZON
  MLcoords <<- coords.subset(A = gpa$coords, group = MLfactor)
  MLspecies <<- classifier$StateProvince[classifier$ML == 'yes']
  ML <<- procD.lm(MLcoords$yes ~ MLspecies, print.progress = FALSE)
  ML <<- as.data.frame(ML$aov.table)
  #anovas <- file.path("summaries", paste("Chrotomys_ANOVAs", view, ".csv", sep = ""))
  #sink(anovas)
  ChrotomysANOVA <<- rbind(CMS, CMW, SW, MS, MW, LS, LW, ML)
  filename = paste("Summaries/ANOVAsChrotomys_", view, ".csv", sep = "")
  write.csv(ChrotomysANOVA, file = filename)
  #sink()
}

########################################
#                                      #
#  ## APPLY FUNCTIONS TO ALL VIEWS ##  #                              
#                                      #
########################################

### You'll see lots of the message:
  # Warning: no geomorph data frame provided.
        # If an error occurs, this might be the reason.
# Individual values were used instead of stitching them together 
# to make data frames; it'll chill

# (!) includes pairs ANOVAs for Chrotomys

views <- c("Dorsal", "Lateral", "Ventral")
for (i in views){
  set_landmarks_and_classifier(view = i)
  run_GPA(view = i)
  run_PCA(view = i)
  run_CVA(view= i)
  run_allspeciesANOVA(view = i)
  run_Crotomys_ANOVAs(view = i)
}

# remember to look at the PCA results and re-run run_PCA for 
# all the PCs of interest (explaining >5% variance)
