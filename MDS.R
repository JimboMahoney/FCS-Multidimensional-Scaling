## For comapring the similarity of multi-dim datasets using Multi-dimensional scaling
## Opens all FCS files in a directory (or asks user) and plots each file on a 2D plot
## More similar files = plotted closer together
## Needs at least 3 files



#########################################################
### Installing and loading required packages
#########################################################

if (!require("svDialogs")) {
  install.packages("svDialogs", dependencies = TRUE)
  library(svDialogs)
}

if (!require("tcltk2")) {
  install.packages("tcltk2", dependencies = TRUE)
  library(tcltk2)
}

if (!require("flowCore")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("flowCore")
    library(flowCore)
}

if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

if (!require("limma")) {
  BiocManager::install("limma")
  library(limma)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("ggrepel")) {
  install.packages("ggrepel", dependencies = TRUE)
  library(ggrepel)
}

# Clear environment
rm(list = ls(all = TRUE))


# library(svDialogs)# Moved to top
# Get user input for file
testfile<-dlg_open()
# Convert to string value
testfile <- capture.output(testfile)[7]

if ((testfile)=="character(0)"){
  stop("File input cancelled")
}else{
  
  #Remove invalid characters from file input location
  testfile <- gsub("[\"]","",testfile)
  testfile<-substring (testfile,5)
  
  #Set file and directory
  filename <- basename (testfile)
  dir <- dirname (testfile)
  
  # Set working directory accoding to file chosen
  setwd(dir)
  
  # List FCS (and fcs!) files in this location
  filesToOpenList <- unique(c(list.files(dir,pattern = ".FCS"),list.files(dir,pattern = ".fcs")))
  
  # Ask which files to open
  filesToOpen <- tk_select.list(filesToOpenList, multiple=TRUE,title="Select files to use.") 
  
  # Load all into a flowset
  #library(flowCore) # Moved to top
  fcs_raw <- read.flowSet(filesToOpen, transformation = FALSE,
                          truncate_max_range = FALSE)
  
  # Create summary table
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  
  # Create list of filenames for the fcs raw data
  filenamelist <- rep(filesToOpen, fsApply(fcs_raw, nrow))
  
  # Extract expression
  expr <- fsApply(fcs_raw, exprs)
  
  
  # Put into dataframe
  ggdf <- data.frame(file_name = filenamelist, expr)
  
  # Preparation work for arcsinh transform (columns is also used later for naming changes)
  # Create list of parameters
  columns<-colnames(ggdf)
  # Remove "Time" column to avoid it being transformed
  columns<-setdiff(columns,"Time")
  # Remove "Cell_Length" and Gaussians column to avoid it being transformed
  columns<-setdiff(columns,"Event_length")
  columns<-setdiff(columns,"Cell_length")
  columns<-setdiff(columns,"Center")
  columns<-setdiff(columns,"Offset")
  columns<-setdiff(columns,"Width")
  columns<-setdiff(columns,"Residual")
  columns<-setdiff(columns,"file_name")
  ## Remove FSC and SSC
  removefscssc<-grep("FSC|SSC",columns,value=TRUE)
  columns<-columns[! columns %in% removefscssc]
  
  
  #Remove unnecessary parameter text
  names(ggdf)[-1] <- sub("Di", "", names(ggdf)[-1])
  names(ggdf)[-1] <- sub("Dd", "", names(ggdf)[-1])
  # Create list of channel / parameter descriptions 
  params<-panel_fcs$desc
  # Add one NA at the front because ggdf also has filename
  params<-c(NA,params)
  # Replace parameters with descriptions, keeping things like Time, Event Length unchanged
  colnames(ggdf)[!is.na(params)] <- na.omit(params)
  ## Remove Time, Event_Length & Gaussian Parameters
  removecolumns <- c("Event_length", "Center", "Offset", "Width", "Residual", "Cell_length")
  ggdf <- ggdf[,!(names(ggdf) %in% removecolumns)]
  
  # Determine whether data is CyTOF or Flow by presence of FSC
  # isflow will be 0 for a CyTOF or greater than 1 if flow
  isflow <-sum(grep("FSC",colnames(ggdf)))
  
  ## Remove FSC and SSC
  # library(tidyverse) # Moved to top
  #ggdf <- ggdf %>% select(-contains("FSC"))
  #ggdf <- ggdf %>% select(-contains("SSC"))
  
    # Remove Time Column
  ggdf<-ggdf[,-2]
  
  
  # Ask user which markers to plot
  markerlist <- tk_select.list(colnames(ggdf), multiple=TRUE,title="Select Markers to use (Cancel to use all).") 
  
  # If user cancels dialog box, use all markers.
  if(length(markerlist)==0 ){
    markerlist <- colnames(ggdf)
    # Remove filename
    markerlist <- markerlist[-grepl("file_name",markerlist)]
  }
  
  # Crop data to that
  ggdf <- ggdf[,c("file_name",markerlist)]
  
  
  
  # apply arcsinh transformation
  # (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
  # cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)
  
  asinh_scale <- 5
  if (isflow>0){
    asinh_scale <- 150
  }
  ggdf[,-1] <- asinh(ggdf[,-1] / asinh_scale)
  
  # Collapse to only medians, ignoring first column (file_name)
  ggdf<-aggregate(ggdf[-1],ggdf[1],median)

# Transpose rows/columns (needed for MDS)
ggdf<-t(ggdf[,-1])  


# Add numbers to filenames
for (i in 1:length(filesToOpen)){
  filesToOpen[i] <-  paste(i,filesToOpen[i])
}

labels <- c(1:length(filesToOpen))

# Add column names 
colnames(ggdf)<-filesToOpen



# Calculate MDS (Multi-dim scaling)
#library(limma)# Moved to top
mds <- limma::plotMDS(ggdf, plot = FALSE)

# Create data frame
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   filename = colnames(ggdf))

# Plot
#library(ggrepel) # Moved to top
ggplot(ggdf, aes(x = MDS1, y = MDS2, color = filename)) +
  geom_point(size = 2, alpha = 0.8) +
  #show.legend = FALSE stops the legend points appearing as "a"
  geom_label_repel(aes(label = labels),size=3,force=10,hjust=1,direction="y",show.legend = FALSE) +
  theme_bw() +
  coord_fixed()#+
  # theme(legend.position = "none")


} # End of file cancel loop