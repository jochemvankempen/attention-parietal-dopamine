# Run the function run_twolines.R multiple times with different parameters

# load the relevant libraries
x<-c("tidyverse", "dplyr", "filenamer", "here")
lapply(x, function(x) {if (!require(x, character.only=T)) {install.packages(x);require(x)}})

# clear all
rm(list=ls(all=TRUE)) 

# set working directory
setwd('<YOUR_LOCAL_PATH>/Thiele-attention-gratc-LIP-pipette/')

# load functions ---------------------------------------------------
source(str_c(getwd(),'/repositories/attention-parietal-dopamine/analyses/external/','twolines.R'))
source(str_c(getwd(),'/repositories/attention-parietal-dopamine/analyses/run_twolines.R'))

# set file directory, load/save files here
filedir = str_c(getwd(),'/data/population/')

# parameter definition, loop over filenames, drugnames and selectivity indices
filenames = c('doseResponse_drugMI_','doseResponse_attAUROC_')
drugnames = c('Dopamine','SCH23390') #
selectivities = c('_none','_att&dru','_att','_dru')

# loop through filenames, drugnames and selectivities
for (filename in filenames) {
  for (drugname in drugnames) {
    for (selectivity in selectivities) {
    
      
      run_twolines(filedir, filename, drugname, selectivity)
    }
  }
}

