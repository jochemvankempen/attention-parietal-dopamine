# load the relevant libraries
x<-c("tidyverse", "dplyr", "filenamer", "here")
lapply(x, function(x) {if (!require(x, character.only=T)) {install.packages(x);require(x)}})

# clear all
rm(list=ls(all=TRUE)) 

# load useful functions ---------------------------------------------------
wddir = here('repositories/attention-parietal-dopamine/analyses/external/')
source(str_c(wddir,'twolines.R'))

# set files and directories

filedir = here('NCL/gratc_DA/population/')
filename = 'doseResponse_drugMI_';
filename = 'doseResponse_attAUROC_';
drugname = 'SCH23390'
drugname = 'Dopamine'
selectivity = '_none'
# selectivity = '_att&dru'
# selectivity = '_att'
# selectivity = '_dru'

# load data
populationdata = read_csv(str_c(filedir, filename, drugname, selectivity,'.csv'))
populationdata <- as.data.frame(populationdata) 

# run twolines
a = twolines(y~x, data = populationdata[complete.cases(populationdata),])

# make output_mat
a <- as.tibble(sapply(a, "[", 1))
output_mat <- a %>%
  dplyr::select(one_of("b1"), one_of("p1"), one_of("b2"), one_of("p2"), one_of("u.sig")) 

output_mat <- output_mat[1,]

output_mat_colnames <- names(data.frame(output_mat))

output_mat_colnames[which(output_mat_colnames=="u.sig")] <- 'usig'

colnames(output_mat) <- output_mat_colnames

Output_mat <- as_tibble(output_mat)

Output_mat$b1  <- as.character(Output_mat$b1)
Output_mat$b2  <- as.character(Output_mat$b2)
Output_mat$p1  <- as.character(Output_mat$p1)
Output_mat$p2  <- as.character(Output_mat$p2)
Output_mat$usig  <- as.character(Output_mat$usig)

# write
write_csv(Output_mat, str_c(filedir, filename, drugname, selectivity, '_R','.csv'))





