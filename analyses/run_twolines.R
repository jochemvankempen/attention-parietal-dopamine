run_twolines <- function(filedir, filename, drugname, selectivity){
  # load csv in filedir specified by filename, drugname and selectivity. Apply "two-lines" regression, write results to file.  
  
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
  
  }



