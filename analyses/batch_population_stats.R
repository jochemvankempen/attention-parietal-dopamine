# Run the population_stats.Rmd, multiple times with different parameters

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
setwd('C:/Jochem/repositories/attention-parietal-dopamine/analyses/')

# parameter definition
drugnames = c('Dopamine','SCH23390') #
dependent_variables = c('rate','FF','gain_log')
# 

# define function to call Rmd script
render_one <- function(drugname,dependent_variable) {

  # Run markdown in new session with specific input defined in params
  xfun::Rscript_call(
    rmarkdown::render,
    list(input="population_stats.Rmd",
         output_format = "github_document",
         output_file = paste0('../results/','population_stats','-',drugname,'-',dependent_variable,'.md'),    
         params = list("drugname"=drugname, "dependent_variable"=dependent_variable),
         envir = new.env()
         )
    )
}

# loop through drugnames and dependent_variables
for (drugname in drugnames) {
  for (dependent_variable in dependent_variables) {
    render_one(drugname, dependent_variable)
  }
}

