# attention-parietal-dopamine

This repository contains the code to run the analyses and reproduce the results described in the paper ["Dopamine influences attentional rate modulation in Macaque posterior parietal cortex"](https://www.biorxiv.org/content/10.1101/2020.05.15.097675v2).

Preprocessed data will be made available [here](https://gin.g-node.org/jochemvankempen/Thiele-attention-gratc-LIP-pipette) upon acceptance of the manuscript.

## Analyses
Code can be found in the folder [./analyses/](https://gitlab.com/JvK/attention-parietal-dopamine/-/tree/master/analyses). 
- Analyses conducted on each recorded unit/session are run in MATLAB, controled from [`batch_analyses.m`](https://gitlab.com/JvK/attention-parietal-dopamine/-/blob/master/analyses/batch_analyses.m).
- The results from these analyses are compiled into population data tables, further analysed and plotted in MATLAB in [`batch_population.mlx`](https://gitlab.com/JvK/attention-parietal-dopamine/-/blob/master/analyses/batch_population.mlx)
- Population statistics are conducted and reported on in R markdown language in [`population_stats.Rmd`](https://gitlab.com/JvK/attention-parietal-dopamine/-/blob/master/analyses/population_stats.Rmd). This script produces a report (html and md) for one combination of parameters (Drug and type of neural activity). In the script [`batch_population_stats.R`](https://gitlab.com/JvK/attention-parietal-dopamine/-/blob/master/analyses/batch_population_stats.R), we loop over different parameters. The resulting reports are stored in [./results/](https://gitlab.com/JvK/attention-parietal-dopamine/-/tree/master/results). 



