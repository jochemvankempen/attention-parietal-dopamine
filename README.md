# attention-parietal-dopamine

This repository contains the code to run the analyses and reproduce the results described in the paper ["Dopamine influences attentional rate modulation in Macaque posterior parietal cortex"](https://www.biorxiv.org/content/10.1101/2020.05.15.097675v2).

Preprocessed data can be found [here](https://gin.g-node.org/jochemvankempen/Thiele-attention-gratc-LIP-pipette).

## Folder setup
If you would like to run the analyses on the provided data, I recommend setting up the folder structure by 
1. downloading the data
2. in the folder `./Thiele-attention-gratc-LIP-pipette`, create a folder named `repositories`
3. clone the following repositories inside the folder `repositories`
    - https://github.com/jochemvankempen/attention-parietal-dopamine
    - https://github.com/jochemvankempen/plotj
        - note that this repository contains submodules that need to be cloned using `git submodule update --init --recursive`
    - https://github.com/jochemvankempen/gain-variability

This leads to the following data/code folder organisation:

```
Thiele_attention_gratc_V1_V4_laminar 
│   LICENCE.txt 
│   README.md 
│   datacite.yml 
│
└─── data
     │
     └─── processed
     └─── analysed (produced when running `batch_analyses`)
     └─── population (produced when running `batch_population.mlx`)
└─── repositories
     │
     └─── attention-parietal-dopamine
          └─── results (produced when running `batch_population_stats.R`)
     └─── gain-variability
     └─── plotj
```

## Analyses
Code can be found in the folder `./analyses/`. 
- Analyses conducted on each recorded unit/session are run in MATLAB, controled from `batch_analyses.m`.
- The results from these analyses are compiled into population data tables, further analysed and plotted in MATLAB in `batch_population.mlx`.
- Population statistics are conducted and reported on in R markdown language in `population_stats.Rmd`. This script produces a report (html and md) for one combination of parameters (Drug and type of neural activity). In the script `batch_population_stats.R`, we loop over different parameters. The resulting reports are stored in `./results/`. 

## Licence
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

