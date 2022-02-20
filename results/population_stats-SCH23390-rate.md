population\_stats: Linear mixed effect model fit to
attention-parietal-dopamine data
================
Jochem van Kempen
20 February, 2022

-   [Summary](#summary)
-   [Packages](#packages)
-   [Params](#params)
-   [set files and directories](#set-files-and-directories)
-   [define the order of the category
    labels](#define-the-order-of-the-category-labels)
-   [load and prep data](#load-and-prep-data)
-   [mixed-level model](#mixed-level-model)
    -   [fit models](#fit-models)
    -   [compare models](#compare-models)
    -   [coefficient confidence
        intervals](#coefficient-confidence-intervals)
    -   [parametric bootstrap model
        fits](#parametric-bootstrap-model-fits)
    -   [Kenward-roger approximation](#kenward-roger-approximation)
    -   [Bayes factor](#bayes-factor)

## Summary

Statistical analyses for the paper:

Dopamine influences attentional rate modulation in Macaque posterior
parietal cortex

Jochem van Kempen, Christian Brandt, Claudia Distler, Mark A. Bellgrove,
Alexander Thiele

Application of mixed effect models to test the relationship between
attention and drug application on firing rates, fano factors (FF) and
gain variability during application of the general dopaminergic agonist
dopamine and the D1-specific antagonist SCH23390.

## Packages

``` r
x<-c("tidyverse","dplyr","lme4","lmerTest","MASS","pbkrtest","parallel","MuMIn","varhandle","BayesFactor","caTools")
lapply(x, function(x) {if (!require(x, character.only=T)) {install.packages(x);require(x)}})

# # remove unnecessary variables
# rm.all.but(keep=c("params"), envir=.GlobalEnv, keep_functions=TRUE, gc_limit=100,
#            regex="auto")
```

## Params

``` r
selectivity_criterium = 'att&dru' # which units to analyse

sprintf('Running population_stats.Rmd with params:\n')
print(params)
## [1] "Running population_stats.Rmd with params:\n"
## $drugname
## [1] "SCH23390"
## 
## $dependent_variable
## [1] "rate"
```

## set files and directories

## define the order of the category labels

``` r
# set predictor levels (their order)
label_attention = c('Attend RF','Attend away')
label_drug_onoff = c('Drug on','Drug off')
```

## load and prep data

``` r
# load data
populationdata <- read_csv(str_c(filedir, filename,'.csv'))
```

    ## Rows: 180 Columns: 27
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr   (7): drug_on, attention, Subject, Task, Drug, unit_class, attention_dr...
    ## dbl  (19): unit, rate, FF, roc, MI, gain, gain_log, EjectCurrent, Weight, pe...
    ## date  (1): Date
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# select units from the recordings with drug "params$drugname"
populationdata <- dplyr::filter(populationdata, Drug==params$drugname)

# set categorical variables
populationdata$Drug <- factor(populationdata$Drug, ordered=FALSE)
populationdata$attention <- factor(populationdata$attention, ordered=FALSE, levels=label_attention)
populationdata$drug_on <- factor(populationdata$drug_on, ordered=FALSE, levels=label_drug_onoff)
populationdata$unit <- factor(populationdata$unit, ordered=FALSE)

# set contrasts to sum (intercept is grand mean, effect-coded contrasts)
populationdata$attention <- C(populationdata$attention, "contr.sum")
populationdata$drug_on <- C(populationdata$drug_on, "contr.sum")
print(populationdata)
```

    ## # A tibble: 56 x 27
    ##    unit  drug_on  attention  rate    FF   roc      MI   gain gain_log Subject
    ##    <fct> <fct>    <fct>     <dbl> <dbl> <dbl>   <dbl>  <dbl>    <dbl> <chr>  
    ##  1 32    Drug off Attend RF 29.2   2.87 0.797  0.207  0.0650   -2.73  W      
    ##  2 34    Drug off Attend RF 18.2   1.79 0.760  0.197  0.0356   -3.34  W      
    ##  3 40    Drug off Attend RF 27.7   3.10 0.291 -0.112  0.0703   -2.66  W      
    ##  4 46    Drug off Attend RF 32.8   1.61 0.837  0.186  0.0177   -4.04  W      
    ##  5 51    Drug off Attend RF 17.2   4.79 0.737  0.266  0.202    -1.60  W      
    ##  6 52    Drug off Attend RF 13.2   5.61 0.746  0.310  0.314    -1.16  W      
    ##  7 54    Drug off Attend RF 17.7   3.82 0.712  0.200  0.150    -1.90  W      
    ##  8 55    Drug off Attend RF 37.4   3.47 0.562  0.0446 0.0555   -2.89  W      
    ##  9 57    Drug off Attend RF  8.70  2.71 0.605  0.113  0.219    -1.52  W      
    ## 10 58    Drug off Attend RF 12.5   4.74 0.737  0.310  0.407    -0.900 W      
    ## # ... with 46 more rows, and 17 more variables: Task <chr>, Date <date>,
    ## #   Drug <fct>, EjectCurrent <dbl>, Weight <dbl>, peak_to_trough_time <dbl>,
    ## #   unit_class <chr>, EjectCurrent_centered <dbl>, s_stim <dbl>, s_att <dbl>,
    ## #   s_dru <dbl>, s_dir <dbl>, `s_att*dru` <dbl>, `s_att*dir` <dbl>,
    ## #   `s_dru*dir` <dbl>, `s_att*dru*dir` <dbl>, attention_drug_on <chr>

## mixed-level model

### fit models

Sequentially add predictors to predict neural activity (defined in
`params$dependent_variable`). We add random intercepts for `unit` to
account for the repeated measurements in the data.

``` r
#intercept only
model.base <- lmerTest::lmer(formula = formula(paste(params$dependent_variable," ~ 1 + (1|unit)")), 
                             data = populationdata, 
                             na.action = na.omit, 
                             REML = FALSE) 

# main effects
model.att <- update(model.base, .~. + attention) 
model.drug <- update(model.att, .~. + drug_on) 

# interaction 
model.att_drug <- update(model.drug, .~. + attention*drug_on) 

summary(model.att_drug)
```

    ## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
    ##   method [lmerModLmerTest]
    ## Formula: rate ~ (1 | unit) + attention + drug_on + attention:drug_on
    ##    Data: populationdata
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    355.6    367.7   -171.8    343.6       50 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.75772 -0.44309 -0.08858  0.42880  2.50496 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  unit     (Intercept) 95.36    9.765   
    ##  Residual             11.08    3.329   
    ## Number of obs: 56, groups:  unit, 14
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)          15.9223     2.6476 14.0000   6.014 3.18e-05 ***
    ## attention1            2.5624     0.4448 42.0000   5.760 8.78e-07 ***
    ## drug_on1             -1.3651     0.4448 42.0000  -3.069  0.00376 ** 
    ## attention1:drug_on1  -0.1640     0.4448 42.0000  -0.369  0.71427    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) attnt1 drg_n1
    ## attention1  0.000               
    ## drug_on1    0.000  0.000        
    ## attntn1:d_1 0.000  0.000  0.000

<!-- ### stepwise regression -->
<!-- Because some variables could be correlated, we used an algorithm for forward/backward step- wise model selection (Venables WN, Ripley BD. 2002. Modern Applied Statistics with S. Fourth ed. New York: Springer) to test whether each predictor indeed explained independent variability that is not explained by any of the other predictors. This procedure can eliminate correlated predictors from the final model. -->
<!-- ```{r stepwise regression, results='hold'} -->
<!-- # -->
<!-- # model.full <- lmerTest::lmer(formula = formula(paste(params$dependent_variable," ~ 1 + attention + drug_on + attention*drug_on + (1|unit)")), -->
<!-- #                              data=populationdata, na.action=na.omit, REML=FALSE) # intercept only -->
<!-- # model.step <- lmerTest::step(model.full, direction='both') -->
<!-- model.step <- lmerTest::step(model.att_drug, direction='both') -->
<!-- model.step -->
<!-- ``` -->

### compare models

Now test whether the addition of a predictor provides a better model fit
(explains more variance) using likelihood ratio tests at every stage.

``` r
# compare model fits
model.comparison <- anova(model.base,
                          model.att,
                          model.drug,
                          model.att_drug)
tidy(model.comparison)
```

    ## Warning in tidy.anova(model.comparison): The following column names in ANOVA
    ## output were not recognized or transformed: npar

    ## # A tibble: 4 x 9
    ##   term            npar   AIC   BIC logLik deviance statistic    df     p.value
    ##   <chr>          <dbl> <dbl> <dbl>  <dbl>    <dbl>     <dbl> <dbl>       <dbl>
    ## 1 model.base         3  379.  385.  -187.     373.    NA        NA NA         
    ## 2 model.att          4  360.  368.  -176.     352.    20.9       1  0.00000492
    ## 3 model.drug         5  354.  364.  -172.     344.     8.47      1  0.00361   
    ## 4 model.att_drug     6  356.  368.  -172.     344.     0.136     1  0.713

### coefficient confidence intervals

Use robust regression (bootstrap) to obtain confidence intervals

``` r
# run bootstrap to produce 95% confidence intervals
knitr::kable(confint(model.att_drug, level = 0.95,
                     method = "boot",
                     nsim = 5000,
                     boot.type = c("perc","basic","norm"),
                     FUN = NULL, quiet = FALSE,
                     oldNames = TRUE,
                     cl = clus
                     )
             )
```

    ## Computing bootstrap confidence intervals ...

|                      |     2.5 % |     97.5 % |
|:---------------------|----------:|-----------:|
| .sig01               |  5.611271 | 12.9596620 |
| .sigma               |  2.497407 |  3.8906731 |
| (Intercept)          | 10.735987 | 21.0216234 |
| attention1           |  1.690242 |  3.4283903 |
| drug\_on1            | -2.248832 | -0.4912110 |
| attention1:drug\_on1 | -1.041981 |  0.6843884 |

### parametric bootstrap model fits

We check the model comparisons using parametric bootstrap to control for
Type I errors, as described in the paper [“A Kenward-Roger Approximation
and Parametric Bootstrap Methods for Tests in Linear Mixed Models The R
Package pbkrtest”](https://www.jstatsoft.org/article/view/v059i09).

``` r
# check whether params$dependent_variable has NA values. Bootstrap does not work with NA values
populationdata_na <- populationdata[complete.cases(populationdata),]
# na.omit(populationdata)

# run models again without NA values
#intercept only
model_na.base <- lmerTest::lmer(formula = formula(paste(params$dependent_variable," ~ 1 + (1|unit)")),
                             data = populationdata_na,
                             na.action = na.omit,
                             REML = FALSE)

# main effects
model_na.att <- update(model_na.base, .~. + attention)
model_na.drug <- update(model_na.att, .~. + drug_on)

# interaction
model_na.att_drug <- update(model_na.drug, .~. + attention*drug_on)

# run on multiple cores
nc <- parallel::detectCores()
clus <- parallel::makeCluster(rep("localhost", nc))

bootstrap.att <- pbkrtest::PBmodcomp(model_na.att, model_na.base, cl = clus)
bootstrap.drug <- pbkrtest::PBmodcomp(model_na.drug, model_na.att, cl = clus)
bootstrap.att_drug <- pbkrtest::PBmodcomp(model_na.att_drug, model_na.drug, cl = clus)

tibble(
  id = c("bootstrap.att", "bootstrap.drug", "bootstrap.att_drug"),
  p = c(bootstrap.att$test$p.value[1], bootstrap.drug$test$p.value[1], bootstrap.att_drug$test$p.value[1]),
  bootp = c(bootstrap.att$test$p.value[2], bootstrap.drug$test$p.value[2], bootstrap.att_drug$test$p.value[2]),
)
```

    ## # A tibble: 3 x 3
    ##   id                         p    bootp
    ##   <chr>                  <dbl>    <dbl>
    ## 1 bootstrap.att      0.0000232 0.000999
    ## 2 bootstrap.drug     0.00480   0.00799 
    ## 3 bootstrap.att_drug 0.938     0.941

### Kenward-roger approximation

We check the model comparisons using the Kenward-roger approximation for
performing F-tests to control for Type I errors, as described in the
paper [“A Kenward-Roger Approximation and Parametric Bootstrap Methods
for Tests in Linear Mixed Models The R Package
pbkrtest”](https://www.jstatsoft.org/article/view/v059i09).

``` r
KR.att <- pbkrtest::KRmodcomp(model.att, model.base)
KR.drug <- pbkrtest::KRmodcomp(model.drug, model.att)
KR.att_drug <- pbkrtest::KRmodcomp(model.att_drug, model.drug)

tibble(
  id = c("KR.att", "KR.drug", "KR.att_drug"),
  F = c(KR.att$stats$Fstat, KR.drug$stats$Fstat, KR.att_drug$stats$Fstat),
  p = c(KR.att$stats$p.value, KR.drug$stats$p.value, KR.att_drug$stats$p.value),
)
```

    ## # A tibble: 3 x 3
    ##   id               F          p
    ##   <chr>        <dbl>      <dbl>
    ## 1 KR.att      26.4   0.00000721
    ## 2 KR.drug      8.94  0.00476   
    ## 3 KR.att_drug  0.126 0.724

### Bayes factor

Compute bayes factor using BayesFactor package. Detailed description can
be found in the paper [An Introduction to Bayesian Hypothesis Testing
for Management
Research](https://journals.sagepub.com/doi/10.1177/0149206314560412).

Interpretation of Bayes Factors is defined
[here](https://journals.sagepub.com/doi/10.1177/0149206314560412) in
[table
1](https://journals.sagepub.com/na101/home/literatum/publisher/sage/journals/content/joma/2015/joma_41_2/0149206314560412/20161119/images/medium/10.1177_0149206314560412-table1.gif).

A nice illustration of the bayes factor is given
[here](https://alexanderetz.com/2015/08/09/understanding-bayes-visualization-of-bf/).

``` r
# number of observations
n = dim(populationdata)[1]

# number of predictors
n_pred = c(dim(coef(model.base)$unit)[2],
           dim(coef(model.att)$unit)[2],
           dim(coef(model.drug)$unit)[2],
           dim(coef(model.att_drug)$unit)[2]
           )

# r2_type = 1 # marginal
r2_type = 2 # conditional

# get r2
r2_model = c(r.squaredGLMM(model.base)[r2_type],
             r.squaredGLMM(model.att)[r2_type],
             r.squaredGLMM(model.drug)[r2_type],
             r.squaredGLMM(model.att_drug)[r2_type]
             )

#compute bayes factor
BF <- as_tibble(R2.to.bf(n=n, nPred=n_pred, r2=r2_model))
colnames(BF) <- c("Model1","Model2","BF")

print('R2 values for each model fit:')
print(r2_model)
# print(BF, nrows = 100)

knitr::kable(BF)
```

    ## [1] "R2 values for each model fit:"
    ## [1] 0.8054289 0.8817419 0.9033728 0.9036848

| Model1 | Model2 |           BF |
|-------:|-------:|-------------:|
|      1 |      2 | 3.221115e+04 |
|      1 |      3 | 4.272048e+05 |
|      1 |      4 | 3.809585e+04 |
|      2 |      3 | 1.326264e+01 |
|      2 |      4 | 1.182691e+00 |
|      3 |      4 | 8.917470e-02 |
