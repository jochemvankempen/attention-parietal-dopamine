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
## [1] "Dopamine"
## 
## $dependent_variable
## [1] "FF"
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

    ## # A tibble: 124 x 27
    ##    unit  drug_on  attention  rate    FF   roc     MI     gain gain_log Subject
    ##    <fct> <fct>    <fct>     <dbl> <dbl> <dbl>  <dbl>    <dbl>    <dbl> <chr>  
    ##  1 2     Drug off Attend RF  8.85  1.43 0.703 0.180  NaN        NaN    W      
    ##  2 5     Drug off Attend RF 29.0   5.24 0.612 0.0881   0.190     -1.66 W      
    ##  3 6     Drug off Attend RF 44.3   7.48 0.731 0.191    0.139     -1.98 W      
    ##  4 7     Drug off Attend RF 24.2   2.60 0.734 0.164    0.0653    -2.73 W      
    ##  5 10    Drug off Attend RF 67     4.37 0.697 0.0858   0.0466    -3.07 W      
    ##  6 11    Drug off Attend RF 11.0   3.43 0.701 0.219    0.247     -1.40 W      
    ##  7 12    Drug off Attend RF  9.4   2.63 0.668 0.167    0.121     -2.11 W      
    ##  8 13    Drug off Attend RF  9.79  1.75 0.722 0.194    0.0717    -2.63 W      
    ##  9 18    Drug off Attend RF 59.7   7.68 0.863 0.313    0.107     -2.23 W      
    ## 10 22    Drug off Attend RF 52.7   4.13 0.784 0.152    0.0600    -2.81 W      
    ## # ... with 114 more rows, and 17 more variables: Task <chr>, Date <date>,
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
    ## Formula: FF ~ (1 | unit) + attention + drug_on + attention:drug_on
    ##    Data: populationdata
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    447.3    464.2   -217.6    435.3      118 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2564 -0.5147 -0.1234  0.2115  3.7006 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  unit     (Intercept) 1.058    1.029   
    ##  Residual             1.379    1.174   
    ## Number of obs: 124, groups:  unit, 31
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)          3.49005    0.21273 31.00000  16.406   <2e-16 ***
    ## attention1           0.10797    0.10547 93.00000   1.024    0.309    
    ## drug_on1            -0.02305    0.10547 93.00000  -0.219    0.827    
    ## attention1:drug_on1  0.09356    0.10547 93.00000   0.887    0.377    
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
    ##   term            npar   AIC   BIC logLik deviance statistic    df p.value
    ##   <chr>          <dbl> <dbl> <dbl>  <dbl>    <dbl>     <dbl> <dbl>   <dbl>
    ## 1 model.base         3  443.  452.  -219.     437.   NA         NA  NA    
    ## 2 model.att          4  444.  455.  -218.     436.    1.03       1   0.309
    ## 3 model.drug         5  446.  460.  -218.     436.    0.0474     1   0.828
    ## 4 model.att_drug     6  447.  464.  -218.     435.    0.784      1   0.376

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

|                      |      2.5 % |    97.5 % |
|:---------------------|-----------:|----------:|
| .sig01               |  0.6419763 | 1.3462841 |
| .sigma               |  0.9840756 | 1.3266569 |
| (Intercept)          |  3.0808821 | 3.9032674 |
| attention1           | -0.0987986 | 0.3168639 |
| drug\_on1            | -0.2230845 | 0.1832643 |
| attention1:drug\_on1 | -0.1061959 | 0.2983938 |

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
    ##   id                     p bootp
    ##   <chr>              <dbl> <dbl>
    ## 1 bootstrap.att      0.118 0.132
    ## 2 bootstrap.drug     0.856 0.863
    ## 3 bootstrap.att_drug 0.586 0.587

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
    ##   id               F     p
    ##   <chr>        <dbl> <dbl>
    ## 1 KR.att      1.03   0.313
    ## 2 KR.drug     0.0463 0.830
    ## 3 KR.att_drug 0.762  0.385

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
    ## [1] 0.4275599 0.4339049 0.4341941 0.4389582

| Model1 | Model2 |        BF |
|-------:|-------:|----------:|
|      1 |      2 | 0.1665668 |
|      1 |      3 | 0.0185137 |
|      1 |      4 | 0.0038675 |
|      2 |      3 | 0.1111489 |
|      2 |      4 | 0.0232191 |
|      3 |      4 | 0.2089011 |
