population\_stats: Linear mixed effect model fit to
attention-parietal-dopamine data
================
Jochem van Kempen
25 June, 2021

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
attention, drug application and unit type on firing rates, fano factors
(FF) and gain variability during application of the general dopaminergic
agonist dopamine and the D1-specific antagonist SCH23390.

## Packages

``` r
x<-c("tidyverse","dplyr","lme4","lmerTest","MASS","pbkrtest","parallel","MuMIn","varhandle")
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
## [1] "rate"
```

## set files and directories

## define the order of the category labels

``` r
# set predictor levels (their order)
label_attention = c('Attend RF','Attend away')
label_drug_onoff = c('Drug on','Drug off')
label_unitclass = c('Narrow','Broad')
```

## load and prep data

``` r
# load data
populationdata <- read_csv(str_c(filedir, filename,'.csv'))
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   .default = col_double(),
    ##   drug_on = col_character(),
    ##   attention = col_character(),
    ##   Subject = col_character(),
    ##   Task = col_character(),
    ##   Date = col_date(format = ""),
    ##   Drug = col_character(),
    ##   unit_class = col_character(),
    ##   attention_drug_on = col_character()
    ## )
    ## i Use `spec()` for the full column specifications.

``` r
# select units from the recordings with drug "params$drugname"
populationdata <- dplyr::filter(populationdata, Drug==params$drugname)

# set categorical variables
populationdata$Drug <- factor(populationdata$Drug, ordered=FALSE)
populationdata$attention <- factor(populationdata$attention, ordered=FALSE, levels=label_attention)
populationdata$drug_on <- factor(populationdata$drug_on, ordered=FALSE, levels=label_drug_onoff)
populationdata$unit_class <- factor(populationdata$unit_class, ordered=FALSE, levels=label_unitclass)
populationdata$unit <- factor(populationdata$unit, ordered=FALSE)

# set contrasts to sum (intercept is grand mean, effect-coded contrasts)
populationdata$attention <- C(populationdata$attention, "contr.sum")
populationdata$drug_on <- C(populationdata$drug_on, "contr.sum")
populationdata$unit_class <- C(populationdata$unit_class, "contr.sum")
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
    ## #   unit_class <fct>, EjectCurrent_centered <dbl>, s_stim <dbl>, s_att <dbl>,
    ## #   s_dru <dbl>, s_dir <dbl>, s_att*dru <dbl>, s_att*dir <dbl>,
    ## #   s_dru*dir <dbl>, s_att*dru*dir <dbl>, attention_drug_on <chr>

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
model.unitc <- update(model.drug, .~. + unit_class) 

# interaction 
model.att_drug <- update(model.unitc, .~. + attention*drug_on) 
model.att_unitc <- update(model.att_drug, .~. + attention*unit_class) 
model.drug_unitc <- update(model.att_unitc, .~. + drug_on*unit_class) 
model.att_drug_unitc <- update(model.att_unitc, .~. + attention*drug_on*unit_class) 

summary(model.att_drug_unitc)
```

    ## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
    ##   method [lmerModLmerTest]
    ## Formula: 
    ## rate ~ (1 | unit) + attention + drug_on + unit_class + attention:drug_on +  
    ##     attention:unit_class + drug_on:unit_class + attention:drug_on:unit_class
    ##    Data: populationdata
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    852.3    880.5   -416.1    832.3      114 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5500 -0.6232  0.0452  0.5310  3.3942 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  unit     (Intercept) 223.36   14.945  
    ##  Residual              18.06    4.249  
    ## Number of obs: 124, groups:  unit, 31
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value Pr(>|t|)
    ## (Intercept)                      18.8819     2.7240 31.0000   6.932 8.94e-08
    ## attention1                        2.6719     0.3834 93.0000   6.969 4.53e-10
    ## drug_on1                         -2.3057     0.3834 93.0000  -6.014 3.53e-08
    ## unit_class1                      -2.5828     2.7240 31.0000  -0.948    0.350
    ## attention1:drug_on1              -0.2715     0.3834 93.0000  -0.708    0.481
    ## attention1:unit_class1           -0.2273     0.3834 93.0000  -0.593    0.555
    ## drug_on1:unit_class1              0.6007     0.3834 93.0000   1.567    0.121
    ## attention1:drug_on1:unit_class1   0.1747     0.3834 93.0000   0.456    0.650
    ##                                    
    ## (Intercept)                     ***
    ## attention1                      ***
    ## drug_on1                        ***
    ## unit_class1                        
    ## attention1:drug_on1                
    ## attention1:unit_class1             
    ## drug_on1:unit_class1               
    ## attention1:drug_on1:unit_class1    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) attnt1 drg_n1 unt_c1 attntn1:d_1 attntn1:n_1 d_1:_1
    ## attention1  0.000                                                     
    ## drug_on1    0.000  0.000                                              
    ## unit_class1 0.097  0.000  0.000                                       
    ## attntn1:d_1 0.000  0.000  0.000  0.000                                
    ## attntn1:n_1 0.000  0.097  0.000  0.000  0.000                         
    ## drg_n1:nt_1 0.000  0.000  0.097  0.000  0.000       0.000             
    ## attn1:_1:_1 0.000  0.000  0.000  0.000  0.097       0.000       0.000

<!-- ### stepwise regression -->
<!-- Because some variables could be correlated, we used an algorithm for forward/backward step- wise model selection (Venables WN, Ripley BD. 2002. Modern Applied Statistics with S. Fourth ed. New York: Springer) to test whether each predictor indeed explained independent variability that is not explained by any of the other predictors. This procedure can eliminate correlated predictors from the final model. -->
<!-- ```{r stepwise regression, results='hold'} -->
<!-- # -->
<!-- # model.full <- lmerTest::lmer(formula = formula(paste(params$dependent_variable," ~ 1 + attention + drug_on + unit_class + attention*drug_on + attention*unit_class + drug_on*unit_class + attention*drug_on*unit_class + (1|unit)")), -->
<!-- #                              data=populationdata, na.action=na.omit, REML=FALSE) # intercept only -->
<!-- # model.step <- lmerTest::step(model.full, direction='both') -->
<!-- model.step <- lmerTest::step(model.att_drug_unitc, direction='both') -->
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
                          model.unitc,
                          model.att_drug,
                          model.att_unitc,
                          model.drug_unitc,
                          model.att_drug_unitc)
tidy(model.comparison)
```

    ## Warning in tidy.anova(model.comparison): The following column names in ANOVA
    ## output were not recognized or transformed: npar

    ## # A tibble: 8 x 9
    ##   term              npar   AIC   BIC logLik deviance statistic    df     p.value
    ##   <chr>            <dbl> <dbl> <dbl>  <dbl>    <dbl>     <dbl> <dbl>       <dbl>
    ## 1 model.base           3  903.  911.  -449.     897.    NA        NA    NA      
    ## 2 model.att            4  876.  887.  -434.     868.    29.2       1     6.44e-8
    ## 3 model.drug           5  847.  861.  -418.     837.    31.1       1     2.44e-8
    ## 4 model.unitc          6  848.  865.  -418.     836.     0.886     1     3.46e-1
    ## 5 model.att_drug       7  849.  869.  -418.     835.     0.552     1     4.58e-1
    ## 6 model.att_unitc      8  851.  873.  -417.     835.     0.341     1     5.59e-1
    ## 7 model.drug_unitc     9  850.  876.  -416.     832.     2.42      1     1.20e-1
    ## 8 model.att_drug_~    10  852.  880.  -416.     832.     0.207     1     6.49e-1

### coefficient confidence intervals

Use robust regression (bootstrap) to obtain confidence intervals

``` r
# run bootstrap to produce 95% confidence intervals
knitr::kable(confint(model.att_drug_unitc, level = 0.95,
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

|                                   |      2.5 % |     97.5 % |
|:----------------------------------|-----------:|-----------:|
| .sig01                            | 10.6762483 | 18.2649089 |
| .sigma                            |  3.4889821 |  4.7310024 |
| (Intercept)                       | 13.6038285 | 24.1722926 |
| attention1                        |  1.9090579 |  3.4466679 |
| drug\_on1                         | -3.0936274 | -1.5451707 |
| unit\_class1                      | -7.9970786 |  2.8470707 |
| attention1:drug\_on1              | -1.0181587 |  0.5068754 |
| attention1:unit\_class1           | -0.9728182 |  0.5496546 |
| drug\_on1:unit\_class1            | -0.1037088 |  1.3443639 |
| attention1:drug\_on1:unit\_class1 | -0.5774732 |  0.9520280 |

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
model_na.unitc <- update(model_na.drug, .~. + unit_class)

# interaction
model_na.att_drug <- update(model_na.unitc, .~. + attention*drug_on)
model_na.att_unitc <- update(model_na.att_drug, .~. + attention*unit_class)
model_na.drug_unitc <- update(model_na.att_unitc, .~. + drug_on*unit_class)
model_na.att_drug_unitc <- update(model_na.att_unitc, .~. + attention*drug_on*unit_class)

# run on multiple cores
nc <- parallel::detectCores()
clus <- parallel::makeCluster(rep("localhost", nc))

bootstrap.att <- pbkrtest::PBmodcomp(model_na.att, model_na.base, cl = clus)
bootstrap.drug <- pbkrtest::PBmodcomp(model_na.drug, model_na.att, cl = clus)
bootstrap.unitc <- pbkrtest::PBmodcomp(model_na.unitc, model_na.drug, cl = clus)
bootstrap.att_drug <- pbkrtest::PBmodcomp(model_na.att_drug, model_na.unitc, cl = clus)
bootstrap.att_unitc <- pbkrtest::PBmodcomp(model_na.att_unitc, model_na.att_drug, cl = clus)
bootstrap.drug_unitc <- pbkrtest::PBmodcomp(model_na.drug_unitc, model_na.att_unitc, cl = clus)
bootstrap.att_drug_unitc <- pbkrtest::PBmodcomp(model_na.att_drug_unitc, model_na.drug_unitc, cl = clus)

tibble(
  id = c("bootstrap.att", "bootstrap.drug", "bootstrap.unitc", "bootstrap.att_drug","bootstrap.att_unitc","bootstrap.drug_unitc","bootstrap.att_drug_unitc"),
  p = c(bootstrap.att$test$p.value[1], bootstrap.drug$test$p.value[1], bootstrap.unitc$test$p.value[1], bootstrap.att_drug$test$p.value[1], bootstrap.att_unitc$test$p.value[1], bootstrap.drug_unitc$test$p.value[1], bootstrap.att_drug_unitc$test$p.value[1]),
  bootp = c(bootstrap.att$test$p.value[2], bootstrap.drug$test$p.value[2], bootstrap.unitc$test$p.value[2], bootstrap.att_drug$test$p.value[2], bootstrap.att_unitc$test$p.value[2], bootstrap.drug_unitc$test$p.value[2], bootstrap.att_drug_unitc$test$p.value[2]),
)
```

    ## # A tibble: 7 x 3
    ##   id                                 p    bootp
    ##   <chr>                          <dbl>    <dbl>
    ## 1 bootstrap.att            0.000000389 0.000999
    ## 2 bootstrap.drug           0.000000214 0.000999
    ## 3 bootstrap.unitc          0.256       0.294   
    ## 4 bootstrap.att_drug       0.487       0.490   
    ## 5 bootstrap.att_unitc      0.770       0.787   
    ## 6 bootstrap.drug_unitc     0.238       0.269   
    ## 7 bootstrap.att_drug_unitc 0.764       0.793

### Kenward-roger approximation

We check the model comparisons using the Kenward-roger approximation for
performing F-tests to control for Type I errors, as described in the
paper [“A Kenward-Roger Approximation and Parametric Bootstrap Methods
for Tests in Linear Mixed Models The R Package
pbkrtest”](https://www.jstatsoft.org/article/view/v059i09).

``` r
KR.att <- pbkrtest::KRmodcomp(model.att, model.base)
KR.drug <- pbkrtest::KRmodcomp(model.drug, model.att)
KR.unitc <- pbkrtest::KRmodcomp(model.unitc, model.drug)
KR.att_drug <- pbkrtest::KRmodcomp(model.att_drug, model.unitc)
KR.att_unitc <- pbkrtest::KRmodcomp(model.att_unitc, model.att_drug)
KR.drug_unitc <- pbkrtest::KRmodcomp(model.drug_unitc, model.att_unitc)
KR.att_drug_unitc <- pbkrtest::KRmodcomp(model.att_drug_unitc, model.drug_unitc)

tibble(
  id = c("KR.att", "KR.drug", "KR.unitc", "KR.att_drug","KR.att_unitc","KR.drug_unitc","KR.att_drug_unitc"),
  F = c(KR.att$stats$Fstat, KR.drug$stats$Fstat, KR.unitc$stats$Fstat, KR.att_drug$stats$Fstat, KR.att_unitc$stats$Fstat, KR.drug_unitc$stats$Fstat, KR.att_drug_unitc$stats$Fstat),
  p = c(KR.att$stats$p.value, KR.drug$stats$p.value, KR.unitc$stats$p.value, KR.att_drug$stats$p.value, KR.att_unitc$stats$p.value, KR.drug_unitc$stats$p.value, KR.att_drug_unitc$stats$p.value),
)
```

    ## # A tibble: 7 x 3
    ##   id                     F            p
    ##   <chr>              <dbl>        <dbl>
    ## 1 KR.att            34.0   0.0000000819
    ## 2 KR.drug           36.2   0.0000000374
    ## 3 KR.unitc           0.841 0.367       
    ## 4 KR.att_drug        0.535 0.466       
    ## 5 KR.att_unitc       0.327 0.569       
    ## 6 KR.drug_unitc      2.32  0.132       
    ## 7 KR.att_drug_unitc  0.194 0.661

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
           dim(coef(model.unitc)$unit)[2],
           dim(coef(model.att_drug)$unit)[2],
           dim(coef(model.att_unitc)$unit)[2],
           dim(coef(model.drug_unitc)$unit)[2],
           dim(coef(model.att_drug_unitc)$unit)[2])

# r2_type = 1 # marginal
r2_type = 2 # conditional

# get r2
r2_model = c(r.squaredGLMM(model.base)[r2_type],
             r.squaredGLMM(model.att)[r2_type],
             r.squaredGLMM(model.drug)[r2_type],
             r.squaredGLMM(model.unitc)[r2_type],
             r.squaredGLMM(model.att_drug)[r2_type],
             r.squaredGLMM(model.att_unitc)[r2_type],
             r.squaredGLMM(model.drug_unitc)[r2_type],
             r.squaredGLMM(model.att_drug_unitc)[r2_type])

#compute bayes factor
BF <- as_tibble(R2.to.bf(n=n, nPred=n_pred, r2=r2_model))
colnames(BF) <- c("Model1","Model2","BF")

print('R2 values for each model fit:')
print(r2_model)
# print(BF, nrows = 100)

knitr::kable(BF)
```

    ## [1] "R2 values for each model fit:"
    ## [1] 0.8627315 0.8997724 0.9282816 0.9282963 0.9287205 0.9289815 0.9308047
    ## [8] 0.9309589

| Model1 | Model2 |           BF |
|-------:|-------:|-------------:|
|      1 |      2 | 6.645246e+06 |
|      1 |      3 | 1.372121e+14 |
|      1 |      4 | 6.431679e+12 |
|      1 |      5 | 4.786734e+11 |
|      1 |      6 | 3.442194e+10 |
|      1 |      7 | 9.720799e+09 |
|      1 |      8 | 7.431187e+08 |
|      2 |      3 | 2.064816e+07 |
|      2 |      4 | 9.678617e+05 |
|      2 |      5 | 7.203245e+04 |
|      2 |      6 | 5.179935e+03 |
|      2 |      7 | 1.462820e+03 |
|      2 |      8 | 1.118271e+02 |
|      3 |      4 | 4.687400e-02 |
|      3 |      5 | 3.488600e-03 |
|      3 |      6 | 2.509000e-04 |
|      3 |      7 | 7.080000e-05 |
|      3 |      8 | 5.400000e-06 |
|      4 |      5 | 7.442430e-02 |
|      4 |      6 | 5.351900e-03 |
|      4 |      7 | 1.511400e-03 |
|      4 |      8 | 1.155000e-04 |
|      5 |      6 | 7.191110e-02 |
|      5 |      7 | 2.030780e-02 |
|      5 |      8 | 1.552500e-03 |
|      6 |      7 | 2.824012e-01 |
|      6 |      8 | 2.158850e-02 |
|      7 |      8 | 7.644630e-02 |
