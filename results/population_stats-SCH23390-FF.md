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
## [1] "SCH23390"
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
```

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
    ## FF ~ (1 | unit) + attention + drug_on + unit_class + attention:drug_on +  
    ##     attention:unit_class + drug_on:unit_class + attention:drug_on:unit_class
    ##    Data: populationdata
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    233.2    253.4   -106.6    213.2       46 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6673 -0.3792 -0.0911  0.2128  5.2559 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  unit     (Intercept) 1.223    1.106   
    ##  Residual             1.919    1.385   
    ## Number of obs: 56, groups:  unit, 14
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error       df t value Pr(>|t|)
    ## (Intercept)                      4.29637    0.42499 14.00000  10.109 8.16e-08
    ## attention1                       0.44330    0.22558 42.00000   1.965   0.0560
    ## drug_on1                         0.08945    0.22558 42.00000   0.397   0.6937
    ## unit_class1                      0.84672    0.42499 14.00000   1.992   0.0662
    ## attention1:drug_on1              0.17175    0.22558 42.00000   0.761   0.4507
    ## attention1:unit_class1           0.24338    0.22558 42.00000   1.079   0.2868
    ## drug_on1:unit_class1            -0.33662    0.22558 42.00000  -1.492   0.1431
    ## attention1:drug_on1:unit_class1 -0.08203    0.22558 42.00000  -0.364   0.7179
    ##                                    
    ## (Intercept)                     ***
    ## attention1                      .  
    ## drug_on1                           
    ## unit_class1                     .  
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
    ## unit_class1 0.571  0.000  0.000                                       
    ## attntn1:d_1 0.000  0.000  0.000  0.000                                
    ## attntn1:n_1 0.000  0.571  0.000  0.000  0.000                         
    ## drg_n1:nt_1 0.000  0.000  0.571  0.000  0.000       0.000             
    ## attn1:_1:_1 0.000  0.000  0.000  0.000  0.571       0.000       0.000

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
    ##   term                  npar   AIC   BIC logLik deviance statistic    df p.value
    ##   <chr>                <dbl> <dbl> <dbl>  <dbl>    <dbl>     <dbl> <dbl>   <dbl>
    ## 1 model.base               3  232.  238.  -113.     226.    NA        NA NA     
    ## 2 model.att                4  231.  239.  -112.     223.     2.24      1  0.134 
    ## 3 model.drug               5  231.  241.  -111.     221.     2.02      1  0.155 
    ## 4 model.unitc              6  230.  242.  -109.     218.     3.49      1  0.0616
    ## 5 model.att_drug           7  231.  245.  -108.     217.     1.27      1  0.260 
    ## 6 model.att_unitc          8  231.  248.  -108.     215.     1.09      1  0.297 
    ## 7 model.drug_unitc         9  231.  250.  -107.     213.     2.16      1  0.141 
    ## 8 model.att_drug_unitc    10  233.  253.  -107.     213.     0.132     1  0.716

### coefficient confidence intervals

Use robust regression (bootstrap) to obtain confidence intervals

``` r
# run bootstrap to produce 95% confidence intervals
confint(model.att_drug_unitc, level = 0.95,
        method = "boot",
        nsim = 5000,
        boot.type = c("perc","basic","norm"),
        FUN = NULL, quiet = FALSE,
        oldNames = TRUE,
        cl = clus)
```

    ## Computing bootstrap confidence intervals ...

    ## 
    ## 59 message(s): boundary (singular) fit: see ?isSingular

    ##                                         2.5 %    97.5 %
    ## .sig01                           0.3112412667 1.5718104
    ## .sigma                           0.9926907078 1.5812407
    ## (Intercept)                      3.4403447191 5.1228886
    ## attention1                      -0.0043023472 0.8729257
    ## drug_on1                        -0.3420812502 0.5311390
    ## unit_class1                     -0.0009630776 1.6826505
    ## attention1:drug_on1             -0.2542204713 0.6158301
    ## attention1:unit_class1          -0.1916724806 0.6713054
    ## drug_on1:unit_class1            -0.7845457529 0.1179412
    ## attention1:drug_on1:unit_class1 -0.5312551298 0.3566060

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
    ##   id                           p bootp
    ##   <chr>                    <dbl> <dbl>
    ## 1 bootstrap.att            0.118 0.136
    ## 2 bootstrap.drug           0.103 0.135
    ## 3 bootstrap.unitc          0.154 0.207
    ## 4 bootstrap.att_drug       0.183 0.200
    ## 5 bootstrap.att_unitc      0.308 0.381
    ## 6 bootstrap.drug_unitc     0.155 0.190
    ## 7 bootstrap.att_drug_unitc 0.878 0.887

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
    ##   id                    F      p
    ##   <chr>             <dbl>  <dbl>
    ## 1 KR.att            2.25  0.141 
    ## 2 KR.drug           1.98  0.168 
    ## 3 KR.unitc          3.40  0.0899
    ## 4 KR.att_drug       1.19  0.281 
    ## 5 KR.att_unitc      0.997 0.324 
    ## 6 KR.drug_unitc     1.96  0.170 
    ## 7 KR.att_drug_unitc 0.113 0.738

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
print(BF)
```

    ## [1] "R2 values for each model fit:"
    ## [1] 0.4015809 0.4329399 0.4598252 0.4610174 0.4771532 0.4906155 0.5163542
    ## [8] 0.5178823
    ## # A tibble: 28 x 3
    ##    Model1 Model2      BF
    ##     <dbl>  <dbl>   <dbl>
    ##  1      1      2 0.528  
    ##  2      1      3 0.294  
    ##  3      1      4 0.0573 
    ##  4      1      5 0.0253 
    ##  5      1      6 0.0107 
    ##  6      1      7 0.00899
    ##  7      1      8 0.00245
    ##  8      2      3 0.557  
    ##  9      2      4 0.109  
    ## 10      2      5 0.0478 
    ## # ... with 18 more rows
