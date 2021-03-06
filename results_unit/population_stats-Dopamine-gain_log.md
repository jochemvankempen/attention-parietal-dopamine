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
## [1] "gain_log"
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
    ## gain_log ~ (1 | unit) + attention + drug_on + unit_class + attention:drug_on +  
    ##     attention:unit_class + drug_on:unit_class + attention:drug_on:unit_class
    ##    Data: populationdata
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    240.2    267.9   -110.1    220.2      108 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.81737 -0.56225  0.09498  0.55010  2.22573 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  unit     (Intercept) 0.5882   0.7670  
    ##  Residual             0.1955   0.4422  
    ## Number of obs: 118, groups:  unit, 31
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error       df t value Pr(>|t|)
    ## (Intercept)                     -1.79651    0.14453 31.20695 -12.430 1.27e-13
    ## attention1                      -0.09972    0.04140 87.38001  -2.409   0.0181
    ## drug_on1                         0.20057    0.04140 87.38001   4.845 5.45e-06
    ## unit_class1                      0.11883    0.14453 31.20695   0.822   0.4172
    ## attention1:drug_on1              0.03329    0.04178 87.87696   0.797   0.4277
    ## attention1:unit_class1          -0.06958    0.04140 87.38001  -1.681   0.0964
    ## drug_on1:unit_class1             0.03481    0.04140 87.38001   0.841   0.4027
    ## attention1:drug_on1:unit_class1  0.06441    0.04178 87.87696   1.542   0.1268
    ##                                    
    ## (Intercept)                     ***
    ## attention1                      *  
    ## drug_on1                        ***
    ## unit_class1                        
    ## attention1:drug_on1                
    ## attention1:unit_class1          .  
    ## drug_on1:unit_class1               
    ## attention1:drug_on1:unit_class1    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) attnt1 drg_n1 unt_c1 attntn1:d_1 attntn1:n_1 d_1:_1
    ## attention1   0.006                                                    
    ## drug_on1     0.006 -0.029                                             
    ## unit_class1  0.102  0.001  0.001                                      
    ## attntn1:d_1 -0.013  0.023  0.023 -0.017                               
    ## attntn1:n_1  0.001  0.144 -0.046  0.006  0.006                        
    ## drg_n1:nt_1  0.001 -0.046  0.144  0.006  0.006      -0.029            
    ## attn1:_1:_1 -0.017  0.006  0.006 -0.013  0.160       0.023       0.023

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
    ##   term                npar   AIC   BIC logLik deviance statistic    df   p.value
    ##   <chr>              <dbl> <dbl> <dbl>  <dbl>    <dbl>     <dbl> <dbl>     <dbl>
    ## 1 model.base             3  255.  263.  -124.     249.    NA        NA  NA      
    ## 2 model.att              4  253.  264.  -123.     245.     3.26      1   7.09e-2
    ## 3 model.drug             5  237.  251.  -113.     227.    18.5       1   1.72e-5
    ## 4 model.unitc            6  238.  255.  -113.     226.     0.730     1   3.93e-1
    ## 5 model.att_drug         7  240.  259.  -113.     226.     0.289     1   5.91e-1
    ## 6 model.att_unitc        8  239.  261.  -112.     223.     2.72      1   9.89e-2
    ## 7 model.drug_unitc       9  241.  265.  -111.     223.     0.631     1   4.27e-1
    ## 8 model.att_drug_un~    10  240.  268.  -110.     220.     2.34      1   1.26e-1

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
| .sig01                            |  0.5327805 |  0.9461775 |
| .sigma                            |  0.3607858 |  0.4929637 |
| (Intercept)                       | -2.0784848 | -1.5157195 |
| attention1                        | -0.1791455 | -0.0199693 |
| drug\_on1                         |  0.1204739 |  0.2819162 |
| unit\_class1                      | -0.1566429 |  0.4081842 |
| attention1:drug\_on1              | -0.0486677 |  0.1135839 |
| attention1:unit\_class1           | -0.1525732 |  0.0116249 |
| drug\_on1:unit\_class1            | -0.0473908 |  0.1183388 |
| attention1:drug\_on1:unit\_class1 | -0.0186709 |  0.1432793 |

### parametric bootstrap model fits

We check the model comparisons using parametric bootstrap to control for
Type I errors, as described in the paper [???A Kenward-Roger Approximation
and Parametric Bootstrap Methods for Tests in Linear Mixed Models The R
Package pbkrtest???](https://www.jstatsoft.org/article/view/v059i09).

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
    ##   id                               p    bootp
    ##   <chr>                        <dbl>    <dbl>
    ## 1 bootstrap.att            0.0444    0.0450  
    ## 2 bootstrap.drug           0.0000531 0.000999
    ## 3 bootstrap.unitc          0.407     0.459   
    ## 4 bootstrap.att_drug       0.527     0.554   
    ## 5 bootstrap.att_unitc      0.255     0.283   
    ## 6 bootstrap.drug_unitc     0.156     0.165   
    ## 7 bootstrap.att_drug_unitc 0.123     0.142

### Kenward-roger approximation

We check the model comparisons using the Kenward-roger approximation for
performing F-tests to control for Type I errors, as described in the
paper [???A Kenward-Roger Approximation and Parametric Bootstrap Methods
for Tests in Linear Mixed Models The R Package
pbkrtest???](https://www.jstatsoft.org/article/view/v059i09).

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
    ##   id                     F         p
    ##   <chr>              <dbl>     <dbl>
    ## 1 KR.att             3.29  0.0734   
    ## 2 KR.drug           20.0   0.0000233
    ## 3 KR.unitc           0.690 0.413    
    ## 4 KR.att_drug        0.281 0.597    
    ## 5 KR.att_unitc       2.64  0.108    
    ## 6 KR.drug_unitc      0.596 0.442    
    ## 7 KR.att_drug_unitc  2.21  0.141

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
    ## [1] 0.6867357 0.6982624 0.7536247 0.7533758 0.7545345 0.7620838 0.7635230
    ## [8] 0.7710614

| Model1 | Model2 |           BF |
|-------:|-------:|-------------:|
|      1 |      2 | 6.039230e-01 |
|      1 |      3 | 8.364080e+03 |
|      1 |      4 | 6.708946e+02 |
|      1 |      5 | 8.535756e+01 |
|      1 |      6 | 5.581543e+01 |
|      1 |      7 | 9.047095e+00 |
|      1 |      8 | 7.053931e+00 |
|      2 |      3 | 1.384958e+04 |
|      2 |      4 | 1.110894e+03 |
|      2 |      5 | 1.413385e+02 |
|      2 |      6 | 9.242144e+01 |
|      2 |      7 | 1.498054e+01 |
|      2 |      8 | 1.168018e+01 |
|      3 |      4 | 8.021140e-02 |
|      3 |      5 | 1.020530e-02 |
|      3 |      6 | 6.673200e-03 |
|      3 |      7 | 1.081700e-03 |
|      3 |      8 | 8.434000e-04 |
|      4 |      5 | 1.272295e-01 |
|      4 |      6 | 8.319550e-02 |
|      4 |      7 | 1.348510e-02 |
|      4 |      8 | 1.051420e-02 |
|      5 |      6 | 6.539014e-01 |
|      5 |      7 | 1.059905e-01 |
|      5 |      8 | 8.263980e-02 |
|      6 |      7 | 1.620895e-01 |
|      6 |      8 | 1.263796e-01 |
|      7 |      8 | 7.796902e-01 |
