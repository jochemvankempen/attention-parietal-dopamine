# load the relevant libraries
x<-c("tidyverse", "dplyr", "lme4")
lapply(x, function(x) {if (!require(x, character.only=T)) {install.packages(x);require(x)}})

# clear all
rm(list=ls(all=TRUE)) 

# set files and directories
filedir = str_c('/Users/jochemvankempen/NCL/gratc_DA/population/')
filename = 'doseResponse_drugMI_';
# filename = 'doseResponse_attAUROC_';
selectivity = 'none'
# selectivity = 'att&dru'
selectivity = 'att'

# load data
populationdata = read_csv(str_c(filedir, filename, selectivity,'.csv'))
populationdata <- as.data.frame(populationdata) 

# set unit index
populationdata$unit <- 1:dim(populationdata)[1]

# define polynomial predictors
populationdata <- within(populationdata, poly_x <- poly(x,2)) # define orthogonal polynomial

# define semi-poly predictor
idx_tmp = populationdata$idx==1
x = populationdata$poly_x[,1]
x[idx_tmp] = populationdata$poly_x[idx_tmp,2]
populationdata$semipoly_x <- x

# model fit
# ---------
# baseline model
model.intercept <- lmer(
  y ~ 1 + (1|idx), # formula
  data=populationdata, # data
  REML=FALSE # 
  ) 

# test linear polynomial
model.poly1 <- update(model.intercept, .~. +  poly_x[,1])

# test quadratic polynomial
model.poly2 <- update(model.poly1, .~. +  poly_x[,2])

# test quadratic polynomial
model.semipoly <- update(model.poly2, .~. +  semipoly_x)

# model comparison
stat <- anova(model.intercept, model.poly1, model.poly2, model.semipoly)





# write
write_csv(Output_mat, str_c(filedir, filename, drugname, selectivity, '_R','.csv'))



model.intercept <- lm(
  y ~ 1, # formula
  data=populationdata, # data
) 

# test linear polynomial
model.poly1 <- update(model.intercept, .~. +  poly_x[,1])

# test quadratic polynomial
model.poly2 <- update(model.intercept, .~. +  poly_x[,2])

# test quadratic polynomial
model.semipoly <- update(model.poly2, .~. +  semipoly_x)

stat <- anova(model.intercept, model.poly2, model.semipoly)







idx_tmp = populationdata$idx==2
fitdata <- populationdata[idx_tmp,]

# drug 1
model.intercept <- lm(
  y ~ 1, # formula
  data=fitdata, # data
) 

# test linear polynomial
model.poly1 <- update(model.intercept, .~. +  poly_x[,1])

# test quadratic polynomial
model.poly2 <- update(model.poly1, .~. +  poly_x[,2])

# test quadratic polynomial
model.semipoly <- update(model.poly2, .~. +  semipoly_x)

stat <- anova(model.intercept, model.poly1, model.poly2, model.semipoly)








idx_tmp = populationdata$idx==1
fitdata <- populationdata[idx_tmp,]

# drug 1
model.intercept <- lm(
  y ~ 1, # formula
  data=fitdata, # data
) 

# test linear polynomial
# model.poly1 <- update(model.intercept, .~. + x)

# test quadratic polynomial
model.poly2 <- update(model.intercept, .~. + I(x^2))

# test quadratic polynomial
model.semipoly <- update(model.poly2, .~. +  semipoly_x)

stat <- anova(model.intercept, model.poly2, model.semipoly)

