#####################
## Mosaic Displays ##
#####################

#########################
## Hair Eye Color Data ##
#########################

data(HairEyeColor)

## Basic Mosaic Display ##

HairEye <- margin.table(HairEyeColor, c(1,2))

mosaicplot(HairEye, main = "Basic Mosaic Display of Hair Eye Color data")

## Hair Mosaic Display with Pearson residuals ##
Hair <- margin.table(HairEyeColor,1)
Hair
mHair <- as.table(rep(mean(margin.table(HairEyeColor,1)),4))
names(mHair) <- names(Hair)
mHair

## Pearson residuals from Equiprobability model ##

resid <- (Hair-mHair)/sqrt(mHair)
resid

## First Step in a Mosaic Display ##

mosaicplot(Hair, res=resid, shade = TRUE, main = "Hair Color Proportions")

## Hair Eye Mosais Display with Pearson residuals ##

mosaicplot(HairEye, shade = TRUE, main = " Hair Eye Color with Pearson residuals")

## Show Pearson Residuals ##

(HairEye-loglin(HairEye, c(1,2), fit=T)$fit)/sqrt(loglin(HairEye, c(1,2), fit=T)$fit)

###################
## UKSoccer Data ##
###################

data(UKSoccer)

## UKSoccer Mosaic Display ##

mosaicplot(UKSoccer, shade = TRUE, main = "UK Soccer Scores")

###############################
## Repeat Victimization Data ##
###############################

data(RepVict)

mosaicplot(RepVict[-c(4,7),-c(4,7)], shade = TRUE, main = "Repeat Victimization Data")


##################
## 3-Way Tables ##
##################

## Hair Eye Sex Mosais Display with Pearson residuals ##
mosaicplot(HairEyeColor,shade = TRUE, main = "Hair Eye Color Sex" )

mosaicplot(HairEyeColor, margin = ~Hair*Eye + Sex, main = "Model: (Hair Eye) (Sex)" )

mosaicplot(HairEyeColor, margin = ~Hair*Sex + Eye*Sex, main = "Model: (Hair Sex) (Eye Sex)")


####################
## Premarital Sex ##
####################

data(PreSex)

par(mfrow=c(1,2))

## Mosaic display for Gender and Premarital Sexual Expirience ##

## (Gender Pre) ##
mosaicplot(margin.table(PreSex,c(3,4)), shade = TRUE, clegend = FALSE, main = "Gender and Premarital Sex")

## (Gender Pre)(Extra) ##
mosaicplot(margin.table(PreSex,c(2,3,4)), clegend = FALSE,
           margin = ~Gender*PremaritalSex + ExtramaritalSex , main = "(PreMaritalSex Gender) (Sex)")

par(mfrow=c(1,2))

## (Gender Pre Extra)(Marital) ##
mosaicplot(PreSex, margin = ~Gender*PremaritalSex*ExtramaritalSex + MaritalStatus,
           clegend = FALSE, main = "(PreMarital ExtraMarital) (MaritalStatus)")

## (GPE)(PEM) ##
mosaicplot(PreSex, margin = ~Gender*PremaritalSex*ExtramaritalSex + MaritalStatus*PremaritalSex*ExtramaritalSex, clegend = FALSE, main = "(G P E) (P E M)")



############################
## Employment Status Data ##
############################

data(Employment)

## Employment Status ##
mosaicplot(Employment, margin = ~LayoffCause*EmploymentLength + EmploymentStatus,
           main = "(Layoff Employment) + (EmployStatus)")


mosaicplot(Employment, margin = ~LayoffCause*EmploymentLength + LayoffCause*EmploymentStatus, main = "(Layoff EmpL) (Layoff EmplS)")

par(mfrow=c(1,2))

## Closure ##
mosaicplot(Employment[,,1], shade = TRUE, main = "Layoff : Closure")

## Replaced ##
mosaicplot(Employment[,,2], shade = TRUE, main = "Layoff : Replaced")


#####################
## Mosaic Matrices ##
#####################

data(UCBAdmissions)

mosaicpairs(PreSex)

mosaicpairs(UCBAdmissions)

mosaicpairs(UCBAdmissions, type="conditional")



