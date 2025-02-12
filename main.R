rm(list = ls())
gc()

home.path <- getwd()
setwd(home.path)

################### Import Librarys ###################
library(readr)
library(Matrix)
library(RSpectra)
library(CVXR)
library(caret)
library(doParallel)
library(ggplot2)
library(reshape2)
library(philentropy)
library(mgcv)
library(nlme)
library(lme4)
library("dplyr")
library(mltools)

################### Import Functions ###################
source("R/helperfunctions.R")
source("R/sflda.R")
source("R/sfldatunerange.R")
source("R/sfldaclassify.R")
source("R/cvsflda.R")
source("R/discriminationPlot.R")
source("R/simulation.R")
source("R/spline.R")
source("R/bs.generator.R")

################### Data Simulation ###################
## Binary Class Simulation
data <- simulation(n.class = 2, case.num = 1) ## choose case.num from 1-4
View(data$train)
View(data$test)

## Multi Class Simulation
data <- simulation(n.class = 3, case.num = 1) ## choose case.num from 1-4
View(data$train)
View(data$test)

## b-spline estimation
train.pred <- spline.prediction(data$train)
View(train.pred$pred.df)
test.pred <- spline.prediction(data$test)
View(test.pred$pred.df)

################### SFLDA-2-class ###################
## Read data
XTrain <- get(load("data/simTrain_spline_c2_case4.rda"))
XTest <- get(load("data/simTest_spline_c2_case4.rda"))

trainX <- XTrain
trainX$group <- as.integer(trainX$group)
trainY <- trainX$group[trainX$time == 1]
trainY <- as.integer(trainY)

testX <- XTest
testY <- testX$group[testX$time == 1]
testX$group <- as.integer(testX$group)
testY <- as.integer(testY)

## Tuning parameter
cvmod <- cvSFLDA(Xdata=trainX,Y=trainY, metrics.choice="CombinedMetrics") ## Accuracy or CombinedMetrics
cvmod$optTau
gc()

## SFLDA with known hyper-parameter
myTau=cvmod$optTau ## can choose your own
sflda.result <- sflda(Xtrain=trainX,Y=trainY,Tau=myTau,Xtestdata=testX,Ytest=testY,plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03)

## find the predicted classes
sflda.result$PredictedClass

## find the discriminant scores
sflda.result$hatalpha

## find the selected variables
sflda.result$varname.selected

## plot the discriminant plots
myplots <- DiscriminantPlots(Xtestdata=testX,Ytest=testY,
myalpha=sflda.result$hatalpha, predicted.class = sflda.result$PredictedClass)
myplot$discriminant.plot
myplot$density.loess
myplot$density.plot

################### SFLDA-3-class ###################
## Read data
XTrain <- get(load("data/simTrain_spline_c3_case4.rda"))
XTest <- get(load("data/simTest_spline_c3_case4.rda"))

## Rest is the same as the 2-class problem
