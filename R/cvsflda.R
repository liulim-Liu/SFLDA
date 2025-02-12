#' @title Cross-Validation for Sparse Functional Linear Discriminant Analysis (SFLDA)
#'
#' @description This function performs cross-validation for Sparse Functional Linear Discriminant Analysis (SFLDA).
#' It allows for optional parallel processing and can return various metrics, including accuracy or combined metrics.
#'
#' @param Xdata A data frame containing the training data in long format, with the first three columns for `ID`, `Time`, and `Group`, followed by variable columns.
#' @param Y A vector representing the group labels for the training data.
#' @param plotIt A logical value indicating whether to generate discriminant plots. Default is FALSE.
#' @param metrics.choice A string specifying the metric to optimize. Options are "Accuracy" or "CombinedMetrics". Default is "Accuracy".
#' @param Xtestdata A data frame containing the test data. If NULL, training data is used for testing.
#' @param Ytest A vector representing the group labels for the test data. Required if Xtestdata is provided.
#' @param isParallel A logical value indicating whether to perform parallel processing. Default is TRUE.
#' @param ncores An integer specifying the number of cores to use for parallel processing. Default is NULL, which uses half the available cores.
#' @param nfolds An integer specifying the number of folds for cross-validation. Default is 5.
#' @param ngrid An integer specifying the number of tuning grid values. Default is 8.
#' @param standardize A logical value indicating whether to standardize the data to have mean zero and variance one for each time point and each variable. Default is TRUE.
#' @param maxiteration An integer specifying the maximum number of iterations for optimization. Default is 20.
#' @param thresh A numeric value indicating the convergence threshold. Default is 1e-03.
#'
#' @return A list containing:
#'   \item{CVOut}{A matrix of cross-validation results for different tuning parameter values.}
#'   \item{sfldaerror.test}{The estimated classification error for the test data.}
#'   \item{sidaerror.train}{The estimated classification error for the training data.}
#'   \item{hatalpha}{The estimated alpha coefficients from the SFLDA.}
#'   \item{PredictedClass}{The predicted class labels for the test data.}
#'   \item{var.selected}{The variables selected by the SFLDA.}
#'   \item{varname.selected}{The names of the selected variables.}
#'   \item{optTau}{The optimal tuning parameter values.}
#'   \item{gridValues}{The grid of tuning parameter values used.}
#'   \item{myDiscPlot}{A ggplot object of the discriminant plot if \code{plotIt} is TRUE; otherwise NULL.}
#'   \item{InputData}{The original input data used.}
#'
#' @examples
#' # Example usage of cvSFLDA
#' result <- cvSFLDA(Xdata, Y, plotIt=TRUE, metrics.choice="Accuracy",
#'                    Xtestdata=test_data, Ytest=test_labels,
#'                    nfolds=5, ngrid=8)
#'
#' # Accessing results
#' print(result$sfldaerror.test)
#' print(result$myDiscPlot)
#'
#' @import foreach
#' @import doParallel
#' @import ggplot2
#' @import caret
#' @import Matrix
#' @import RSpectra
#' @import CVXR
#' @import reshape2
#' @import mgcv
#' @import nlme
#' @import lme4
#' @import philentropy
#' @import mltools
#'
#' @export
cvSFLDA=function(Xdata=Xdata,Y=Y,plotIt=FALSE, metrics.choice="Accuracy",
                Xtestdata=NULL,Ytest=NULL,isParallel=TRUE,ncores=NULL,
                nfolds=5,ngrid=8,standardize=TRUE,maxiteration=20, thresh=1e-03){

  ###################### set-up #########################
  starttimeall=Sys.time()

  XdataOrig=Xdata
  XtestdataOrig=Xtestdata
  YOrig=Y
  YtestOrig=Ytest
  D = 1

  #If testing data are not provided, the default is to use training data
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
    Ytest=Y
  }

  #check inputs for testing data
  ntestsizes=lapply(Xtestdata, function(x) dim(x)[1])

  if(is.null(plotIt)){
    plotIt=FALSE
  }

  if(is.null(standardize)){
    standardize=TRUE
  }


  # #standardize if true (for train)
  Xstand=list()
  nTime=length(unique(Xdata$time))
  if(standardize==TRUE){
    for(j in 1:nTime){
      myX=scale(as.matrix(Xdata[Xdata$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xstand[[j]]=cbind.data.frame(as.matrix(Xdata[Xdata$time==j,c(1:3)]), myX)
    }
    Xdata=do.call(rbind.data.frame,Xstand)
  }


  # #standardize if true (for test)
  Xteststand=list()
  ntestTime=length(unique(Xtestdata$time))
  if(standardize==TRUE){
    for(j in 1:ntestTime){
      myX=scale(as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xteststand[[j]]=cbind.data.frame(as.matrix(Xtestdata[Xtestdata$time==j,c(1:3)]), myX)
    }
    Xtestdata=do.call(rbind.data.frame,Xteststand)
  }

  if(is.null(isParallel)){
    isParallel=TRUE
  }

  if(is.null(nfolds)){
    nfolds=5
  }

  if(is.null(ngrid)){
    ngrid=8
  }


  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(thresh)){
    thresh=1e-03
  }

  ###################### split the folds ######################################
  set.seed(1234)
  nK=length(unique(as.vector(Y))) -1

  nc=length(unique(as.vector(Y)))
  Nn=mat.or.vec(nc,1)
  foldid=list()
  for(i in 1:nc)
  {
    Nn[i]=sum(Y==i)
    mod1=Nn[i]%%nfolds
    if(mod1==0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds)),Nn[i])
    }else if(mod1> 0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds), 1:(Nn[i]%%nfolds)),Nn[i])
    }
  }

  foldid=unlist(foldid)

  #print(Y)
  #print(foldid)


  #################### obtain tuning range common to all K #######################
  starttimetune=Sys.time()
  print('Getting tuning grid values')
  myTauvec=sfldatunerange(Xdata,Y,ngrid,standardize)
  endtimetune=Sys.time()
  print('Completed at time')
  print(endtimetune-starttimetune)

  #define the grid
  mygrid=expand.grid(do.call(cbind,myTauvec))
  gridcomb=dim(mygrid)[1]
  gridValues=mygrid

  gc()

  ################### CV ####################################
  starttimeCV=Sys.time()
  CVOut=matrix(0, nfolds, nrow(gridValues))

  #cross validation
  if(isParallel==TRUE){
    cat("Begin", nfolds,"-folds cross-validation", "\n")
    registerDoParallel()
    if(is.null(ncores)){
      ncores=parallel::detectCores()
      ncores=ceiling(ncores/2)}
    cl=makeCluster(ncores)
    registerDoParallel(cl)
    CVOut=matrix(0, nrow(gridValues), nfolds)

    ## .export=c('minv','myfastLDAnonsparse','mysqrtminv','sflda','sfldaclassify','sfldainner', 'sfldatunerange'),
    mycv=foreach(i = 1:nrow(gridValues), .combine='rbind',.export=c('minv','myfastLDAnonsparse','DiscriminantPlots','mysqrtminv','sflda','sfldaclassify','sfldainner', 'sfldatunerange'),
                 .packages=c('readr','caret', 'Matrix', 'RSpectra', 'CVXR','reshape2','mgcv','nlme','lme4','philentropy','ggplot2','mltools')) %dopar% {
      myTau=gridValues[i,]

      #cat("Begin CV-fold", i, "\n")

      CVOut[i,]= sapply(1:nfolds, function(j){
        testInd0=which(foldid==j)
        testInd = subset(Xdata, time == 1)$id[testInd0]

        testX=subset(Xdata, Xdata$id %in% testInd)
        testY=Y[testInd0]
        trainX=subset(Xdata, !(Xdata$id %in% testInd))
        trainY=Y[-testInd0]

        mysflda=sflda(Xtrain=trainX,Y=trainY,Tau=myTau,Xtestdata=testX,Ytest=testY,
                      plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03)

        if (metrics.choice=="Accuracy"){return(mysflda$AverageError)}
        if (metrics.choice=="CombinedMetrics"){return(mysflda$Metrics)}

      } )
    }
    CVOut=t(mycv)
    stopCluster(cl)
  }else if(isParallel==FALSE){
    cat("Begin", nfolds,"-folds cross-validation", "\n")
    CVOut=matrix(0, nfolds, nrow(gridValues))
    for (j in 1:nfolds){
      testInd0=which(foldid==j)
      testInd = subset(Xdata, time == 1)$id[testInd0]

      testX=subset(Xdata, Xdata$id %in% testInd)
      testY=Y[testInd0]
      trainX=subset(Xdata, !(Xdata$id %in% testInd))
      trainY=Y[-testInd0]

      cat("Begin CV-fold", j, "\n")

      CVOut[j,]= sapply(1:nrow(gridValues), function(itau){
        myTau=gridValues[itau,]
        #print(itau)
        mysflda=sflda(Xtrain=trainX,Y=trainY,Tau=myTau,Xtestdata=testX,Ytest=testY,
                      plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03)

        if (metrics.choice=="Accuracy"){return(mysflda$AverageError)}
        if (metrics.choice=="CombinedMetrics"){return(mysflda$Metrics)}
      } )
    }
  }

  View(CVOut)
  gc()

  ############################## results #######################################
  endtimeCV=Sys.time()
  print('Cross-validation completed at time')
  print(endtimeCV-starttimeCV)

  print('Getting Results......')
  #compute average classification error
  if (metrics.choice=="Accuracy"){
    minEorrInd=max(which(colMeans(CVOut, na.rm = TRUE)==min(colMeans(CVOut, na.rm = TRUE))))
    optTau=gridValues[minEorrInd,]
  } else if (metrics.choice=="CombinedMetrics"){
    maxScoreInd=max(which(colMeans(CVOut, na.rm = TRUE)==max(colMeans(CVOut, na.rm = TRUE))))
    optTau=gridValues[maxScoreInd,]
  } else {
    print("please input correct CV target")
  }

  #Apply on testing data
  moptTau= optTau
  print(moptTau)
  mysfldaTest=sflda(Xtrain=Xdata,Y=Y,Tau=moptTau,Xtestdata=Xtestdata,
                    Ytest=Ytest,plotIt=FALSE,standardize,maxiteration,thresh)

  #Apply on training data
  mysfldaTrain=sflda(Xtrain=Xdata,Y=Y,Tau=moptTau,Xtestdata=Xdata,
                    Ytest=Y,plotIt=FALSE,standardize,maxiteration,thresh)

  #print out some results
  cat("Estimated Test Classification Error is", mysfldaTest$AverageError, "\n")
  cat("Estimated Train Classification Error is", mysfldaTrain$AverageError, "\n")

  endtimeall=Sys.time()
  print("Total time used is")
  print(endtimeall-starttimeall)
  gc()

  ############################## plots #########################################
  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
    if (nk == 2){
      myalpha_plot = mysfldaTest$hatalpha
    } else if (nk == 3){
      myalpha_plot = mysfldaTest$hatalpha[[1]]
    }
    myDiscPlot <- DiscriminantPlots(Xtestdata = Xtestdata,Ytest = Ytest, myalpha = myalpha_plot,
                                    predicted.class = mysfldaTest$PredictedClass)
  }else{
    myDiscPlot=NULL
  }

  ############################## outputs #######################################
  result=list(CVOut=CVOut,sfldaerror.test=mysfldaTest$AverageError,sidaerror.train=mysfldaTrain$AverageError,
              hatalpha=mysfldaTest$hatalpha,PredictedClass=mysfldaTest$PredictedClass,
              var.selected = mysfldaTest$var.selected, varname.selected = mysfldaTest$varname.selected,
              optTau=moptTau,gridValues=gridValues, myDiscPlot = myDiscPlot,
              InputData=XdataOrig)
  class(result)="SFLDA"

  return(result)
}
