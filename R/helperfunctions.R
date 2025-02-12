myfastLDAnonsparse=function(X, Y){


  #   %--------------------------------------------------------------------------
  #   %myfastLDAnonsparse.R: function to obtain nonsparse solution to sparse functional linear
  # discriminant problem
  # %and to obtain matrix needed in constraints
  # %--------------------------------------------------------------------------
  #
  # X is is in long format with a variable name time
  #first column is ID, second column is time, third column is group, fourth-end is
  #variable names
  # if (is.list(X)) {
  #   Xdata2 = X
  # } else {
  #   Xdata2=list(X)
  # }

  #D = length(Xdata2)

  Y=as.vector(Y)

  nc=max(unique(as.vector(Y)))

  Crxd=list()
  Sbx=list()
  myalphaold1=list()
  myalphaold2=list()
  myalphaoldmat=list()
  rmyalphaoldmat=list()
  sqrtminvmat=list()
  tildealphamat=list()
  tildelambda=list()
  Ux=list()
  Swx=list()
  myeigenT=list()
  sqrtminvStrxSbrx=list()

  Xdata=X
  nTime=length(unique(Xdata$time))
  #nTime=5

  print("start fastLDA iteration")

  for (j in 1:nTime){
    #print(j)
    myX=as.matrix(Xdata[Xdata$time==j,-c(1:3)])
    Y=as.matrix(Xdata[Xdata$time==j,3])
    n=dim(myX)[1]
    p=dim(myX)[2]

    mysvd=svd(t(myX));
    Ux[[j]]=mysvd$u;
    V=mysvd$v;
    W=diag(mysvd$d)
    R=W%*%t(V)


    rdata2=cbind(Y,t(R))
    rdata=t(rdata2)
    mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)
    mr=rowMeans(rdata[-1,])
    nc=max(unique(Y))
    C=list()
    for(i in 1:nc)
    {
      C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
    }
    C=as.matrix(do.call(rbind,C))
    Swrx=t(C)%*%C /(n-1)

    Crx=R-rowMeans(R)
    Srx=Crx%*%t(Crx)/(n-1)

    Sbrx=Srx-Swrx
    Sbx[[j]]=Sbrx

    lambda=sqrt(log(p)/n)
    if(n<p){
      Strx=Swrx + lambda*diag(n)

    } else if(n >= p){
      Strx=Swrx
    }

    Swx[[j]]=Strx;

    sqrtminv= mysqrtminv(Strx)$sqrtminv;
    sqrtminvmat[[j]]=sqrtminv;

    #form matrix
    sqrtminvStrxSbrx[[j]]=Ux[[j]]%*%sqrtminv%*%Sbrx%*%sqrtminv%*%t(Ux[[j]]) #uXSw^(-1/2)SbSw^(-1/2)uX^t

    myeigenT[[j]]=eigs_sym(sqrtminv%*%Sbrx%*%sqrtminv,1,which="LM") #pick the first eigenvalue-vector pair

    myalphaold2[[j]]=Ux[[j]]%*%myeigenT[[j]]$vectors
    tildealphamat[[j]]=myalphaold2[[j]]/norm(myalphaold2[[j]],'2')
    tildelambda[[j]]=myeigenT[[j]]$values
    gc()

  }
  print("end fastLDA iteration")
  result=list(tildealphamat=tildealphamat,tildelambda=tildelambda,sqrtminvmat=sqrtminvmat,
              Sbx=Sbx, Ux=Ux, Swx=Swx,SqrtmSwSbSqrtmSw=sqrtminvStrxSbrx);
  return(result)
}


mysqrtminv=function(W){
  #W is symmetric, positive definite
  mysvd=svd(W);
  d=diag(mysvd$d^(-0.5))
  out=mysvd$u%*%d%*%t(mysvd$u)
  result=list(sqrtminv=out)
  return(result)
}

minv=function(X){
  mysvd=svd(X)
  if(length(mysvd$d)==1){
    d=diag(as.matrix(mysvd$d^(-1)))
  }else{
    d=diag(mysvd$d^(-1))
  }
  out=mysvd$v%*%d%*%t(mysvd$u)
  return(out)
}


sfldainner = function(Xtrain,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold, tildelambda,myTau){

  print("start inner")
  if (is.list(myTau)) {
    Tau = myTau
  } else {
    Tau=list(myTau)
  }

  nK=length(unique(as.vector(Y))) -1

  tildelambdaT=do.call(rbind,tildelambda)
  Ux1=as.matrix(do.call(rbind,Ux))

  myhatalpha=list()
  myalphamat=list()

  nTime=length(unique(Xtrain$time))
  #nTime=5
  D=1

  Xdata1=as.matrix(Xtrain[Xtrain$time==1,-c(1:3)])
  p=dim(Xdata1)[2]

  for(d in 1:D){
    for(ii in 1:nK){
      print(paste("solveing ",ii, "-th LD"))

      if(ii==2){
        alphamat = myalphamat[[1]]
        Xnew = data.frame()

        for(j in 1:nTime){
          Xmeta=as.matrix(Xtrain[Xtrain$time==j,1:3])
          Xdata1=as.matrix(Xtrain[Xtrain$time==j,-c(1:3)])
          p=dim(Xdata1)[2]

          ProjmX=Xdata1%*%(as.matrix(alphamat)%*%minv(t(as.matrix(alphamat))%*%as.matrix(alphamat)+0.001)%*%t(as.matrix(alphamat))) #0.001*diag(ii-1)
          ProjmX[is.nan(ProjmX)]=0
          Xn1=Xdata1-ProjmX
          Xn=cbind.data.frame(Xmeta,Xn1)
          Xnew = rbind(Xnew, Xn)
        }
        myfastlda=myfastLDAnonsparse(Xnew, Y)
        Ux1=as.matrix(do.call(rbind,myfastlda$Ux))
        SqrtmSwSbSqrtmSw=myfastlda$SqrtmSwSbSqrtmSw

        tildelambda = myfastlda$tildelambda
        myalphaold=as.matrix(do.call(cbind,myfastlda$tildealphamat))
        print("succuess solving")
      }

      Alphai=list()
      Separation=SqrtmSwSbSqrtmSw
      tildelambdaT=list()

      for(j in 1:nTime){
        Alphaij=Variable(p,1)
        Alphaij2=vec(Alphaij) #vectorize
        Objx=sum(norm2(Alphaij,axis=1))

        Xdata1=as.matrix(Xtrain[Xtrain$time==j,-c(1:3)])
        p=dim(Xdata1)[2]
        tildelambdaT[[j]]=tildelambda[[j]]*matrix(1,nrow=p,ncol=1)

        tryCatch(
          expr = {
            constraints=list(norm_inf(sum_entries(abs(Separation[[j]] %*% as.matrix(myalphaold[,j]) - tildelambdaT[[j]]*Alphaij2), axis=1 ))<= Tau[[d]])
            prob=Problem(Minimize(Objx),constraints)
            result=solve(prob,solver="ECOS")

            alphaj=result$getValue(Alphaij)

            Alphai[[j]]=alphaj
          },
          error = function(e){
            Alphai[[j]]=0
          },
          finally = {
            gc()
          }
        )
      }

      gc()
      alphai = do.call(cbind,Alphai)

    alphai[abs(alphai) <=10^-10]=0

    if((sum(sum(abs(alphai)))==0)){
      myalpha=alphai
    }else{
      numerator <- matrix(rep(t(sqrt(colSums(alphai*alphai))),times=p),ncol=ncol(alphai),byrow=TRUE)
      numerator[,colSums(numerator) == 0] = 1
      myalpha=alphai/numerator
    }

    myalphamat[[ii]]=myalpha
    }

  myhatalpha[[d]]=myalphamat
  }

  print("end inner")
  gc()
  result=list(hatalpha=myhatalpha)
  return(result)
}


