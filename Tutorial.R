set.seed(12345)
n <- 1000
U <- rnorm(n)
hist(U)

u <- 2
Ua <- c(U, u)
hist(Ua)
abline(v=u)

print(round(mean(Ua > u), 3))

nEval <- 200
u.candidate <- seq(from=min(U), to=max(U), length = nEval)
Cbounds <- c(u.candidate[match(5, apply(as.matrix(u.candidate),
                                        1,
                                        function(x){floor(2.5+100*(mean(c(U, x) < x)))}))],
             rev(u.candidate)[match(95, rev(apply(as.matrix(u.candidate),
                                                  1,
                                                  function(x){ceiling(-2.5+100*(mean(c(U, x) < x)))})))])
print(Cbounds)

hist(U)
abline(v=Cbounds)

X <- rnorm(n)
Y <- X + U
regData <- data.frame(X,Y)
fitlm <- lm(Y~X, data=regData)
eVec <- abs(fitlm$residuals)
hist(eVec)

plot(X,Y)
Xnew <- -1.5
abline(fitlm)
abline(v=Xnew, lty="dashed")
muHat <- predict(fitlm, newdata=data.frame(X=Xnew))
points(Xnew, muHat, pch=19, col="red")
C.X <- c(muHat-quantile(eVec, .975), muHat+quantile(eVec, .975))
points(rep(Xnew, 2), C.X, type="l", col="red", lwd=2)

nEval <- 200
yCand <- seq(from=min(Y), to=max(Y), length=nEval)

confPredict <- function(y, Xin){
  nData <- nrow(regData)  
  regData.a <- rbind(regData,c(Xin, y))
  fitlm.a <- lm(Y~X, data=regData.a)
  resOut <- abs(fitlm.a$residuals)
  resOut_new <- resOut[length(resOut)]
  pi.y <- mean(apply(as.matrix(resOut),
                     1,
                     function(x){x<=resOut_new}))
  testResult <- pi.y*(nData+1) <= ceiling(.975*(nData+1))
  return(testResult)
}

Cxa <- range(yCand[sapply(yCand, confPredict, Xin=Xnew)])

plot(X,Y)
abline(fitlm)
abline(v=Xnew, lty="dashed")
points(rep((Xnew+.05), 2), C.X, type="l", col="red", lwd=2)
points(rep(Xnew, 2), Cxa, type="l", col="blue", lwd=3)

splitConfPredict <- function(Xin){
  nData <- nrow(regData)
  regData$index <- 1:nData
  regData$split <- 1
  regData$split[sample(regData$index, floor(nrow(regData)/2), replace=F)] <- 2
  fitlm.spl <- lm(Y~X, data=subset(regData, split==1))
  resOut <- abs(subset(regData, split==2)$Y - predict(fitlm.spl,
                                                      newdata=subset(regData, split==2)))
  kOut <- ceiling(((nData/2)+1)*(.975))
  resUse <- resOut[order(resOut)][kOut]
  Y.hat <- predict(fitlm.spl, newdata=data.frame(X=Xin))
  C.split <- c(Y.hat-resUse,Y.hat+resUse)
  return(C.split)
}

plot(X,Y)
abline(fitlm)
abline(v=Xnew, lty="dashed")
points(rep((Xnew+.05), 2), C.X, type="l", col="red", lwd=2)
points(rep(Xnew, 2), Cxa, type="l", col="blue", lwd=3)
points(rep(Xnew-.05, 2), splitConfPredict(Xnew), type="l", col="green", lwd=3)

aug.over.X <- function(Xval){range(yCand[sapply(yCand, 
                                                confPredict, 
                                                Xin=Xval)])}
Xvals <- -4:4
augBands <- matrix(unlist(lapply(Xvals, aug.over.X)), ncol=2, byrow=T)
plot(X,Y, type="n")
polygon(c(Xvals, rev(Xvals)),
        c(augBands[,1],rev(augBands[,2])),
        border=F,
        col="gray")
points(X,Y)
abline(fitlm)

Uhet <- rnorm(n, sd=exp(.5*X))
Y <- X + Uhet
regData <- data.frame(X,Y)
plot(X, Y)

fitlm.het <- lm(Y~X, data=regData)
regData$absRes <- abs(fitlm.het$residuals)
locFit <- locpol(absRes~X, 
                 data=regData,
                 bw=2)  
plot(regData$X, regData$absRes)
points(locFit$lpFit[,locFit$X],locFit$lpFit[,locFit$Y],type="l",col="blue")

nEval <- 200
yCand <- seq(from=min(Y), to=max(Y), length=nEval)

confPredict.lw <- function(y, Xin, bwArg){
  nData <- nrow(regData)  
  regData.a <- rbind(regData,c(Xin, y))
  fitlm.a <- lm(Y~X, data=regData.a)
  regData.a$resOut.sc <- abs(fitlm.a$residuals)
  xevalVec <- c(rev(seq(from=Xin, to=min(X), by= -0.25)), 
                seq(from=Xin+.25, to=max(X), by=.25))
  locFit <- locpol(resOut.sc~X, 
                   data=regData.a,
                   bw=bwArg,
                   xeval=xevalVec)  
  resScale <- locFit$lpFit[,2][locFit$lpFit[,1]==Xin]
  resOut <- regData.a$resOut.sc*resScale
  resOut_new <- regData.a$resOut.sc[length(resOut)]
  pi.y <- mean(apply(as.matrix(resOut),
                     1,
                     function(x){x<=resOut_new}))
  testResult <- pi.y*(nData+1) <= ceiling(.975*(nData+1))
  return(list(testResult,resScale))
}

aug.over.X.lw <- function(Xval){
  Cout.lw <- range(yCand[unlist(sapply(yCand,
                                       confPredict.lw,
                                       Xin=Xval,
                                       bwArg=2)[1,])])
  return(Cout.lw)  
}

Xvals <- -3:3

augBands.lw <- matrix(unlist(lapply(Xvals, aug.over.X.lw)), ncol=2, byrow=T)
plot(X,Y, type="n")
polygon(c(Xvals, rev(Xvals)),
        c(augBands.lw[,1],rev(augBands.lw[,2])),
        border=F,
        col="gray")
points(X,Y)
abline(fitlm.het)

