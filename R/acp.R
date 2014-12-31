.poissonlik<-function(theta,y,x){ 
n<-nrow(x)
k<-ncol(x)
beta<-theta[1:k]
mu<-exp(x%*%beta)
logl<-(y*(x%*%beta))-mu-log(factorial(y))
sum(-logl)
}

.poissonEst <- function(x, y)
{
k<-ncol(x)
## perform optimization
p<-optim(array(0.1,dim=c(1,k)),.poissonlik,y=y,x=x,method="BFGS",hessian=T)
## compute coefficients
coef<-c(p$par)
## degrees of freedom and standard deciation of residuals 
df <- nrow(x)-ncol(x)
sigma2 <- sum((y - exp(x%*%coef))^2)/df
## compute covariance matrix
vcov<-solve(p$hessian)
colnames(vcov) <- rownames(vcov) <- colnames(x)
list(coefficients = coef,
vcov = vcov,
sigma = sqrt(sigma2),
df = df)
}

.poisson.default <- function(x, y, ...)
{
x <- as.matrix(x)
y <- as.numeric(y)
est <- .poissonEst(x, y)
est$fitted.values <- as.vector(exp(x %*% est$coefficients))
est$residuals <- y - est$fitted.values
est$call <- match.call()
class(est) <- "acp"
est
}


.acpliknc<-function(theta,y){ 
n<-nrow(y)
lambda <- matrix(NA,n,1)
lambda[1] <- theta[1]+theta[2]*mean(y)+theta[3]*mean(y)
lambda[2:n] <- filter(theta[1] + theta[2] * y[1:(n-1)], filter = theta[3], init=lambda[1], method = "recursive")
mu<-lambda
logl<-(y*log(mu))-mu-lfactorial(y)
sum(-logl)
}

.acplikc<-function(theta,y,x){ 
n<-nrow(x)
k<-ncol(x)
beta<-theta[1:k]
lambda <- matrix(NA,n,1)
lambda[1] <- theta[k+1]+theta[k+2]*mean(y)+theta[k+3]*mean(y)
lambda[2:n] <- filter(theta[k + 1] + theta[k + 2] * y[1:(n-1)], filter = theta[k + 3], init=lambda[1], method = "recursive")
mu<-exp(x%*%beta)*lambda
logl<-(y*log(mu))-mu-lfactorial(y)
sum(-logl)
}


.acpEst <- function(x, y, startval,varopt)
{
n<-nrow(x)

if(sum(x)==0){
if(is.null(startval)){
#arma regression to detect initial values for autoregressive parameters
ar<-arma(y, order = c(1, 1))
asval <- matrix(NA,1,3)
asval[3]<-abs(coef(ar)[2])
asval[2]<-abs(abs(coef(ar)[1])-asval[3])
asval[1]<-(mean(y)*abs(1-asval[2]-asval[3]))
## perform optimization
p<-optim(c(asval),.acpliknc,y=y,method="BFGS",hessian=varopt)
}
else{
p<-optim(startval,.acpliknc,y=y,method="BFGS",hessian=varopt)
}

## compute coefficients
coef<-p$par
## compute conditional mean
lambda <- matrix(NA,n,1)
lambda[1] <- coef[1]+coef[2]*mean(y)+coef[3]*mean(y)
for (t in 2:n) {lambda[t] <- coef[1]+coef[2]*y[t-1]+coef[3]*lambda[t-1]}
mu<-lambda
## compute standardized residuals
pres<-(y-mu)/sqrt(mu)
## degrees of freedom and standard deviation of residuals 
df <- nrow(x)-3
sigma2 <- sum(pres^2)/df
if(is.null(p$hessian))
{
vcov<- matrix(NA,1,3)
## create vector of parameter names
namevec <- matrix(NA,1,3)
namevec[1]<-"a"
namevec[2]<-"b"
namevec[3]<-"c"
list(coefficients = coef, vcov = vcov,
df = df)
}
else
{
## compute covariance matrix
vcov<-solve(p$hessian)
## create vector of parameter names
namevec <- matrix(NA,1,3)
namevec[1]<-"a"
namevec[2]<-"b"
namevec[3]<-"c"
colnames(vcov) <- rownames(vcov) <- namevec
list(coefficients = coef,
vcov = vcov,
sigma = sqrt(sigma2),
df = df)
}
}
else{
k<-ncol(x)
if(is.null(startval)){
## poisson regression to detect initial values for covariates
pr<-optim(array(0.1,dim=c(1,k)),.poissonlik,y=y,x=x,method="BFGS",hessian=F)
#arma regression to detect initial values for autoregressive parameters
ar<-arma(y, order = c(1, 1))
asval <- matrix(NA,1,3)
asval[3]<-abs(coef(ar)[2])
asval[2]<-abs(abs(coef(ar)[1])-asval[3])
asval[1]<-(mean(y)*abs(1-asval[2]-asval[3]))
## perform optimization
p<-optim(c(pr$par,asval),.acplikc,y=y,x=x,method="BFGS",hessian=varopt)
}
else{
p<-optim(startval,.acplikc,y=y,x=x,method="BFGS",hessian=varopt)
}
## compute coefficients
coef<-p$par
## compute conditional mean
lambda <- matrix(NA,n,1)
lambda[1] <- coef[k+1]+coef[k+2]*mean(y)+coef[k+3]*mean(y)
for (t in 2:n) {lambda[t] <- coef[k+1]+coef[k+2]*y[t-1]+coef[k+3]*lambda[t-1]}
mu<-exp(x%*%coef[1:k])*lambda
## compute standardized residuals
pres<-(y-mu)/sqrt(mu)
## degrees of freedom and standard deviation of residuals 
df <- nrow(x)-ncol(x)-3
sigma2 <- sum(pres^2)/df
if(is.null(p$hessian))
{
vcov<- matrix(NA,1,k+3)
## create vector of parameter names
namevec <- matrix(NA,1,k+3)
namevec[,1:k]<-colnames(x)
namevec[,k+1]<-"a"
namevec[,k+2]<-"b"
namevec[,k+3]<-"c"
list(coefficients = coef, vcov = vcov,
df = df)
}
else
{
## compute covariance matrix
vcov<-solve(p$hessian)
## create vector of parameter names
namevec <- matrix(NA,1,k+3)
namevec[,1:k]<-colnames(x)
namevec[,k+1]<-"a"
namevec[,k+2]<-"b"
namevec[,k+3]<-"c"
colnames(vcov) <- rownames(vcov) <- namevec
list(coefficients = coef,
vcov = vcov,
sigma = sqrt(sigma2),
df = df)
}
}
}

acp <- function(x, ...) UseMethod("acp")



acp.default <- function(x, y, startval, varopt,...)
{
n<-nrow(x)
x <- as.matrix(x)
y <- as.matrix(y)
if(sum(x)==0){
n<-nrow(x)
est <- .acpEst(x, y ,startval, varopt)
lambda <- matrix(NA,n,1)
lambda[1] <- est$coefficients[1]+est$coefficients[2]*mean(y)+est$coefficients[3]*mean(y)
for (t in 2:n) {lambda[t] <- est$coefficients[1]+est$coefficients[2]*y[t-1]+est$coefficients[3]*lambda[t-1]}
mu<-lambda
est$fitted.values <- as.vector(mu)
est$residuals <- (y - est$fitted.values)/sqrt(est$fitted.values)
}
else{
k<-ncol(x)
est <- .acpEst(x, y ,startval, varopt)
lambda <- matrix(NA,n,1)
lambda[1] <- est$coefficients[k+1]+est$coefficients[k+2]*mean(y)+est$coefficients[k+3]*mean(y)
for (t in 2:n) {lambda[t] <- est$coefficients[k+1]+est$coefficients[k+2]*y[t-1]+est$coefficients[k+3]*lambda[t-1]}
mu<-exp(x%*%est$coefficients[1:k])*lambda
est$fitted.values <- as.vector(mu)
est$residuals <- (y - est$fitted.values)/sqrt(est$fitted.values)
}
est$call <- match.call()
class(est) <- "acp"
est
}

print.acp <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}


summary.acp <- function(object, ...)
{
if(is.na(sum(object$vcov)))
{
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")
print(object$coefficients)
}
else
{
se <- sqrt(diag(object$vcov))
tval <- coef(object) / se
TAB <- cbind(Estimate = coef(object),
StdErr = se,
t.value = tval,
p.value = 2*pt(-abs(tval), df=object$df))
res <- list(call=object$call,
coefficients=TAB)
class(res) <- "summary.acp"
res
}
}

print.summary.acp <- function(x, ...)
{
if(is.na(sum(x$vcov)))
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}
else
{
cat("Call:\n")
print(x$call)
cat("\n")
printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}
}


acp.formula <- function(formula, data=list(), startval=NULL, varopt=T, family="acp",...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)
if (family=="acp"){
est <- acp.default(x, y, startval, varopt,...)
}
else{
est <- .poisson.default(x, y, ...)
}
est$call <- match.call()
est$formula <- formula
est$data <- data
est
}

predict.acp <- function(object, newxdata=NULL, newydata=NULL,...)
{
if(is.null(newxdata))
y <- fitted(object)
else{
if(!is.null(object$formula)){
## model has been fitted using formula interface
x <- model.matrix(object$formula, newxdata)
}
else{
x <- newxdata
}
if(sum(x)==0){
n<-nrow(x)
lambda <- matrix(NA,n,1)
lambda[1] <- 1
for (t in 2:n) {lambda[t] <- coef(object)[1]+coef(object)[2]*newydata[t-1]+coef(object)[3]*lambda[t-1]}
y <- lambda
}
else{
n<-nrow(x)
k<-ncol(x)

if(sum(x[,1])==nrow(x)){
y <-as.vector(exp(x %*% coef(object)))
}
else{
lambda <- matrix(NA,n,1)
lambda[1] <- 1
for (t in 2:n) {lambda[t] <- coef(object)[k+1]+coef(object)[k+2]*newydata[t-1]+coef(object)[k+3]*lambda[t-1]}
y <- as.vector(exp(x%*%coef(object)[1:k])*lambda)
}
}
}
y
}

forecast <- function(object, sample, ydata,...)
{
newxdata<-object$data
x<-model.matrix(object$formula, newxdata)
n<-nrow(x)
if(sum(x)==0){

vpar<-matrix(NA,3,(n-sample+1))
fordata<-newxdata[c(1:(sample)),]
modfor <- acp(object$formula,data=fordata,varopt=F)
vpar[,1]<-coef(modfor)
lambdatempv<-matrix(NA,n-sample,1)
for (t in 1:(n-sample)){
fordata<-newxdata[c(1:(sample+t)),]
modfor <- acp(object$formula,data=fordata,coef(modfor),varopt=F)
vpar[,t+1]<-coef(modfor)
lambdatemp <- matrix(NA,sample+t,1)
lambdatemp[1] <- mean(ydata)
for (tt in 2:(sample+t)) {lambdatemp[tt] <- coef(modfor)[1]+coef(modfor)[2]*ydata[tt-1]+coef(modfor)[3]*lambdatemp[tt-1]}
lambdatempv[t]<-lambdatemp[sample+t]
}
lambda <- matrix(NA,n-sample,1)
lambda[1]<- mean(ydata)
for (tttt in 2:(n-sample)) {
lambda[tttt] <- vpar[1,tttt]+vpar[2,tttt]*ydata[sample+tttt-1]+vpar[3,tttt]*lambdatempv[tttt-1]
}
y <- as.vector(lambda)

}
else{

if(sum(x[,1])==nrow(x)){

k<-ncol(x)
vpar<-matrix(NA,k,(n-sample+1))
fordata<-newxdata[c(1:(sample)),]
modfor <- acp(object$formula,data=fordata,varopt=F,family="poisson")
vpar[,1]<-coef(modfor)
for (t in 1:(n-sample)){
fordata<-newxdata[c(1:(sample+t)),]
modfor <- acp(object$formula,data=fordata,coef(modfor),varopt=F,family="poisson")
vpar[,t+1]<-coef(modfor)
}
vexpbf<-matrix(NA,n-sample,1)
for (ttt in 1:(n-sample)) {
expbf<-exp(vpar[1:k,ttt]%*%x[sample+ttt,])
vexpbf[ttt]<-expbf
}
y <- as.vector(vexpbf)

}
else{

k<-ncol(x)
vpar<-matrix(NA,k+3,(n-sample+1))
fordata<-newxdata[c(1:(sample)),]
modfor <- acp(object$formula,data=fordata,varopt=F)
vpar[,1]<-coef(modfor)
lambdatempv<-matrix(NA,n-sample,1)
for (t in 1:(n-sample)){
fordata<-newxdata[c(1:(sample+t)),]
modfor <- acp(object$formula,data=fordata,coef(modfor),varopt=F)
vpar[,t+1]<-coef(modfor)
lambdatemp <- matrix(NA,sample+t,1)
lambdatemp[1] <- mean(ydata)
for (tt in 2:(sample+t)) {lambdatemp[tt] <- coef(modfor)[k+1]+coef(modfor)[k+2]*ydata[tt-1]+coef(modfor)[k+3]*lambdatemp[tt-1]}
lambdatempv[t]<-lambdatemp[sample+t]
}
vexpbf<-matrix(NA,n-sample,1)
for (ttt in 1:(n-sample)) {
expbf<-exp(vpar[1:k,ttt]%*%x[sample+ttt,])
vexpbf[ttt]<-expbf
}
lambda <- matrix(NA,n-sample,1)
lambda[1]<- 1
for (tttt in 2:(n-sample)) {
lambda[tttt] <- vpar[k+1,tttt]+vpar[k+2,tttt]*ydata[sample+tttt-1]+vpar[k+3,tttt]*lambdatempv[tttt-1]
}
y <- as.vector(lambda*vexpbf)

}

}
y
}


### function for nonrandomized PIT histogram 
###
### input: 
###   x    observed data 
###   Px   CDF at x 
###   Px1  CDF at x-1 

.pit <- function(x, Px, Px1, n.bins=10, y.max=2.75, my.title="PIT Histogram")
  {
  a.mat <- matrix(0,n.bins,length(x))
  k.vec <- pmax(ceiling(n.bins*Px1),1)
  m.vec <- ceiling(n.bins*Px)
  d.vec <- Px-Px1
  for (i in 1:length(x))
      {
      if (k.vec[i]==m.vec[i]) {a.mat[k.vec[i],i]=1}
      else 
        { 
        a.mat[k.vec[i],i]=((k.vec[i]/n.bins)-Px1[i])/d.vec[i]
        if ((k.vec[i]+1)<=(m.vec[i]-1))
           {for (j in ((k.vec[i]+1):(m.vec[i]-1))) {a.mat[j,i]=(1/(n.bins*d.vec[i]))}}
        a.mat[m.vec[i],i]=(Px[i]-((m.vec[i]-1)/n.bins))/d.vec[i]     
        }
      }
  a <- apply(a.mat,1,sum)
  a <- (n.bins*a)/(length(x))
  p <- (0:n.bins)/n.bins
  PIT <- "Probability Integral Transform"
  RF <- "Relative Frequency"
  plot(p, p, ylim=c(0,y.max), type="n", xlab=PIT, ylab=RF, main=my.title) 
  temp1 <- ((1:n.bins)-1)/n.bins
  temp2 <- ((1:n.bins)/n.bins)
  o.vec <- rep(0,n.bins)
  segments(temp1,o.vec,temp1,a)
  segments(temp1,a,temp2,a)
  segments(temp2,o.vec,temp2,a)
  segments(0,0,1,0)
  }
  
evaluation <- function(ydata, yhatdata,...)
{

## compute Pearson standardized residuals
pres<-(ydata-yhatdata)/sqrt(yhatdata)

### parameter settings for computing scores

kk <- 100000                            ### cut-off for summations 
my.k <- (0:kk) - 1                      ### to handle ranked probability score

lambdahat<-yhatdata

pois.Px <- ppois(ydata,lambdahat)                        ### cumulative probabilities
pois.Px1 <- ppois(ydata-1,lambdahat)
pois.px <- dpois(ydata,lambdahat)                        ### probabilities 


pois.logs <- - log(pois.px)
pois.norm <- sum(dpois(my.k,lambdahat)^2) 
pois.qs <- - 2*pois.px + exp(-2*lambdahat)*besselI(2*lambdahat,0)
pois.sphs <- - pois.px / sqrt(exp(-2*lambdahat)*besselI(2*lambdahat,0))
i.cumsum <- cumsum(ppois(my.k,lambdahat)^2)
ii.sum <- sum((ppois(my.k,lambdahat)-1)^2)
ii.cumsum <- cumsum((ppois(my.k,lambdahat)-1)^2)
pois.rps <- (i.cumsum[ydata+1] + ii.sum - ii.cumsum[ydata+1]) 
pois.dss <- (ydata-lambdahat)^2/lambdahat + log(lambdahat)
pois.ses <- (ydata-lambdahat)^2
pois.mae <- abs(ydata-lambdahat)

scores <- matrix(c(round(mean(pois.logs),2),round(mean(pois.qs),3),round(mean(pois.sphs),3),round(mean(pois.rps),2),round(mean(pois.dss),2),round(mean(pois.ses),1),round(mean(pois.mae),1),round(mean(pois.ses)^0.5,1)),ncol=1,byrow=TRUE)
rownames(scores) <- c("logarithmic score","quadratic score","spherical score","ranked probability score","Dawid-Sebastiani score","squared error score","mean absolute error score","root squared error score")
colnames(scores) <- c("Scores")
scores <- as.table(scores)


###create charts

par(mfrow=c(3,3))
plot(lambdahat, type="l",lwd=2, col="red", main="Model fit" )
lines(ydata,lwd=2,col="yellow")
z <-ppois(ydata-1,lambdahat)+(ppois(ydata,lambdahat)-ppois(ydata-1,lambdahat))*runif(length(ydata))
.pit(ydata, ppois(ydata,lambdahat) , ppois(ydata-1,lambdahat))
acf(z,demean=TRUE,main="ACF of PIT")
acf(z^2,demean=TRUE,main="ACF of Squared PIT")
acf(pres,main="ACF of residuals")
acf(pres^2,main="ACF of Squared residuals")
acf(pres^3,main="ACF of Cubic residuals")
scores
}

Berkowitz <- function(ydata, yhatdata, rep, ...)
{
##Perform Berkowitz test
lambdahat<-yhatdata
berkmat <- matrix(NA,rep,1)
for (r in 1:rep) {
z <-ppois(ydata-1,lambdahat)+(ppois(ydata,lambdahat)-ppois(ydata-1,lambdahat))*runif(length(ydata))
u <-qnorm(z)
ar<-arma(u, order = c(1, 0))
arpar <- coef(ar)
res <- residuals(ar)
arssr <- res[2:length(ydata)]%*%res[2:length(ydata)]
s2hat <- (arssr)/ (length(ydata)-2)
s2ut <- s2hat/(1-(arpar[1]^2))
mut <- arpar[2]/(1-arpar[1])
loglR <- -0.5*( u%*%u )
loglUR <- -0.5*log(s2ut)-(((u[1]-mut)^2)/(2*s2ut)) -(length(res)/2)*log(s2hat)-(1/(2*s2hat))*(arssr)
LRberk <- -2*(loglR-loglUR)
berkmat[r] <- LRberk
}
berkowitz <- matrix(c(pchisq(mean(berkmat), df=3, lower.tail=FALSE)),ncol=1,byrow=TRUE)
rownames(berkowitz) <- c("P-Likelihood ratio")
colnames(berkowitz) <- c("Berkowitz test")
berkowitz <- as.table(berkowitz)
berkowitz
}





