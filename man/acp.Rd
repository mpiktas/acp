\name{acp}
\alias{acp}
\alias{acp.default}
\alias{acp.formula}
\alias{print.acp}
\alias{predict.acp}
\alias{summary.acp}
\alias{print.summary.acp}
\title{Autoregressive Conditional Poisson (ACP) Regression}
\description{Fit an ACP(1,1) regression model.}
\usage{
acp(x, ...)
\method{acp}{default}(x, y, startval, varopt,...)
\method{acp}{formula}(formula, data=list(), startval=NULL, varopt=T, family="acp",...)
\method{print}{acp}(x, ...)
\method{summary}{acp}(object, ...)
\method{predict}{acp}(object, newxdata=NULL, newydata=NULL,...)
}
\arguments{
\item{x}{ a numeric design matrix for the model. }
\item{y}{ a numeric vector of responses. }
\item{startval}{ a numeric vector of starting values. If not provided the package 
will obtain starting values for the covariate parameters from a poisson regression 
and for the autoregressive parameters from an arma(1,1) regression.}
\item{family}{ A description of the specification to be used. If family="acp" or not provided an Autoregressive Poisson regression 
will be estimated whereas if family="poisson" a plain Poisson regression is provided.}
\item{formula}{ a symbolic description of the model to be fit. }
\item{data}{ an optional data frame containing the variables in the model. }
\item{varopt}{ an optional logical operator T (TRUE) or F (FALSE) determining whether the covariance 
matrix will be calculated (T) or not (F).}
\item{object}{ an object of class \code{"acp"}, i.e., a fitted model. }
\item{newxdata}{ a data frame containing the covariates data upon which a static forecast will be performed.} 
\item{newydata}{ a data frame containing the dependent variable upon which a static forecast will be performed.} 
\item{\dots}{ not used. }
}
\details{
This model has been proposed by Heinen (2003) for cases of count data exhibiting autoregressive behaviour.
As pointed by Cameron and Trivedi (1998), when a count data set exhibits time dependence the plain Poisson 
regression is not adequate. Heinen (2003) proposed the ACP model in close analogy to the Autoregressive Conditional Duration model (ACD) of 
Engle and Russel (1998) and the GARCH model of Bollerslev (1986). The model can be also found in the international bibliography 
as Integer GARCH (Fokianos and Fried, 2010). The most commonly used variant of the model is the ACP(1,1) which is implemented in the package.
}
\value{
An object of class \code{logreg}, basically a list including elements
\item{coefficients}{ a named vector of coefficients }
\item{vcov}{ covariance matrix of coefficients }
\item{fitted.values}{ fitted values }
\item{residuals}{ residuals }
}
\author{Siakoulis Vasileios}
\references{
\itemize{
\item {Bollerslev, T., 1986. Generalized autoregressive conditional heteroskedasticity. Journal of Econometrics 31, 307-327.}
\item {Cameron, A., and Trivedi, P., 1998. Regression Analysis of Count Data. New York: Cambridge University Press.}
\item {Engle, RF., Russell, JR., 1998. Autoregressive conditional duration: a new model for irregularly spaced transaction data. Econometrica 66, 1127-1162.}
\item {Fokianos, K., Fried, R., 2010. Interventions in INGARCH processes. Journal of Time Series Analysis 31/3, 210-225.}
\item {Heinen, A., 2003. Modelling Time Series Count Data: An Autoregressive Conditional Poisson Model. University Library of Munich.MPRA paper, 8113.}
}
}
\examples{

data(polio)

#Create time trend and seasonality variables
trend=(1:168/168)
cos12=cos((2*pi*(1:168))/12)
sin12=sin((2*pi*(1:168))/12)
cos6=cos((2*pi*(1:168))/6)
sin6=sin((2*pi*(1:168))/6)

polio_data<-data.frame(polio, trend , cos12, sin12, cos6, sin6)

#Autoregressive Conditional Poisson Model with explaning covariates
mod1 <- acp(polio~-1+trend+cos12+sin12+cos6+sin6,data=polio_data)
summary(mod1)

##Static out-of-sample fit example
train<-data.frame(polio_data[c(1: 119),])
test<-data.frame(polio_data[c(120:nrow(polio_data)),])
vpolio<-polio_data[[1]]
yfor<-vpolio[120:nrow(polio_data)]
predict(mod1,test,yfor)

#Poisson Model with explaning covariates
mod2 <- acp(polio~trend+cos12+sin12+cos6+sin6,data=polio_data,family="poisson")
summary(mod2)

#Default method for ACP regression
X<-cbind(trend , cos12, sin12, cos6, sin6)
mod3<-acp(X,polio,startval=NULL,varopt=TRUE)
print(mod3)
summary(mod3)
residuals(mod3)
mod3$vcov
}
\keyword{regression}

