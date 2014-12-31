\name{Berkowitz}
\alias{Berkowitz}
\title{Berkowitz test}
\description{Implements Berkowitz test (2001) for density evaluation.}
\usage{
Berkowitz(ydata, yhatdata, rep, ...)
}
\arguments{
\item{ydata}{ a data frame containing the real values of the dependent varible. }
\item{yhatdata}{ a data frame containing the fitted values of the dependent varible. }
\item{rep}{ number of uniform distirbution drawings. }
\item{\dots}{ not used. }
}
\details{
Diebold et al. (1998) proposed a density evaluation method which consists in computing the sequence of cumulative probability
of the observed counts under the assumed forecast distribution (Probability Transform Integral-PIT). If the density fit is adequate this sequence will be uniformly distributed and will have no-autocorrelation 
left neither in level nor when raised to integer powers. For this purpose intuitive graphical methods such as correlograms on the basis of the usual Bartlett confidence intervals, histograms and quantile-quantile (QQ) plots
are used. In the case of discrete data Heinen et al. (2007) propose the use of a uniform zero-one continued extension as suggested by Denuit and Lambert (2005). 
Finally instead of using graphical tools for detecting uniformity and independence, Berkowitz (2001) applied a formal test for normality and independence of the inverse standard cumulative normal transform of the PIT sequence 
through the estimation of an AR(1) specification and the use of an LR test to the coefficients. 
}
\value{
P-value of the Likelihood Ratio test statistic based on the chi-square distribution with 3 degress of freedom.
}
\author{Siakoulis Vasileios}
\references{
\itemize{
\item {Berkowitz, J., 2001. Testing density forecasts with applications to risk management. American Statistical Association.Journal of Business and Economics Statistics, 19, 4.}
\item {Denuit , M., and Lambert, P., 2005. Constraints on concordance measures in bivariate discrete data. Journal of Multivariate Analysis, 93, 40-57.}
\item {Diebold, F., Gunther, T., and Tay, A., 1998. Evaluating density forecasts with applications financial to risk management. International Economic Review,39, 4, 863-883.}  
\item {Heinen,A., Rengifo, E., 2007. Multivariate autoregressive modeling of time series count data using copulas. Journal of empirical finance 14 (2007) 564-583.}
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
mod1 <- acp(polio~-1+trend+cos12+sin12+cos6+sin6,data=polio_data)
summary(mod1)
Berkowitz(polio_data[[1]],fitted(mod1),50)
}
\keyword{Berkowitz}



