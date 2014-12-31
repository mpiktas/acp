\name{evaluation}
\alias{evaluation}
\title{acp Evaluation}
\description{Evaluation of an acp regression model.}
\usage{
evaluation(ydata, yhatdata,...)
}
\arguments{
\item{ydata}{ a data frame containing the real values of the dependent varible. }
\item{yhatdata}{ a data frame containing the fitted values of the dependent varible. }
\item{\dots}{ not used. }
}
\details{
Diebold et al. (1998) proposed a density evaluation method which consists in computing the sequence of cumulative probability
of the observed counts under the assumed forecast distribution (Probability Transform Integral-PIT). If the density fit is adequate this sequence will be uniformly distributed and will have no-autocorrelation 
left neither in level nor when raised to integer powers. For this purpose intuitive graphical methods such as correlograms on the basis of the usual Bartlett confidence intervals, histograms and quantile-quantile (QQ) plots
are used. In the case of discrete data Heinen et al. (2007) propose the use of a uniform zero-one continued extension of the PIT as suggested by Denuit and Lambert (2005). 
}
\value{
 A group of scores for count model evaluation proposed by Czado et al (2009) along with a series of evaluation plots. 
 More precisely the measures calculated are logarithmic score, quadratic score, spherical score, ranked probability score, 
 Dawid-Sebastiani score, squared error score, mean absolute error score and root squared error score. Relatively to the
 graphical evaluation, sub-plot 1 depicts the predicted relatively to the real values, sub-plot 2 the non-randomized PIT histogram
 (Czado et al, 2009), sub-plots 3 and 4 the first two powers of the demeaned randomized PIT and sub-plots 5 to 7 the first
 three powers of the Pearson standardized residuals. 
}
\author{Siakoulis Vasileios}
\references{
\itemize{
\item {Czado,C., Erhardt, V., Min, A., and Wagner, S., 2007. Zero-inflated generalized Poisson models with regression effects on the mean,dispersion and zero-inflation level applied to patent outsourcing rates. Statistical Modelling, 2007, 125.}
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
evaluation(polio_data[[1]],fitted(mod1))
}
\keyword{evaluation}



