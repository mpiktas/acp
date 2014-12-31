\name{forecast}
\alias{forecast}
\title{acp Forecast}
\description{ Recursive forecasting of an acp model.}
\usage{
forecast(object, sample, ydata,...)
}
\arguments{
\item{object}{ an object of class \code{"acp"}, i.e., a fitted model. }
\item{sample}{ the point of sample from which the recursive forecast process will begin. }
\item{ydata}{ a data frame containing the real values of the dependent varible. }
\item{\dots}{ not used. }
}
\value{
a series of forecasted values
}
\author{Siakoulis Vasileios}
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
fitfor<-forecast(mod1,158,polio_data[[1]])

}
\keyword{forecast}



