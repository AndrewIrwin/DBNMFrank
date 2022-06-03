\name{DBrank}
\alias{DBrank}
\title{Rank Selection for Non-Negative Matrix Factorization}
\usage{
DBrank(data,k,alpha,distn,sz,inisz)
}
\arguments{
\item{data}{Matrix. The non-negative data. Its rows are different observations and columns are variables.}
\item{k}{Optional. The value where the hypothesis test start.}
\item{alpha}{Optional. The significance level. Default is 0.1.}
\item{distn}{Character. The distribution of the non-negative data. It should be either \emph{"Normal"} or \emph{"Poisson"}.}
\item{sz}{Optional. The bootstrap size.}
\item{inisz}{Optional. The number of initial values used to obtain the true maximum likelihood for NMF.}
}
\value{
\item{rank}{The NMF rank selected by the function.}
\item{pvalue}{The pvalue for the estimated rank.}
}
\description{
The package estimates the rank parameter for Non-negative Matrix Factorization given the non-negative data and its disitribution. The method is based on hypothesis testing, using a deconvolved bootstrap distribution to assess the significance level accurately despite the large amount of optimization error. The distribution of the non-negative data can be either Normal distributed or Poisson distributed.
}
\details{
Our rank selection for NMF is based on sequentially performing the following hypothesis test:

$H_0$: the rank of the feature matrix is $k$.

$H_a$: the rank of the feature matrix is at least $k+1$.

After applying the goodness-of-fit test, if $H_0$ is rejected by significance level 'alpha', let $k=k+1$ and repeat the test until the pvalue is greater than 'alpha'. For our hypothesis test, the test statistic is the likelihood rato. 'inisz' different initial values are used to get the maximum likelihood for rank 'k' NMF and rank 'k+1' NMF. We use a deconvolved parametric bootstrap to obtain the null distribution of the test statistic. The bootstrap size is 'sz'.
}
\examples{
\donttest{
library(NMF)
set.seed(45217)
########generate a rank 2 Poisson NMF data
x=syntheticNMF(50,2,30)
est.rank=DBrank(t(x),k=2,sz=50,inisz=6)
}
}
