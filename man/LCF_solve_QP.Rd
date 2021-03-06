% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LCF_solve_QP.R
\name{LCF_solve_QP}
\alias{LCF_solve_QP}
\title{Linear combination fitting solve function}
\usage{
LCF_solve_QP(LCF.stds, LCF.samp)
}
\arguments{
\item{LCF.stds}{Standards for LCF}

\item{LCF.samp}{Sample for LCF}
}
\description{
Quadratic programming solution function for linear combination fitting (LCF)
}
\examples{
data(stdmix)
corr.spec.standards  <- initial_load(specdat[1:4], 
  corr.norm = c(-36, -15, 37, 58))
corr.spec.samples    <- initial_load(specdat[5:8], 
  corr.norm = c(-36, -15, 37, 58))
fit.standards <- std_df(sample = corr.spec.samples[[1]], 
  all.standards = corr.spec.standards)
corr.spec <- bkg_corr(raw.spec = corr.spec.samples[[1]], 
  corr.norm = c(-36, -15, 37, 58))
## set fitting range parameters relative to E zero
E.zero <- corr.spec.samples[[1]]$data$E0
LC.pre <- -14
LC.post <- 46
## find ranges that have to be fitted
abs.pre <- abs(corr.spec[["energy"]]-(E.zero+LC.pre))
abs.post <- abs(corr.spec[["energy"]]-(E.zero+LC.post))
range.pre <- which(abs.pre == min(abs.pre))
range.post <- which(abs.post == min(abs.post))
## extract standards and sample in given range
LC.sample <- corr.spec["cor.absorption"][range.pre:range.post,]
LC.standards <- fit.standards[range.pre:range.post,]
## actual fitting
fit.result <- LCF_solve_QP(LCF.stds = LC.standards, LCF.samp = LC.sample)
print(fit.result)
}
\keyword{correction}
\keyword{normalization,}
