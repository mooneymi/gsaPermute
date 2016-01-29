\name{gsaPermute-package}
\alias{gsaPermute-package}
\alias{gsaPermute}
\docType{package}

\title{
%What the package does (short line)
%~~ package title ~~
Gene permutation procedure for unbiased gene set analyses.
}

\description{
%More about what it does (maybe more than one line)
%~~ A concise (1-5 lines) description of the package ~~
This package can be used to create random gene sets of approximately the same composition (same number of genes with the same number of SNPs) as an a priori defined gene set (e.g. from a publicly available pathway database). As part of a gene set (pathway) analysis, this procedure can be used to calculate a null distribution of association statistics for each pathway being tested, which enables adjustment of p-values for pathway size and linkage disequilibrium (correlation between SNPs). 
}

\details{
\tabular{ll}{
Package: \tab gsaPermute\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-25\cr
License: \tab GPLv3\cr
Copyright: \tab (C) 2014 Michael Mooney
}
%~~ An overview of how to use the package, including the most important functions ~~
}

\author{
%Who wrote it
Michael Mooney

Maintainer: Michael Mooney <mooneymi@ohsu.edu>
%Maintainer: Who to complain to <yourfault@somewhere.net>
%~~ The author and/or maintainer of the package ~~
}

%\references{
%~~ Literature or other references for background information ~~
%}

%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
%\keyword{ package }

\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
\code{\link{makeRandomGS}}, \code{\link{snpsToGenes}}, \code{\link{getGeneSizes}}, \code{\link{getGeneSetComposition}}, \code{\link{makeGeneBins}}
}

\examples{
%~~ simple examples of the most important functions ~~
data(sampleGSA)

# Get gene sizes, and bin genes by size
snpCounts = getGeneSizes(snpMap)
geneBins = makeGeneBins(snpCounts, 25)
  
# Get the composition of the target pathway
gsComp = getGeneSetComposition(snpCounts, pathways[[1]])

# Create random pathways (gene sets) that approximately match the composition of the target pathway
numPerms = 100
x = lapply(c(1:numPerms), function(x){gsComp})
names(x) = c(1:numPerms)
randomPathways = lapply(x, function(y){makeRandomGS(y, geneBins=geneBins)})
}
