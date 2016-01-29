\name{makeRandomGS}
\alias{makeRandomGS}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
makeRandomGS
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function creates a random gene set of approximately the same composition (same number of genes with the same number of SNPs) as the specified gene set. As part of a gene set (pathway) analysis, this procedure can be used to calculate a null distribution of association statistics for each pathway being tested. 
}

\usage{
makeRandomGS(geneSetComp, geneBins)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{geneSetComp}{
  %%     ~~Describe \code{geneSetComp} here~~
  A table providing the number of genes of each size in the gene set, as created with the function \code{\link{getGeneSetComposition}}.
  }

  \item{geneBins}{
  %%     ~~Describe \code{geneBins} here~~
  A list containing gene bins and bin boundaries, as created with the function \code{\link{makeGeneBins}}.
  }
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a vector of gene IDs randomly selected to approximately match the composition of the specified gene set.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
%%  ~~who you are~~ 
Michael Mooney <mooneymi@ohsu.edu>
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{makeGeneBins}}, \code{\link{getGeneSetComposition}}
}

\examples{
# This example uses the GRASS algorithm implemented in the SNPath package to calculate association statistics

\dontrun{
library(SNPath)
library(snow)
# Create a SNOW cluster for parallel processing
cl <- makeCluster(c("localhost","localhost"), type = "SOCK")}

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

\dontrun{
# Calculate association p-values for each random gene set (see SNPath documentation)
randPathPvals = grass(cl=cl, snp.info=snp.info, gene.info=gene.info, gene.set=randomPathways, y=y, gene.def="abs", dist=1000, nominal.p=T) 
stopCluster(cl)

# Assuming the p-value of the target pathway has already been calculated, 
# the permutation p-value is calculated as the proportion of random pathways 
# with a p-value less than the p-value of the target pathway

permPval = (sum(randPathPvals < targetPathPval)+1)/(numPerms+1)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
