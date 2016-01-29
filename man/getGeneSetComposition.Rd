\name{getGeneSetComposition}
\alias{getGeneSetComposition}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
getGeneSetComposition
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function retrieves the size of each gene in the gene set and returns a table specifying the number of genes of each size.
}

\usage{
getGeneSetComposition(geneSnpCounts, geneSet)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{geneSnpCounts}{
  %%     ~~Describe \code{geneSnpCounts} here~~
  A list of gene sizes (number of SNPs), as created with the function \code{\link{getGeneSizes}}.
  }
  \item{geneSet}{
  %%     ~~Describe \code{geneSet} here~~
  A vector of gene IDs.
  }
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a table specifying the number of genes of each size.
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
\code{\link{getGeneSizes}}
}

\examples{
data(sampleGSA)

# Get gene sizes, and bin genes by size
snpCounts = getGeneSizes(snpMap)
  
# Get the composition of the target pathway (gene set)
gsComp = getGeneSetComposition(snpCounts, pathways[[1]])
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
