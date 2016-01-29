\name{getGeneSizes}
\alias{getGeneSizes}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
getGeneSizes
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function counts the number of SNPs assigned to each gene, and returns a named list containing the counts.
}

\usage{
getGeneSizes(snpGeneMap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpGeneMap}{
  %%     ~~Describe \code{snpGeneMap} here~~
  A named list containing vectors of SNP IDs mapped to each gene. 
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
Returns a list of SNP counts for each gene.
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
\code{\link{snpsToGenes}}
}

\examples{
data(sampleGSA)

# Get gene sizes
snpCounts = getGeneSizes(snpMap)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
