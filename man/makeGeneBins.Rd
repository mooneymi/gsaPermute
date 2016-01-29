\name{makeGeneBins}
\alias{makeGeneBins}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
makeGeneBins
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function bins genes based on their size (number of SNPs assigned to the gene), and then combines bins to create final bins containing (approximately) a minimum number of genes as specified by the user. 
}

\usage{
makeGeneBins(geneSnpCounts, minBinSize = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneSnpCounts}{
  %%     ~~Describe \code{geneSnpCounts} here~~
  A list specifying the number of SNPs assigned to each gene, as created with the function \code{\link{getGeneSizes}}.
  }
  \item{minBinSize}{
  %%     ~~Describe \code{minBinSize} here~~
  The minimum number of genes to be placed in each bin. The default is 25 genes per bin. 
  }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
It is guaranteed that each bin will contain at least \code{0.75*minBinSize} genes. If \code{0.75*minBinSize} is greater than the number of genes in \code{geneSnpCounts} an error will be returned.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with two components: simpleBins and binBounds.
  \item{simpleBins}{
  A list containing vectors of geneIDs for each gene size (number of SNPs).
  }
  \item{binBounds}{
  A matrix with two columns specifying the final bin boundaries. Each row of the matrix specifies the gene size limits for each final bin. For instance if the values of a row are 50 and 55, there should be at least \code{0.75*minBinSize} genes containing between 50 and 55 SNPs. 
  }
}

%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%%  ~~who you are~~
%}

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
snpCounts = getGeneSizes(snpMap)
geneBins = makeGeneBins(snpCounts, 25)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
