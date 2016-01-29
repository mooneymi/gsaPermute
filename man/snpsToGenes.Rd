\name{snpsToGenes}
\alias{snpsToGenes}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
snpsToGenes
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function maps SNPs to genes based on their genomic location. All SNPs within a specified distance from a gene are assigned to that gene. 
}

\usage{
snpsToGenes(snpLocs, geneLocs, dist = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpLocs}{
  %%     ~~Describe \code{snpLocs} here~~
  A dataframe or matrix specifying SNP locations. It should contain 3 columns: snpName, chr, position.
  }
  \item{geneLocs}{
  %%     ~~Describe \code{geneLocs} here~~
  A dataframe or matrix specifying gene locations. It should contain 4 columns: geneName, chr, start, end. 
  }
  \item{dist}{
  %%     ~~Describe \code{dist} here~~
  An integer specifying a distance from the gene boundaries (number of bases).
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
Returns a list of vectors containing SNP IDs, one for each gene. Genes with no assigned SNPs are not included in the list. 
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


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
# Load SNP locations from a Plink .map file
snpLocs = read.delim(plink_map_file, sep="", header=F, stringsAsFactors=F)
snpLocs = snpLocs[,c(2,1,4)]
colnames(snpLocs) = c("snpName", "chr", "pos")

# Map SNPs to genes
snpMap = snpsToGenes(snpLocs, geneLocs, dist=1000)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
