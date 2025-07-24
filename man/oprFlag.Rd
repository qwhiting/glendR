\name{oprFlag}
\alias{oprFlag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{oprFlag: a unction to flag the OPR samples and the associated samples
%%  ~~function to do ... ~~
}
\description{
This function flags the OPR samples (FLS) and the associated samples based on the analytical batch
}
\usage{
oprFlag(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
the data frame
}
}
\details{
the df must have the resultQcIdentifier, analyteType, analysisDateTime, Sample.Name, and recovery columns. OPR samples must be numbered with a 1 or 2 depending on if they were part of the first or second analytical batch.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Quinn Whiting
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%add later
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
