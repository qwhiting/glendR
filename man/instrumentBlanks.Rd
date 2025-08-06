\name{instrumentBlanks}
\alias{instrumentBlanks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{instrumentBlanks: a function to flag IBs and associated samples
%%  ~~function to do ... ~~
}
\description{
This function first flags the IBs if they are >0.5* MDL, it then associates each IB to the corresponding samples and flags them if they are not 5x greater than the IB. It also adds the IB comment.
}
\usage{
IBflags(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
the dataframe:
}
}
\details{
Some column names have to be in the original SCIEX format while others must be named in the GLENDA format, they include the following: Sample.Name, Component.Name, resultQcIdentifier, analyteType, results, weightVolumeAnalyzed, analysisDateTime, and mdl. Additionally, results, weightVolumeAalyzed, and mdl must be in a numeric format and analysisDateTime must be in a time format (POSIXct).
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
%add in later
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
