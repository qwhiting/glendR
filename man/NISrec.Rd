\name{NISrec}
\alias{NISrec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{NISrec: a function to calculate the NIS recovery in Quality Control, Unknown, and Blank samples.
%%  ~~function to do ... ~~
}
\description{Uses the average area in the calibration to determine the NIS recovery in samples. THIS WILL ONLY RETURN NIS ANALYTES.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NISrec(df, NIS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe of the entire dataset from the SCIEX OS export.}
  \item{NIS}{a vector of the NIS names.}
%%     ~~Describe \code{x} here~~
}

\details{Make sure the NIS vector matches the NIS names exactly. This only retuen the NIS analytes so make sure to make a new df with this function.}
%%  ~~ If necessary, more details than the description above ~~

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
\author{Quinn Whiting
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
% add in later
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
