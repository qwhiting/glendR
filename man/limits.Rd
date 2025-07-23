\name{limits}
\alias{limits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{limits: a function to flag results based on limits
%%  ~~function to do ... ~~
}
\description{
Flags based on calibration limits, mdl limits, loq limits, and non-detects.
}
\usage{
limits(df, cal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
the main dataframe
}
  \item{cal}{
  a dataframe of the calibration (Component.Name, LLOD, and ULOD in ng/mL, this is created with detectionLimits function).
  }
}
\details{
Component.Name is the column name of the analytes
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
