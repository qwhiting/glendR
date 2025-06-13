\name{detectionLimits}
\alias{detectionLimits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{detectionLimmits: calculates the instrumental LLOD and ULOD for each analyte.
%%  ~~function to do ... ~~
}
\description{This function allows for the lower limmit of detetion (LLOD) and upper limmit of detection (ULOD) to be calcualted from the SCIEX processed data. It uses the 'used' column from SCIEX so make sure you uncheck all the unused calibration standards.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
detectionLimits(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a data frame containig at minimum all of the calibration detections. All rows and all columns should be exported and used here. If the data frame contains other samples (i.e., quality conytrols, unknowns, blanks) that is okay.
%%     ~~Describe \code{x} here~~
}
}
\details{Please make sure that you created the calibration curve before exporting and using this function.
%%  ~~ If necessary, more details than the description above ~~
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
