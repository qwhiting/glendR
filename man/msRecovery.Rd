\name{msRecovery}
\alias{msRecovery}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{msRecovery: a function that will return the matrix spike recovery in a new column called "MSrec"
%%  ~~function to do ... ~~
}
\description{This function subtracts away the associated sample's concentration of each matching analyte then calculates the recovery based on the spikeLevel.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
msRecovery(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe with all the data present, the LMS samples and the associated sampels. Names must match eachother and the LMS sample must end in "_MS"
%%     ~~Describe \code{x} here~~
}
}
\details{be careful to match sample names because that is how this function gets the associated sample. Columns must also be named as follows: "Sample.Name", "Component.Name", "analyteType", "Sample.Type", "result", "spikeLevel", "weightVolumeAnalyzed"
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
