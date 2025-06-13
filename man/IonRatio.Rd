\name{IonRatio}
\alias{IonRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{IonRatio: function that calculates the area reatio between the quantifer and qualifier ions
%%  ~~function to do ... ~~
}
\description{this function calculates the ion ratio for each analyte that has quantifiers and qualifiers and then calculates the 'recovery' based on the peak areas in the cal standards and the peak areas in the sample.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IonRatio(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe contaning all samples and analytes exported from SCIEX OS. This will return only Quantifiers and Internal Standards.
%%     ~~Describe \code{x} here~~
}
}
\details{this works because the concentrations are the same in the calset and the sampels. Only use for flagging when quantifier is > 2x the MDL.}
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
%fill in later
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
