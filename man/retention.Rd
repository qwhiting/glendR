\name{retention}
\alias{retention}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{retention: calculates the devaition in retention time from the average calset for each PFAS
%%  ~~function to do ... ~~
}
\description{This function used the "used" calset standards and calaculated RT for all PFAS (Quantifiers, Qualifiers, and Internal Standards) for Quality Control, Unknown, and Blank samples.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
retention(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe contaning all data and analytes imported from SCIEX OS}
%%     ~~Describe \code{x} here~~
}

\details{apply a flag if deltaRT is >0.25
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
%add this later
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
