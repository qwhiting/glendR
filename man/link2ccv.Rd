\name{link2ccv}
\alias{link2ccv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{link2ccv: function that will associate each CCV from the batch to the correct samples.
%%  ~~function to do ... ~~
}
\description{this function will assocaite the samples ran before the CCV to the correct CCV, useful for associating the flags of one CCV but not another.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
link2ccv(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe contaning all samples and analytes. Recovery must be calculated and the injection time must be a column.
%%     ~~Describe \code{x} here~~
}
}
\details{this function needs to be ran after the recoveries are already calculated and the falgs for the CCVs are assigned to the CCV but not the associated samples.}
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
