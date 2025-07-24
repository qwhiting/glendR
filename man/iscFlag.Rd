\name{iscFlag}
\alias{iscFlag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{iscFlag: function to flag the ISC (FVM) and associated samples
%%  ~~function to do ... ~~
}
\description{
This function flags the LVM sample if it fails the recovery and then flags all associated samples (excludes CAL and CLC sampels).
}
\usage{
iscFlag(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
the df
}
}
\details{
Only flags if above MDL.
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
%later
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
