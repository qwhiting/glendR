\name{calPlot}
\alias{calPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calPlot: a plotting function to view the calibration curves for each PFAS
%%  ~~function to do ... ~~
}
\description{This plotting function wraps around each PFAS and shows the R2 for each fit.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calPlot(df, Title, Subtitle)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe of the imported SCIEX data that includes the calibration standards}
  \item{Title}{a character string of the title of the plot}
  \item{Subtitle}{a characher string of the subtitle of the plot, usually a batch ID or the date it was ran}
%%     ~~Describe \code{x} here~~
}

\details{A 1/x weighting factor is applied to the linear fit
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
