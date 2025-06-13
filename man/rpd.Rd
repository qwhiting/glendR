\name{rpd}
\alias{rpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{rpd: relative percent difference calculation for LD samples
%%  ~~function to do ... ~~
}
\description{a function that will calculate the RPD and and create a new column with the RPD values and comments for the corresponding LD1 and LD2 samples. LD samples must be within the same dataframe.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rpd(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe with LD1 and LD2 samples with their detections and mdls}
}
\details{Column names must be "result", "mdl", "Sample.Type", "analyteType", "Sample.Name", and "Component.Name".
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
