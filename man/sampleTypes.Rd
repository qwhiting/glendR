\name{sampleTypes}
\alias{sampleTypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sampleTypes: this function gets all the OCIDs for each sample type
%%  ~~function to do ... ~~
}
\description{This labels the "sample.type" with the corresponding 3 digit code. You have to specify which samples are what with vectors
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleTypes(df, MS.Samples, LD1.Samples, LD2.Samples, Travel.Blanks, Field.Blanks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{a dataframe woth all samples}
  \item{MS.Samples}{a vector or column in dataframe with the names of the matrix spike samples (typically end in _MS)}
  \item{LD1.Samples}{a vector or dataframe column that has the names of the samples ran in the first duplicate}
  \item{LD2.Samples}{a vector or dataframe column that has the names of the samples ran in the second duplicate}
  \item{Travel.Blanks}{a vector that has the names of the travel blanks (usually "TravelBlank#")}
  \item{Field.Blanks}{a vector with the field blanks' names}
}
\details{Make sure the names match what was in the SCIEX output, mismatching can occur due to typos!
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
