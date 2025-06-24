\name{loq}
\alias{loq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{loq: a function to calculate the LOQ based from the MDL for GLENDA PFAS submission
%%  ~~function to do ... ~~
}
\description{This function uses the mdl colum to calculate the loq, then rounds it to the nearest 1,2,5 x10^n. It will return the dataframe with the additional 'qunatificationLimit' column.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loq(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
A dataframe that must contain the mdl in a column called "mdl"
}
}
\details{
if this function does not work, make sure that the mdl column name is nammed correctly, this function relys on the fact the the column name is "mdl" exactly.
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
