\name{EISrec}
\alias{EISrec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EISrec: calculated the recovery for each EIS standard.
%%  ~~function to do ... ~~
}
\description{This function first calculated the response factor between EIS and NIS analytes in the calibration standards. Then it calculated the EIS recovery based off of the average RF.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EISrec(df, EIS, NIS, matchDF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{the dataframe contaning all the samples and analytes}
  \item{EIS}{a vector of the EIS analytes, this is used to subset the dataframe for all EIS analytes, names must match EXACTLY.}
  \item{NIS}{a vector of the NIS analytes, this is used to subset the dataframe for all NIS analytes, names must match EXACTLY.}
  \item{matchDF}{a dataframe with the first column being the names of the EIS analytes and the second column being the names of the NIS analytes that are used to calculate the RF between the EIS and NIS. Names must match EXACTLY.}
%%     ~~Describe \code{x} here~~
}

\details{matchDF must have the columns names be "NIS" and 'EIS" for the corresponding EIS and NIS analytes. If errors occur, chack that the names match the names in the original dataframe. also make sure that the correct EIS an NIS analytes are being used. This only works if the processing method on the SCIEX OS already puts in the actual concentrations of each EIS/NIS for each sample. This only returns EIS analytes so make sure to create a new df with this function.
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
