\name{cooclassifier}
\alias{cooclassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
COO Classification of DLBCL Gene Expression Profile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to classify diffuse large B cell lymphoma(DLBCL) samples into activated B cell like(ABC) or germinal center B cell like(GCB) subtypes based on gene expression profiles. 
}
\usage{
cooclassifier(data, n=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
a gene-by-sample matrix of gene expression profile of DLBCL.The matrix should have gene symbol as row names and sample id as column names. 
}
  \item{n}{
%%     ~~Describe \code{n} here~~
the number of top most differentially expressed genes in the training set to be used in the predictor. Default is 50.   
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
returns a vector indicating the classified subtype of each sample("ABC","GCB" or "unclassified"). 
}
\references{
%% ~put references to the literature/web site here ~
Wright, G. et al. A gene expression-based method to diagnose clinically distinct subgroups of diffuse large B cell lymphoma. Proc Natl Acad Sci U S A 100, 9991-9996, doi:10.1073/pnas.1732008100 (2003).
}
\author{
%%  ~~who you are~~
Wei Du <weidu.pku@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(staudt)
result <- cooclassifier(staudt)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classification }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
