\name{subsets}
\alias{subsets}
\title{
  Choose Sets Of Size k From The n-Elements Of a Vector
}
\description{
Return a matrix where each row is a set of size k chosen from the 
n elements in vector v.  Optional to allow repeated elements.
}
\usage{
subsets(n, k, v=1:n, allow.repeat=F)
}
\arguments{
<!--move the above two lines to just above the first optional argument-->

\item{n}{
Length of element vector.  
}
\item{k}{
Number of elements chosen in the subsets.  
}
\item{v}{
Vector including elements from which to choose.  If not specified, 
integers 1 to n are used.
}
\item{allow.repeat}{
Logical, if (T)rue, include repeats of the same integer.
}
}
\value{
Matrix with k columns and varying  number of rows.   Each row is the
size=k subset of integers 1:n.   If no repeats allowed, it will have 
nCk = n!/((n-k)!k!) rows.  If repeats allowed, (n+k)!/((n-1)!k!) rows.
}
\section{Side Effects}{
<!--describe any side effects if they exist-->
}
\details{
Optional to allow repeats in the subsets.  Sets all sorted lowest 
element  to  highest.   This function  does  almost  the  same as
combinations(), but this function uses recursion, which allows it 
to work very fast.   In using recursion, it runs into problems in 
the amount of data frames  needed for the recursive calls.   This 
limit is reached with an n near 20, and worse with a larger k.
}
\section{References}{
Venebles, W.N. and Ripley, B.D., "Statistics and Computing", New York: 
Springer-Verlag.  2000.  49+
}
\seealso{
combine
}
\examples{
\dontrun{
sets <- subsets(5, 3, v=6:10)
sets.w.rep_subsets(5,3,v=6:10, allow.repeat=T)
  ## sets has 10 rows and sets w/repeats has 35 rows.
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
