\name{symmetricize}
\alias{symmetricize}
\title{Make a matrix symmetric}
\usage{
  symmetricize(matrix,
    method = c("max", "min", "avg", "ld", "ud"))
}
\arguments{
  \item{matrix}{The matrix to make symmatric}

  \item{method}{The method to use to make the}
}
\value{
  The symmetric matrix
}
\description{
  Make the matrix symmetric by making all "mirrored"
  positions consistent. A variety of methods are provided
  to make the matrix symmetrical.
}
\author{
  Jeffrey D. Allen \email{Jeffrey.Allen@UTSouthwestern.edu}
}

