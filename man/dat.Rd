\name{dat}
\alias{dat}
\docType{data}
\title{ A simulated scRNA-Seq data. }
\description{
 This data contains expression levels (log-transformed and centered) for 50 cells and 2000 genes.
 The 50 cells are randomly assigned to two cell types and three cell-cycle stages.
 400 genes are assigned as cell-cycle genes, and the other 1600 genes are control genes.
 For descriptions of how we generated this data, please refer to the paper.
}
\usage{data(dat)}
\format{
  A list that contains the following attributes (only \code{x} and \code{if.cc} are used by ccRemover.main.)
  \describe{
		\item{\code{x}}{the data matrix. rows are genes, and columns are cells. These should be treated as log-transformed and
		centered (each row has mean 0) expression levels.}
		\item{\code{if.cc}}{a vector of values FALSE's or TRUE's, denoting whether the genes are cell-cycle related or control.}
		\item{\code{n}}{the number of cells. n=ncol(x).}
		\item{\code{p}}{the number of genes. p=nrow(x).}
		\item{\code{pc}}{the number of cell-cycle genes. pc=sum(if.cc).}
		\item{\code{ct}}{cell types. a vector of values 1 and 2.}
		\item{\code{cc}}{cell-cycle stages. a vector of values 1, 2, or 3.}
  }
}
\value{A simulated dataset used to demonstrate the application of 
ccRemover}
\keyword{datasets}
