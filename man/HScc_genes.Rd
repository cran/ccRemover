\name{HScc_genes}
\alias{human_cell_cycle_genes}
\docType{data}
\title{Homo Sapien genes which are annotated to the cell-cycle}
\description{
  This data set contains Homo Sapien genes which are annotated to the cell-cycle.
  These genes were retrieved from biomart and are intended for use with the 
  "gene_indexer" function. The data set contains the gene names in four different 
  formats, Ensemble Gene IDs (1838 values), HGNC symbols (1740 values), 
  Entrez Gene IDs (1744 values) and Unigene IDs (1339).

}
\usage{data("HScc_genes")}
\format{
  A data set that contains with the following attributes
  \describe{
    \item{\code{human_cell_cycle_genes}}{A data frame with four columns corresponding 
      to each of the different ID formats.}
  }
}
\value{A data set containing genes annotated to the cell-cycle in different ID formats}
\keyword{datasets}
