\name{MMcc_genes}
\alias{mouse_cell_cycle_genes}
\docType{data}
\title{Mus Musculus genes which are annotated to the cell-cycle}
\description{
  This data set contains Mus Musculus genes which are annotated to the cell-cycle.
  These genes were retrieved from biomart and are intended for use with the 
  "gene_indexer" function. The data set contains the gene names in three different 
  formats, Ensemble Gene IDs (1433 values), MGI symbols (1422 values),
  Entrez Gene IDs (1435 values) and Unigene IDs (1102 values).
  
}
\usage{data("MMcc_genes")}
\format{
  A data set that contains with the following attributes
  \describe{
    \item{\code{mouse_cell_cycle_genes}}{A data frame with four columns corresponding 
      to each of the different ID formats.}
  }
}
\value{A dataset containing genes annotated to the cell-cycle in different ID formats}
\keyword{datasets}
