#' Removes the effect of the cell-cycle
#'
#' \code{ccRemover} returns a data matrix with the effects of the cell-cycle
#' removed.
#'
#' Implements the algorithm described in Barron, M. & Li, J.
#' "Identifying and removing the cell-cycle effect from scRNA-Sequencing data"
#' (2016), Scientific Reports. This function takes a normalized,
#' log-transformed and centered matrix of scRNA-seq  expression data
#' and a list of genes which are known to be related to the cell-cycle effect.
#' It then captures the main sources of variation in the data and determines
#' which of these are related to the cell-cycle before removing those that are.
#' Please see the original manuscript for further details.
#'
#' @param dat A list containing a data frame , \code{x}, that contains gene expression
#' measurements with each column representing a sample and each row
#' representing a gene and a logical vector, \code{if_cc}, that indicates
#' which of the genes/rows are related to the cell-cycle or factor of interest.
#'
#' It is recommended that the elements of x are log-transformed and centered
#' for each gene. For example if \code{x} contains TPM measurements then we
#' suggest the following two-steps:
#' Step 1: \code{dat$x <- log(dat$x + 1)}
#' Step 2: \code{dat$x} - rowMeans(dat$x)
#' ccRemover requires that the samples have been properly normalized for
#' sequencing depth and we recommend doing so prior to applying the above steps.
#'
#' The \code{if_cc} vector must be the same length as the number of rows in
#' \code{x} and have elements equal to \code{TRUE} for genes which are related
#' to the cell-cycle and and elements equal to \code{FALSE} for genes which
#' are unrelated to the cell-cycle.
#' @param cutoff The significance cutoff for identifying sources of variation
#' related to the cell-cycle. The default value is 3, which roughly corresponds
#' to a p-value of 0.01.
#' @param max_it The maximum number of iterations for the algorithm. The
#' default value is 4.
#' @param nboot The number of bootstrap repititions to be carried out on each
#' iteration to determine the significance of each component.
#' @param ntop The number of components considered tested at each iteration as
#' cell-cycle effects. The default value if 10
#' @param bar Whether to display a progress bar or not. The progress bar will
#' not work in R-markdown enviornments so this option may be turned off. The
#' default value is \code{TRUE}.
#'
#' @return A data matrix with the effects of the cell-cycle removed.
#'
#' @examples
#' set.seed(10)
#' # Load in example data
#' data(t.cell_data)
#' head(t.cell_data[,1:5])
#' # Center data and select small sample for example
#' t_cell_data_cen <- t(scale(t(t.cell_data[,1:20]), center=TRUE, scale=FALSE))
#' # Extract gene names
#' gene_names <- rownames(t_cell_data_cen)
#' # Determine which genes are annotated to the cell-cycle
#' cell_cycle_gene_indices <- gene_indexer(gene_names,
#' species = "mouse", name_type = "symbol")
#' # Create "if_cc" vector
#' if_cc <- rep(FALSE,nrow(t_cell_data_cen))
#' if_cc[cell_cycle_gene_indices] <- TRUE
#' # Move data into list
#' dat <- list(x=t_cell_data_cen, if_cc=if_cc)
#' # Run ccRemover
#' \dontrun{
#'  xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 10)
#' }
#' # Run ccRemover with reduced bootstrap repetitions for example only
#' xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 20, ntop = 10)
#' head(xhat[,1:5])
#' # Run ccRemover with more compoents considered
#' \dontrun{
#' xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 15)
#'  }
#'
ccRemover <- function(dat, cutoff=3, max_it=4, nboot=200, ntop=10, bar=TRUE)
{
  ## check arguments
  if (!is.list(dat)) stop("dat has to be a list!")
  if (is.null(dat$x)) stop("dat has to have an element x!")
  if (is.null(dat$if_cc)) stop("dat has to have an element if.cc!")

  if (!is.vector(dat$if_cc)) stop("dat$if_cc has to be a vector!")
  if (!is.logical(dat$if_cc)) stop("dat$if_cc has to be a vector of
                                   TRUE or FALSE!")

  if (!is.matrix(dat$x) & !is.data.frame(dat$x)) stop("dat$x has to
                                                      be a matrix or data frame")
  if (nrow(dat$x) != length(dat$if_cc)) stop("The number of rows of dat$x has
                                             to be the same as the length of
                                             dat$if.cc!")

  if (sum(dat$x < 0) == 0) warning("All the values of dat$x are
                                   nonnegative--did you perform log(x) or
                                   log(x+1) transformation? This transformation
                                   is ususally necessary.")
  if (sum(rowSums(dat$x) > 1 | rowSums(dat$x) < -1) > 0) warning("It is highly
                                                                 recommended
                                                                 that the rows
                                                                 of dat$x are
                                                                 centered.")
  if (sum(dat$if_cc) < 1) stop("No genes are cell-cycle genes! Please check
                               your dat$if_cc vector")
  if(sum(!dat$if_cc) < 1) stop("All genes are cell-cycle genes! Please check
                               your dat$if_cc vector")
  percentage_cc <- sum(dat$if_cc)/length(dat$if_cc)
  cat(percentage_cc, " of genes are cell-cycle genes")
  if(percentage_cc <= 0.01) warning("Few genes are cell-cycle genes. This may
                                    lead to unreliable estimation of the
                                    cell-cycle effect. Are you sure your
                                    dat$if_cc vector is correct?")
  if(percentage_cc >= 0.99) warning("Almost all genes are cell-cycle genes.
                                    This may lead to unreliable estimation of
                                    the cell-cycle effect. Are you sure your
                                    dat$if_cc vector is correct?")
  ## begin calculation
  xy <- dat$x[dat$if_cc, ]
  xn <- dat$x[!dat$if_cc, ]

  for (i in 1 : max_it)
  {
    cat("\nIteration ", i, "...", fill=TRUE)

    ## calculate the test statistic and the t statistic
    res_boot <- bootstrap_diff(xy=xy, xn=xn, nboot=nboot, bar)
    cat("The bootstrap results on the top", ntop, "components are:")
    print(res_boot[1 : ntop, ])

    ## decide which components to remove
    which_cc <- which(res_boot$t_load_boot[1 : ntop] >= cutoff)

    ## when there is nothing to remove, end the iteration
    if (length(which_cc) == 0)
    {
      cat("No more cell-cycle effect is detected.", fill=TRUE)
      break
    }

    ## remove the cell-cycle effect
    cat("The follow components are removed:", which_cc, fill=TRUE)
    xn_pca <- stats::prcomp(xn, scale.=FALSE)
    xn_hat <- xn
    xy_hat <- xy
    for (i in 1 : length(which_cc))
    {
      xn_hat <- xn_hat - (xn %*% xn_pca$rotation[, which_cc[i]]) %*%
        t(xn_pca$rotation[, which_cc[i]])
      xy_hat <- xy_hat - (xy %*% xn_pca$rotation[, which_cc[i]]) %*%
        t(xn_pca$rotation[, which_cc[i]])
    }
    xn <- xn_hat
    xy <- xy_hat
  }

  xhat <- dat$x
  xhat[dat$if_cc, ] <- xy
  xhat[!dat$if_cc, ] <- xn

  return(xhat)
}


#' Calculates the difference in the loading score for cell-cycle and control
#' genes
#'
#' This function is only used internally inside ccRemover. The function
#' calcualtes the average load difference on the cell-cycle and control genes.
#' Bootstrap resampling is then used to provide a score for each component.
#' Please see the original manuscript for the mathematical details.
#'
#' @param xy The data for the genes which are annotated to the cell-cycle, i.e.
#' those genes for which "if_cc" is \code{TRUE}.
#' @param xn The data for the genes which are not annotated to the cell-cycle,
#' control genes, genes for which "if_cc" is \code{FALSE}
#' @param nboot The number of bootstrap repititions to be carried out on each
#' iteration to determine the significance of each component.
#' @param bar Whether to display a progress bar or not. The progress bar will
#' not work in R-markdown enviornments so this option may be turned off. The
#' default value is \code{TRUE}.
#'
#'
#' @return A data frame containing the loadings for each component on the
#' cell-cycle and control genes as well as the difference between the loadings
#' and the bootstrapped statistic for each loading.
#'
bootstrap_diff <- function(xy, xn, nboot=200, bar=TRUE)
{
  res0 <- get_diff(xy, xn)
  val <- min(ncol(xn), nrow(xn))
  diff_load_boot <- matrix(NA, nrow=val, ncol=nboot)
  cat("Bootstrapping...")
  if (bar == TRUE){
    pb <- utils::txtProgressBar(min = 1, max = nboot, style = 3)
  }
  for (i in 1 : nboot)
  {
    if (bar == TRUE){
      utils::setTxtProgressBar(pb, i)
    }
    res <- get_diff(xy[sample(1 : nrow(xy), nrow(xy), replace=TRUE), ],
                    xn[sample(1 : nrow(xn), nrow(xn), replace=TRUE), ])
    diff_load_boot[, i] <- res$diff_load
  }
  if(bar == TRUE){
    close(pb)
  }
  sd1 <- apply(diff_load_boot, 1, stats::sd)
  res0$t_load_boot <- res0$diff_load / sd1

  return(res0)
}


#' Calculates the average load difference between the cell-cycle genes and
#' control genes for each component.
#'
#' This is an interal function for use by "bootstrap_diff" only.
#'
#' @param xy The data for the genes which are annotated to the cell-cycle, i.e.
#' those genes for which "if_cc" is \code{TRUE}.
#' @param xn The data for the genes which are not annotated to the cell-cycle,
#' control genes, genes for which "if_cc" is \code{FALSE}
#'
#' @return A data frame containing the loadings for each component on the
#' cell-cycle and control genes.

get_diff <- function(xy, xn)
{
  xn_pca <- stats::prcomp(xn, scale.=FALSE)
  xn_proj <- xn %*% xn_pca$rotation
  xy_proj <- xy %*% xn_pca$rotation

  xn_load <- sqrt(colMeans(xn_proj ^ 2))
  xy_load <- sqrt(colMeans(xy_proj ^ 2))

  return(data.frame(xn_load=xn_load, xy_load=xy_load,
                    diff_load=xy_load-xn_load))
}
