#' Identifies genes annotated to the cell-cycle
#'
#' Determines which of the genes contained in the dataset are annotated ti the
#' cell-cycle. This is a preprocessing function for ccRemover. Genes can be
#' either mouse or human and either official gene symbols, Ensembl, Entrez or
#' Unigene IDs.
#'
#' @param gene_names A vector containing the gene names for the dataset.
#' @param species The species which the gene names are from. Either
#' \code{"human"} or \code{"mouse"}.
#' @param name_type The type of gene name considered either, Ensembl gene IDS
#' (\code{"ensembl"}), offical gene symbols (\code{symbol}), Entrez gene IDS
#' (\code{"entrez"}), or Unigene IDS (\code{unigene}).
#'
#' @return A vector containg the indices of genes which are annotated to the cell-cycle
#'
#' @examples
#' set.seed(10)
#' # Load in example data
#' data(t.cell_data)
#' head(t.cell_data[,1:5])
#' # Center example data
#' t_cell_data_cen <- t(scale(t(t.cell_data), center=TRUE, scale=FALSE))
#' # Extract gene names
#' gene_names <- rownames(t_cell_data_cen)
#' # Determine which genes are annotated to the cell-cycle
#' cell_cycle_gene_indices <- gene_indexer(gene_names = gene_names,
#' species = "mouse", name_type = "symbol")
#' # Create "if_cc" vector
#' if_cc <- rep(FALSE,nrow(t_cell_data_cen))
#' if_cc[cell_cycle_gene_indices] <- TRUE
#'
#' # Can allow the function to automatically detect the name type
#' cell_cycle_gene_indices <- gene_indexer(gene_names = gene_names,
#' species = NULL, name_type = NULL)
#'
gene_indexer <- function(gene_names, species=NULL,
                         name_type=NULL){
  possible_species <- c("human", "mouse")
  possible_names <- c("ensembl", "symbol", "entrez",
                      "unigene")
  if (!is.null(species)){
    if (!species %in% possible_species){
      cat("Invalid species input. Switching to NULL")
      species <- NULL
    }
  }
  if (!is.null(name_type)) {
    if (!name_type %in% possible_names){
      cat("Invalid name type input. Switching to NULL")
      name_type <- NULL
    }
  }
  if (is.null(species) & is.null(name_type)){
    cat("No species or name format input.", fill=TRUE)
    cat("Checking to see if match can be found:", fill=TRUE)
    match_vec <- rep(NA, 8)
    match_vec[1] <- sum(gene_names %in% human_cell_cycle_genes[,1])
    match_vec[2] <- sum(gene_names %in% human_cell_cycle_genes[,2])
    match_vec[3] <- sum(gene_names %in% human_cell_cycle_genes[,3])
    match_vec[4] <- sum(gene_names %in% human_cell_cycle_genes[,4])
    match_vec[5] <- sum(gene_names %in% mouse_cell_cycle_genes[,1])
    match_vec[6] <- sum(gene_names %in% mouse_cell_cycle_genes[,2])
    match_vec[7] <- sum(gene_names %in% mouse_cell_cycle_genes[,3])
    match_vec[8] <- sum(gene_names %in% mouse_cell_cycle_genes[,4])
    if (max(match_vec) == 0){
      stop("No matches found. Please convert gene names to a suitable type")
    } else if (which.max(match_vec) <= 4){
      cat("Human genes detected:", fill=TRUE)
      guess <- possible_names[which.max(match_vec)]
      cat("Best guess is ", guess , " IDs",
          fill = TRUE)
      matches <- match_vec[which.max(match_vec)]
      cat(matches, " matches out of a possible ",
          length(gene_names), fill = TRUE)
      species <- "human"
      name_type <- possible_names[which.max(match_vec)]
    } else if (which.max(match_vec) >= 5){
      cat("Mouse genes detected:" , fill = TRUE)
      cat("Best guess is ", possible_names[which.max(match_vec)-4], " IDs",
          fill = TRUE)
      cat(match_vec[which.max(match_vec)]," matches out of a possible ",
          length(gene_names) , fill = TRUE)
      species <- "mouse"
      name_type <- possible_names[which.max(match_vec)-4]
    }
  } else if (!is.null(species) & is.null(name_type)) {
    cat("No name format input.", fill=TRUE)
    cat("Checking to see if match can be found:" , fill = TRUE)
    if (species == "human"){
      match_vec <- rep(NA, 4)
      match_vec[1] <- sum(gene_names %in% human_cell_cycle_genes[,1])
      match_vec[2] <- sum(gene_names %in% human_cell_cycle_genes[,2])
      match_vec[3] <- sum(gene_names %in% human_cell_cycle_genes[,3])
      match_vec[4] <- sum(gene_names %in% human_cell_cycle_genes[,4])
      if (max(match_vec) == 0){
        stop("No matches found. Please convert gene names to a suitable type")
      } else if ( max(match_vec) > 0){
        cat("Best guess is ", possible_names[which.max(match_vec)], " IDs",
            fill = TRUE)
        cat(match_vec[which.max(match_vec)]," matches out of a possible ",
            length(gene_names), fill = TRUE)
        name_type <- possible_names[which.max(match_vec)]
      }
    } else if (species == "mouse"){
      match_vec <- rep(NA, 4)
      match_vec[1] <- sum(gene_names %in% mouse_cell_cycle_genes[,1])
      match_vec[2] <- sum(gene_names %in% mouse_cell_cycle_genes[,2])
      match_vec[3] <- sum(gene_names %in% mouse_cell_cycle_genes[,3])
      match_vec[4] <- sum(gene_names %in% mouse_cell_cycle_genes[,4])
      if (max(match_vec) == 0){
        cat("No matches found. Please convert gene names to a suitable type",
            fill = TRUE)
      } else if ( max(match_vec) > 0){
      cat("Best guess is ", possible_names[which.max(match_vec)], " IDs",
          fill = TRUE)
      cat(match_vec[which.max(match_vec)]," matches out of a possible ",
          length(gene_names), fill = TRUE)
      name_type <- possible_names[which.max(match_vec)]
      }
    }
  } else if (is.null(species) &!is.null(name_type)){
    cat("No species format input.", fill = TRUE)
    cat("Checking to see if match can be found:", fill = TRUE)
    match_vec <- rep(NA,2)
    match_vec[1] <- sum(gene_names %in%
                          human_cell_cycle_genes[,which(possible_names ==
                                                          name_type)])
    match_vec[2] <- sum(gene_names %in%
                          mouse_cell_cycle_genes[,which(possible_names ==
                                                          name_type)])
    if (max(match_vec) == 0){
      stop("No matches found. Please convert gene names to a suitable type")
    } else if (which.max(match_vec) == 1){
      cat("Human genes detected:", fill = TRUE)
      cat(match_vec[1]," matches out of a possible ",
          length(gene_names), fill = TRUE)
      species <- "human"
    } else if (which.max(match_vec) == 2){
      cat("Mouse genes detected:", fill = TRUE)
      cat(match_vec[2]," matches out of a possible ",
          length(gene_names), fill = TRUE)
      species <- "mouse"
    }
  }
  cc_gene_indices <- rep(NA, 1)
  if (species == "human"){
    human_cell_cycle_genes <- NULL
    utils::data("HScc_genes", envir = environment())
    if (name_type == "ensembl"){
      cc_gene_indices <- which(gene_names %in% human_cell_cycle_genes[,1])
    } else if (name_type == "symbol"){
      cc_gene_indices <- which(gene_names %in% human_cell_cycle_genes[,2])
    } else if (name_type == "entrez"){
      cc_gene_indices <- which(gene_names %in% human_cell_cycle_genes[,3])
    } else if (name_type == "unigene"){
      cc_gene_indices <- which(gene_names %in% human_cell_cycle_genes[,4])
    }
  } else if (species == "mouse") {
    mouse_cell_cycle_genes <- NULL
    utils::data("MMcc_genes", envir = environment())
    if (name_type == "ensembl"){
      cc_gene_indices <- which(gene_names %in% mouse_cell_cycle_genes[,1])
    } else if (name_type == "symbol"){
      cc_gene_indices <- which(gene_names %in% mouse_cell_cycle_genes[,2])
    } else if (name_type == "entrez") {
      cc_gene_indices <- which(gene_names %in% mouse_cell_cycle_genes[,3])
    } else if (name_type == "unigene") {
      cc_gene_indices <- which(gene_names %in% mouse_cell_cycle_genes[,4])
    }
  }
  return(cc_gene_indices)
}
