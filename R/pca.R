#' Create variable groups for PCA analysis
#' @param data A data frame containing survey data
#' @param group_definitions A list where names are group categories and elements are character vectors of variable names
#' @return A list containing the grouped variables that exist in the dataset
#' @export
create_variable_groups <- function(data, group_definitions) {
  existing_cols <- colnames(data)
  grouped_vars <- lapply(group_definitions, function(vars) {
    vars[vars %in% existing_cols]
  })
  grouped_vars[sapply(grouped_vars, length) > 0]
}

#' Perform PCA analysis with contribution calculation
#' @param data A data frame containing survey data
#' @param vars Character vector of variables to include
#' @param scale Logical indicating whether to scale variables
#' @param n_components Number of components to retain
#' @return List containing PCA results and contributions
#' @export
perform_pca_analysis <- function(data, vars, scale = TRUE, n_components = NULL) {
  if (is.null(n_components)) {
    n_components <- min(5, ncol(data[vars]))
  }
  
  # Perform PCA
  pca_result <- prcomp(data[vars], scale. = scale)
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
  
  # Calculate variable contributions
  loadings <- abs(pca_result$rotation[, 1:n_components])
  contributions <- sweep(loadings^2, 2, var_explained[1:n_components], "*")
  
  list(
    pca = pca_result,
    var_explained = var_explained,
    contributions = contributions,
    total_contribution = rowSums(contributions)
  )
}

#' Calculate group contributions to PCA components
#' @param pca_results Results from perform_pca_analysis
#' @param variable_groups List of variable groups from create_variable_groups
#' @return Data frame of group contributions
#' @export
calculate_group_contributions <- function(pca_results, variable_groups) {
  group_contribs <- lapply(names(variable_groups), function(group) {
    vars <- variable_groups[[group]]
    contrib <- colSums(pca_results$contributions[vars, , drop = FALSE])
    data.frame(
      Group = group,
      Component = paste0("PC", 1:length(contrib)),
      Contribution = contrib
    )
  })
  do.call(rbind, group_contribs)
}

#' Generate correlation analysis for PCA variables
#' @param data A data frame containing survey data
#' @param vars Character vector of variables to analyze
#' @param threshold Correlation threshold for highlighting strong correlations
#' @return List containing correlation matrix and strong correlations
#' @export
analyze_correlations <- function(data, vars, threshold = 0.5) {
  cor_matrix <- cor(data[vars])
  
  # Find strong correlations
  strong_cors <- which(abs(cor_matrix) > threshold & cor_matrix != 1, arr.ind = TRUE)
  strong_cors_df <- data.frame(
    Var1 = rownames(cor_matrix)[strong_cors[,1]],
    Var2 = colnames(cor_matrix)[strong_cors[,2]],
    Correlation = cor_matrix[strong_cors]
  )
  
  list(
    correlation_matrix = cor_matrix,
    strong_correlations = strong_cors_df[order(-abs(strong_cors_df$Correlation)),]
  )
}

#' Select optimal variables for PCA based on contributions and correlations
#' @param pca_results Results from perform_pca_analysis
#' @param correlation_results Results from analyze_correlations
#' @param min_contribution Minimum total contribution to retain variable
#' @param max_correlation Maximum correlation between retained variables
#' @return Character vector of selected variables
#' @export
select_optimal_variables <- function(pca_results, correlation_results, 
                                   min_contribution = 5, max_correlation = 0.8) {
  # Order variables by contribution
  vars_ordered <- names(sort(pca_results$total_contribution, decreasing = TRUE))
  
  # Initialize selected variables with highest contributor
  selected_vars <- vars_ordered[1]
  
  # Evaluate remaining variables
  for (var in vars_ordered[-1]) {
    # Check correlations with already selected variables
    cors <- correlation_results$correlation_matrix[var, selected_vars, drop = FALSE]
    if (all(abs(cors) < max_correlation) && 
        pca_results$total_contribution[var] >= min_contribution) {
      selected_vars <- c(selected_vars, var)
    }
  }
  
  selected_vars
}

#' Create a summary report of PCA analysis
#' @param pca_results Results from perform_pca_analysis
#' @param group_contributions Results from calculate_group_contributions
#' @param selected_variables Results from select_optimal_variables
#' @return List containing summary statistics and recommendations
#' @export
create_pca_report <- function(pca_results, group_contributions, selected_variables) {
  n_significant_pcs <- sum(pca_results$var_explained > 5)
  
  list(
    n_components_retained = n_significant_pcs,
    variance_explained = pca_results$var_explained[1:n_significant_pcs],
    cumulative_variance = cumsum(pca_results$var_explained)[1:n_significant_pcs],
    group_importance = aggregate(
      Contribution ~ Group, 
      group_contributions, 
      sum
    ),
    selected_variables = selected_variables,
    recommendations = list(
      n_components = n_significant_pcs,
      key_groups = group_contributions$Group[
        order(group_contributions$Contribution, decreasing = TRUE)[1:3]
      ]
    )
  )
}
