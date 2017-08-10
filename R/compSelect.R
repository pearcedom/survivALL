#' Given a dataframe of phenotypic information, use a variable (i.e. a single 
#' column) to define a patient subset of given proportion
#' @param pheno a phenotypic dataframe. Sample IDs must be assigned as rownames
#' @param column the name of the column used to define the subset - e.g. "grade"
#' @param values the values within column that you are defining by These must 
#' be categorical - e.g. c("+", "-")  
#' @param props The number of how many of each value in column to be returned. 
#' - e.g. c(50, 50). Note, be careful not to ask for more samples of a 
#' particular value than are available in the dataset
#' @return A dataframe, which is the subset of pheno, with a specified 
#' proportion of each value found in column
#' @examples
#' library(survivALL)
#' data(nki_subset)
#' library(Biobase)
#' pheno <- pData(nki_subset)
#'
#' compSelect(pheno, "grade", values = c(1, 2, 3), props = c(10, 10, 5))
#'
#' #to manufacture composition from a continuous measure, first translate into a
#' #categorical equivalent, e.g.
#' age <- pheno$age
#' pheno$age_group <- ifelse(age < 40, "<40", ifelse(age < 50, "40-50", ">=50"))
#' compSelect(pheno, 
#'     "age_group", 
#'     values = c("<40", "40-50", ">=50"), 
#'     props = c(2, 5, 10))
#' @export
compSelect <- function(pheno, column, values, props) {
    do.call(rbind, 
            lapply(1:length(values), function(x) {
                       samples <- 
                           sample(row.names(pheno)[which(pheno[[column]] == 
                                                             values[x])],
                                  props[x])
                       pheno[samples,]
    })) 
}
