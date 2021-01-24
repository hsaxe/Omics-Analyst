#' Takes an OmicsAnalyst formatted data.frame, replaces any missing value(s) with zeros, and removes variables (rows) with constant values.
#' @param x a data.frame
#' @export
pre_process <- function(x){

  ifelse(any(is.na(x)), print(paste(sum(is.na(x)), "missing value(s) replaced with zeros")), print("Data has no missing values"))

  if(any(is.na(x))){x[is.na(x)] <- 0}

  rawT <- data.table::transpose(x, make.names = 1)

  rem_con <- janitor::remove_constant(rawT, quiet = T)

  rem_con$Sample <- as.character(colnames(x)[2:length(colnames(x))])

  t <- data.table::transpose(rem_con, keep.names = "Gene.ID")

  colnames(t) = colnames(x)

  print(paste(nrow(x)-nrow(t), "features with constant values removed"))

  return(t)
}

#' Takes an OmicsAnalyst formatted data.frame and removes varbiables(rows) with IQR below a certain threshold
#' @param x a data.frame
#' @param Threshold a value indicating the minimum IQR to be included in the output. Varibles having IQR below this will be excluded from the output
#' @export
iqr_filter <- function(x, Threshold = 0.5){

  x$iqrs = c("NA", matrixStats::rowIQRs(as.matrix(sapply(x[-1,-1], as.numeric))))

  filtered <- x[x$iqrs >= Threshold,]

  filtered[1,1] <- x[1,1]

  final <- filtered[,-ncol(filtered)]

  print(paste("Removed", ((nrow(x)-1) - nrow(filtered)-1), "features with IQR below", Threshold))

  return(filtered)
}

#' Takes an OmicsAnalyst formatted data.frame and converts it to long format for plotting descriptive statistics and individual features
#' @param x a data.frame
#' @export

format_to_plot <- function(x){

  a <- reshape2::melt(x[-1,-ncol(x)], id.vars = 1, variable.name = "Sample", value.name = "value")

  a$Class = gsub("\\.\\d$", "", a$Sample)

  a$value = as.numeric(a$value)

  return(a)

}


#' Takes a plotting data.frame (typically long format) and returns the data.frame/data.table with a new column grouped by "group", transformed by "fun", and named by "name" arguments.
#' @param x a data.frame
#' @param group column name of groups to group by
#' @param variable column name of variable that will be grouped by "groups"
#' @param fun function to apply to grouped variable
#' @param name desired name for new grouped variable produced by "fun"
#' @export
group_stats <- function(x, group = "", variable = "", fun, name = ""){

  d <- data.table::data.table(x)

  a <- d[, fun(get(variable)), by = group]

  colnames(a)[ncol(a)] <- name

  m <- merge(d, a, by = group)

  return(m)
}
