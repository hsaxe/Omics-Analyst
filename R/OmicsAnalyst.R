#' Takes a data.frame, replaces any missing value(s) with zeros, optionally performs counts per million (CPM) normalization, calculates the row-wise statistic provided by FilterFun, and filters the data with above the FilterThreshold or keeps only the top n provided by RankThreshold. Also, the output of expression filter is a two element list. 'final' is a data frame showing the filtered and if specified, transformed, data and a column with whatever statistic the user specified to filter it by. This is to visualize what features had what filter/rank statistic. 'dat' is the raw data (not CPM normalized). If specified, 'dat' will be a DGEList.
#' @param dat A data.frame
#' @param CPM Logical. Whether or not to normalize columns by CPM
#' @param FilterFUN Row-wise function to use for filtering (mean, max, median, IQR, etc..).
#' @param FilterThreshold Threshold value. Features with FilterFUN output higher than this value will be kept.
#' @param RankThreshold Threshold value. Features will be ranked with FilterFUN and the top n values provided by RankThreshold will be kept.
#' @param DGEList Logical. Is input data a DGEList? If TRUE, input data is handled as a list with a data.matrix/frame in it and not a data.frame alone. Also, CPM normalization will be performed if DGEList is set to TRUE
#' @param CPH Logical. If there are very few mapped reads, your depth may not be in millions, but in hundreds or thousands. Therefore CPM normalization will skew expression very high. CPH normalization, or normalizing by hundreds factor, may provide better results.
#' @import dplyr
#' @importFrom dplyr filter
#' @importFrom rlang :=
#' @importFrom stats prcomp reorder
#' @importFrom utils head
#' @export
expression_filter <- function(dat,
                              DGEList = FALSE,
                              CPM = TRUE,
                              CPH = FALSE,
                              FilterFUN = mean,
                              FilterThreshold = NULL,
                              RankThreshold = NULL){

  message('Make sure no numeric identifiers are in data as this will negatively impact filtering')

  iqrs <- . <- ID <- Group <- NULL

  if(DGEList == T) {

    a = as.data.frame(dat$counts)

  } else {

    a <- as.data.frame(dat)

  }

  if(any(sapply(a, is.character)) == T) {
    stop('Characters detected in at least one column. Only numeric values allowed. The ID column (GeneID, proteinID, etc.) must be moved to rownames. Also make sure the data is of type "numeric."')
  }

  if(any(is.na(a))){

    print(paste(sum(is.na(a)), "missing value(s) replaced with zeros"))

  } else {

   print("Data has no missing values")

  }

  if(any(is.na(a))){a[is.na(a)] <- 0}

  if(CPH == T) {

    cpm = a %>%
      # column_to_rownames('GeneID') %>%
      apply(2, FUN = function(x){
        x/(sum(x)/100)
      }) %>%
      as.data.frame()

  } else {

    if(DGEList == T) {

      cpm = dat %>% edgeR::cpm() %>% as.data.frame()

    } else {

      if(CPM == T) {cpm = a %>% edgeR::cpm() %>% as.data.frame()}

    }

  }





  FilterFUN = enquo(FilterFUN)

  if(is.null(FilterThreshold) & is.null(RankThreshold) | !is.null(FilterThreshold) & !is.null(RankThreshold)) {
    stop('Please provide EITHER FilterThreshold OR RankThreshold')
  }

  if(!is.null(FilterThreshold)) {

    final = cpm %>%
      mutate(!!FilterFUN := apply(., 1, !!FilterFUN)) %>%
      filter(!!FilterFUN >= FilterThreshold)

    filter = list(final = final)

    # dat = append(dat, filter)

    print(paste0("Removed ", (nrow(a)) - (nrow(final)), " features based on ", as_label(FilterFUN), " threshold of ", FilterThreshold, ', ', nrow(final), ' remaining'))


    a = a %>%
      tibble::rownames_to_column(var = 'ID') %>%
      filter(ID %in% rownames(final)) %>%
      tibble::column_to_rownames(var = 'ID')

    if(DGEList == T) {

      dat$counts = as.matrix(a)

      return(dat)

    } else {

      return(a)

    }



  } else {

    if(!is.null(RankThreshold)) {

      final = cpm %>%
        mutate(!!FilterFUN := apply(., 1, !!FilterFUN)) %>%
        arrange(desc(!!FilterFUN)) %>%
        head(n = RankThreshold)

      filtered = list(final = final,
                      dat = dat)

      # dat = append(dat, filter)

      print(paste0("Removed ", (nrow(a)) - (nrow(final)), " features based on ", as_label(FilterFUN), " rank threshold of ", RankThreshold, ', ', nrow(final), ' remaining'))


      a = a %>%
        tibble::rownames_to_column(var = 'ID') %>%
        filter(ID %in% rownames(final)) %>%
        tibble::column_to_rownames(var = 'ID')

      if(DGEList == T) {

        filtered$dat$counts = as.matrix(a)

        return(filtered)

      } else {

        return(a)

      }

    }

  }

}


#' Multiple options for plotting descriptive statistics from PCA. Takes a numeric data frame as input.
#' @param dat A data.frame
#' @param metadata Optional data frame with metadata for filling, coloring, or plotting PCA results against.
#' @param join_by_name Name of column in metadata that matches rownames/colnames of data. Used to join metadata with data. This is usually sample names. For example, sample1.1, sample1.2, etc.
#' @param plotting_factors_name Name for factor you want to plot. Is made by removing two characters from the end of 'join_by_name' variable. For example, sample1.1, sample1.2 and sample2.1, sample2.2 in this variable become sample1, sample1 and sample2, sample2. Used for coloring, filling, and grouping of plot. Default name is Group
#' @param plotting_factors_in Location of plotting factors. Options are either 'col_names' or 'row_names'.
#' @param x What do you want to plot on the x-axis? Default is 'PC1'
#' @param y What do you want to plot on the y-axis? Default is 'PC2'
#' @param scale Logical. Scale data? Default is TRUE
#' @param center Logical. Center data? Default is TRUE
#' @param color Which variable to color by? Default is 'Group'.
#' @param fill Which variable to fill by? Default is 'Group'.
#' @param plot_type One of three options: '2D', 'boxplot', or 'scatter'. Default is '2D'.
#' @param summarise_for_scatter Logical. Plotting factors can sometimes contain psuedoreplication which inflates the p-value of this scatterplot. This option will summarize the plotting factors by mean, removing psuedoreplication for a more realistic p-value. Default is TRUE.
#' @param points_label Which variable to label points by? Default is 'NULL' in which case no labeling will occur.
#' @param sep Desired separator for sample names in data to match that of metadata
#' @import dplyr
#' @import tibble
#' @import ggpubr
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom rlang :=
#' @importFrom stats prcomp reorder
#' @importFrom utils head
#' @export
plot_pca = function(dat,
                    metadata = NULL,
                    join_by_name = 'Sample',
                    plotting_factors_name = 'Group',
                    plotting_factors_in = 'col_names',
                    x = 'PC1',
                    y = 'PC2',
                    scale = T,
                    center = T,
                    color = 'Group',
                    fill = 'Group',
                    plot_type = '2D',
                    points_label = NULL,
                    summarise_for_scatter = T,
                    sep = '-') {

  . <- ID <- Group <- NULL


  message('Make sure no numeric identifiers are in data as this will drastically impact PCA')

  if(any(sapply(dat, is.character)) == T) {
    stop('Characters detected in at least one column. Only numeric values allowed. The ID column (GeneID, proteinID, etc.) must be moved to rownames. Also make sure the data is of type "numeric."')
  }

  plotting_factors_name = enquo(plotting_factors_name)

  # join_by_name = enquo(join_by_name)

  if(plotting_factors_in == 'col_names') {
    dat = t(dat) %>%
      as.data.frame() %>%
      mutate(across(everything(), as.numeric)) %>%
      scale(scale = scale, center = center) %>%
      as.matrix()
  } else {
    dat = dat %>%
      # as.data.frame() %>%
      column_to_rownames(var = rlang::quo_name(plotting_factors_name)) %>%
      # select_if(is.character) %>%
      mutate(across(everything(), as.numeric)) %>%
      scale(scale = scale, center = center) %>%
      as.matrix()
  }


  dat = prcomp(dat, scale = F, center = F)

  loadings = dat$rotation

  scores = dat$x

  eigs = dat$sdev^2

  var_exp = paste('(', formatC((eigs/sum(eigs))*100, format = 'f', digits = 2), '%', ')', sep = '')

  names(var_exp) = colnames(scores)

  plot_list = list()

  plot_list$loadings = loadings %>%
    as.data.frame() %>%
    # arrange(PC1) %>%
    rownames_to_column(var = 'Variable')

  plot_list$var_exp = var_exp

  plot_list$loadings_plot = ggplot(plot_list$loadings,
                                   aes(PC1, PC2, color = Variable))+
    geom_point()+
    geom_label_repel(aes(label = Variable))+
    ggtitle('Loadings')+
    theme(plot.title = element_text(hjust = 0.5))

  if(is.null(metadata)){

    plot_dat = scores %>%
      data.frame() %>%
      rownames_to_column(var = 'var') %>%
      mutate(var = gsub('X', '', var) %>%
               gsub('\\.|-', sep, .)) %>%
      rename(!!join_by_name := var) %>%
      mutate(!!plotting_factors_name := gsub('\\.\\d+$', '', get(join_by_name)))

    plot_list$plot_dat <- plot_dat


  } else {

    plot_dat = scores %>%
      data.frame() %>%
      rownames_to_column(var = 'var') %>%
      mutate(var = gsub('X', '', var) %>%
                  gsub('\\.|-', sep, .)) %>%
      rename(!!join_by_name := var) %>%
      mutate(!!plotting_factors_name := gsub('\\.\\d+$', '', get(join_by_name))) %>%
      left_join(metadata, by = join_by_name)
      # left_join(metadata)

    plot_list$plot_dat <- plot_dat

  }



  if(plot_type == '2D') {
    if(is.null(points_label)) {

      p = ggplot(plot_dat, aes(get(x), get(y), color = get(color), group = get(fill), fill = get(fill)))+
        geom_point()+
        stat_ellipse(geom = 'polygon', alpha = 0.5, level = 0.65)+
        # ggforce::geom_mark_ellipse(aes(fill = get(fill), label = !!plotting_factors_name))+
        labs(x = paste(x, var_exp[x]), y = paste(y, var_exp[y]), fill = color, color = color)

      plot_list$scores_plot <- p

      return(plot_list)

    } else {

      p = ggplot(plot_dat, aes(get(x), get(y), color = get(color), group = get(fill), fill = get(fill)))+
        geom_point()+
        stat_ellipse(geom = 'polygon', alpha = 0.5, level = 0.65)+
        # ggforce::geom_mark_ellipse(aes(fill = get(fill), label = !!plotting_factors_name))+
        geom_text(aes(label = get(points_label)), color = 'black', size = 2.5)+
        labs(x = paste(x, var_exp[x]), y = paste(y, var_exp[y]), fill = color, color = color)

      plot_list$scores_plot <- p

      return(plot_list)

    }

  }

  if(plot_type == 'boxplot') {

    p = ggplot(plot_dat, aes(reorder(get(x), get(y)), get(y), fill = get(fill)))+
      geom_boxplot()+
      labs(x = plotting_factors_name, y = paste(y, var_exp[y]), fill = fill)

    plot_list$scores_boxplot <- p

    return(plot_list)
  }

  if(plot_type == 'scatter') {
    if(summarise_for_scatter == T) {

      plot_dat = plot_list$plot_dat %>%
        group_by(!!plotting_factors_name) %>%
        summarise_if(is.numeric, mean)

      plot_list$plot_dat <- plot_dat

      p = ggscatter(plot_dat, x = x, y = y, add = 'reg.line', color = color,
                    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                    cor.coeff.args = list(method = "pearson", label.x = 0, label.y = max(plot_dat[y])*1.1),
                    cor.coef.size = 5)+
        # stat_cor(label.x = 0, label.y = max(plot_dat[y])*1.1)+
        stat_smooth(method = 'lm')+
        labs(x = paste(x, var_exp[x]))

      plot_list$scatter_plot <- p

      return(plot_list)

    } else {
      p = ggscatter(plot_dat, x = x, y = y, add = 'reg.line', color = get(color))+
        stat_cor(label.x = 0, label.y = max(plot_dat[y])*1.1)+
        stat_smooth(method = 'lm')+
        labs(x = paste(x, var_exp[x]))

      plot_list$scatter_plot <- p

      return(plot_list)
    }

  }

}


