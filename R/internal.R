#' Create cellcount table for calculating percentage parent
#' @param x A data frame containing the Timer Angle variable with the column name Angle.
#' @param percentTimer A logical value for whether Percent Timer data is produced. Default is FALSE and produces Percent Parent data.
#' @return A data frame for the cell number of each Tocky Locus. NA is returned for Timer negative cells.
#' @examples
#' \dontrun{
#' cellcounttable <- cellcount_table(x, percentTimer = FALSE)
#'}
#' @keywords internal
#' @importClassesFrom TockyPrep TockyPrepData

cellcount_table <- function(x, percentTimer = TRUE) {
  sampletable <- x@sampledef$sampledef
  data <- x@Data
  cell_counts <- x@cell_counts
  if(percentTimer){
      data <- data[!is.na(data$Angle),]
  }
  
  locus_levels <- c("New", "NPt", "Persistent", "PAt", "Arrested")
  data$Locus <- factor(data$Locus, levels = locus_levels)
  counts <- table(factor(data$file), data$Locus)
  counts_df <- as.data.frame.matrix(counts)
  counts_df <- cbind(counts_df, data.frame(sample = rownames(counts_df)))
  counts_df <- merge(counts_df, cell_counts, by = 'sample')

  if(percentTimer){
      percentages<- counts_df[locus_levels] / rowSums(counts_df[locus_levels]) * 100
      
      }else{
          percentages <- counts_df[locus_levels] / counts_df$cell_count * 100
      }
  

  percentages[['file']] <- counts_df[['sample']]
  out <- merge(percentages, sampletable, by = "file")

  return(out)
}






#' Calculate mean percentage for each locus using cellcount_table
#' @param x A TockyPrepData or the output of cellcount_table
#' @param sampledef a sample_definition data frame
#' @return A summary table as a data frame.
#' @examples
#' \dontrun{
#' p <- locus_mean(x)
#'}
#' @keywords internal
#' @importFrom stats sd
#' @importClassesFrom TockyPrep TockyPrepData

locus_mean <- function(x, sampledef = sampledef) {
    
    if(inherits(x, "TockyPrepData")){
        cellcounttable <- x@Stats$cellcounttable
        if(length(x@sampledef$sampledef) == 0){
            stop("Apply sample_definition to your TockyPrepData. \n")
        }
        sampledef <- x@sampledef$sampledef
        
    }

    if(inherits(x, "data.frame")){
        cellcounttable <- x
        if(is.null(sampledef)){
            stop("A data frame for sample definition is required. \n")
        }

    }
    group_column <- 'group'
    locus_levels <- c("New", "NPt", "Persistent", "PAt", "Arrested")
    data <- data.frame(
      Percent = c(t(cellcounttable[, locus_levels])),
      Locus = rep(locus_levels, nrow(cellcounttable)),
      File = rep(cellcounttable$file, each = 5),
      Group = rep(cellcounttable[,group_column], each = 5)
    )
    
    data$Locus <- factor(data$Locus, levels = locus_levels)

    group <- as.factor(sampledef$group)

    df <- data.frame(Locus = character(),
                   Group = character(),
                   Percent = numeric(),
                   sd = numeric(),
                   stringsAsFactors = FALSE)
  
    for(i in unique(data$Locus)) {
        for(j in unique(group)) {
            subset_data <- data[(data$Locus == i) & (data$Group == j), ]
            if (nrow(subset_data) > 0) {
                df <- rbind(df, data.frame(Locus = i,
                                    Group = j,
                                    Percent = mean(subset_data$Percent),
                                    sd = sd(subset_data$Percent)))
      }
    }
  }
  return(df)
}






