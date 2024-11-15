#' Classify Angle Data into Tocky Loci
#'
#' This function classifies each value in a vector of angles into specific loci categories
#' based on the angle's value. It categorizes angles into 'New', 'NPt', 'Persistent',
#' 'PAt', 'Arrested', or `NA` for Timer negative cells.
#'
#' @param angle A `A numeric vector` of angles to classify.
#' @return Tocky Locus data as a character vector.
#' @examples
#' \dontrun{
#' angle <- c(0, 29, 45, 65, 90)
#' locus <- angle_to_locus(angle)
#'}
#' @keywords internal

angle_to_locus <- function(angle){
    new  <-  angle == 0 & !is.na(angle)
    npt  <-  angle > 0 & angle < 30 & !is.na(angle)
    pers  <-  angle >= 30 & angle < 60 & !is.na(angle)
    pat  <-  angle >= 60 & angle < 90 & !is.na(angle)
    arr  <-  angle == 90 & !is.na(angle)
    output  <-  ifelse(new,"New", ifelse(npt, "NPt", ifelse(pers, "Persistent", ifelse(pat, "PAt", ifelse(arr, "Arrested", NA)))))
    return(output)
}

#' Calculate Tocky Locus using Timer Angle
#' @param x A data.frame object or a TockyPrepData object
#' @return Input data frame object including Tocky Locus data.
#' @export
#' @examples
#' \dontrun{
#' x <- TockyLocus(x)
#'}
#'
#' @importClassesFrom TockyPrep TockyPrepData

TockyLocus <- function(x){
    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData. \n")
        
    }
    
    if(inherits(x, "TockyPrepData")){
        data <- x@Data
        sampledef  <- x@sampledef
        
    }
        
    locus  <- angle_to_locus(data$Angle)
    
    if(is.null(data$Locus)){
        locus  <- angle_to_locus(data$Angle)
        data <- cbind(data.frame(Locus=locus), data)
        x@Data <- data
    }else{

         data <- x@Data

    }
    
    cellcounttable <- cellcount_table(x, percentTimer = FALSE)

    sampledef_df <- sampledef$sampledef
    locus_levels <- c("New", "NPt", "Persistent", "PAt", "Arrested")
    
    df <- data.frame(
        Percent = c(t(cellcounttable[, locus_levels])),
        Locus = rep(c('New', 'NPt', 'Persistent', 'PAt', 'Arrested'), nrow(cellcounttable)),
        File = rep(cellcounttable$file, each = 5),
        Group = rep(cellcounttable[, 'group'], each = 5)
    )
    
    df$Locus <- factor(df$Locus, levels = locus_levels)
    df_summary <- locus_mean(cellcounttable, sampledef = sampledef_df)

    data <-  data[!is.na(data$Angle),]
    cellcounttable <- cellcount_table(x, percentTimer = TRUE)
    df_timer <- data.frame(
        Percent = c(t(cellcounttable[, locus_levels])),
        Locus = rep(c('New', 'NPt', 'Persistent', 'PAt', 'Arrested'), nrow(cellcounttable)),
        File = rep(cellcounttable$file, each = 5),
        Group = rep(cellcounttable[, 'group'], each = 5)
    )
    
    df_timer$Locus <- factor(df_timer$Locus, levels = locus_levels)
    df_timer_summary <- locus_mean(cellcounttable, sampledef = sampledef_df)

    
    TockyLocusStats = list(
        percentTimer = list(data = df_timer, summary =df_timer_summary),
        percentParent = list(data = df, summary =df_summary)
    )
    
    x@Tocky[['TockyLocusStats']] <- TockyLocusStats
    
    return(x)
    
}



#' Produce scatter plots of percentages of cells in each Tocky Locus.
#' @param x A TockyPrepData object
#' @param percentTimer A logical value for whether Percent Timer data is produced. Default is FALSE and produces Percent Parent data.
#' @param p_adjust_method A method for p-value adjustment in statistical tests.
#' @param group_order The order of groups (optional).
#' @param group_by A logical value for whether different groups are plotted in different panels.
#' @param stats A logical value for whether to produce statistical outputs. This is effective only for two-group analysis.
#' @param ylim (Optional) the range of y values to be displayed.
#' @param locus_colours (optional) to choose colours for Tocky Loci.
#' @param group_colors (optional) to choose colours for groups.
#' @param verbose Logical indicating whether to print Tocky Locus stats. Default is `TRUE`.
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' plotTockyLocus(x)
#'}
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang sym
#' @importClassesFrom TockyPrep TockyPrepData

plotTockyLocus <- function(x, percentTimer = FALSE, group_order = NULL, locus_colours = NULL, group_colors = NULL, group_by = TRUE, p_adjust_method = 'fdr', ylim = NULL, stats = TRUE, verbose = TRUE){
    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    
    if(length(x@Data$Locus) == 0){
        stop("Apply TockyLocus to your TockyPrepData. \n")
    }
    
    
    if(!is.null(group_order)){
        group <- as.vector(x@sampledef$sampledef$group)
        group <- factor(group, levels = group_order)
    }else{
        group <- as.factor(x@sampledef$sampledef$group)
    }
    num_groups <- length(levels(group))
    
    if(percentTimer){
        df <- x@Tocky$TockyLocusStats$percentTimer$data
        df_summary <- x@Tocky$TockyLocusStats$percentTimer$summary
        
    }else{
        df <- x@Tocky$TockyLocusStats$percentParent$data
        df_summary <- x@Tocky$TockyLocusStats$percentParent$summary
    }
    
    locus_levels <- c('New', 'NPt', 'Persistent', 'PAt', 'Arrested')
    
    if(is.null(locus_colours)){
        locus_colours <- c('blue', '#2B7BB4', 'purple', '#BF4080', 'red')
        names(locus_colours) <- locus_levels
    }
    if(!group_by & is.null(group_colors)){
        if (num_groups <= 8) {
            group_colors <- brewer.pal(num_groups, "Set1")
        } else {
            group_colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_groups)
        }
    }
    
    
    df$Locus <- factor(df$Locus, levels = locus_levels)
    df_summary$Locus <- factor(df_summary$Locus, levels = locus_levels)
    
    if(!is.null(group_order)){
        df$Group <- factor(df$Group, levels = group_order)
        df_summary$Group <- factor(df_summary$Group, levels = group_order)
        
    }else{
        df$Group <- as.factor(df$Group)
        df_summary$Group <- as.factor(df_summary$Group)
        
    }
    
    if(verbose){
        cat("plotting \n")
    }
    
    if(percentTimer){
        title  <-  "%Timer"
    }else{
        title  <-  "%Parent population"
    }
    
    if(group_by & length(levels(group)) > 1){
        p <- ggplot() +
        geom_jitter(data = df,
        aes(x = !!sym('Locus'), y = !!sym('Percent'), colour = !!sym('Locus')),
        position = position_jitter(width = 0.02)
        ) +
        geom_errorbar(data = df_summary,
        aes(x = !!sym('Locus'), ymin = !!sym('Percent') - !!sym('sd'), ymax = !!sym('Percent') + !!sym('sd'),
        group = !!sym('Group'),colour = !!sym('Locus')
        ), width = 0.1
        ) +
        scale_color_manual(values = locus_colours) +
        facet_grid(cols = vars(!!sym('Group')), scales = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(title) +
        ylab('Percentage (%)') + geom_line(data = df_summary,
        aes(x = !!sym('Locus'), y = !!sym('Percent'), group = !!sym('Group'), colour = 'black'),
        linewidth = 0.6
        )
        
        
    } else {
        p <- ggplot() +
        geom_jitter(data = df,
        aes(x = !!sym('Locus'), y = !!sym('Percent'), colour = !!sym('Group')),
        position = position_jitter(width = 0.02)
        ) +
        geom_errorbar(data = df_summary,
        aes(x = !!sym('Locus'), ymin = !!sym('Percent') - !!sym('sd'), ymax = !!sym('Percent') + !!sym('sd'),
        group = !!sym('Group'), colour = !!sym('Group')), width = 0.1) +
        geom_line(
        data = df_summary,
        aes(x = !!sym('Locus'), y = !!sym('Percent'), group = !!sym('Group'), colour = !!sym('Group')),
        linewidth = 0.6
        ) +
        scale_color_manual(values = group_colors) +  #
        theme_bw() +
        theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle(title) +
        ylab('Percentage (%)')
        
    }
    
    if(!is.null(ylim)){
        p <- p + ylim(ylim)
    }

    return(p)
}




#' Perform Statistical Tests for Tocky Locus Analysis
#'
#' This function performs statistical tests on Tocky Locus data, allowing for different methods and p-value adjustments.
#'
#' @param x A \code{TockyPrepData} object containing Tocky Locus data.
#' @param percentTimer Logical. If \code{TRUE}, the percentages of Timer-positive cells will be used; if \code{FALSE}, percentages of parent cells will be used.
#' @param p_adjust_method Character string specifying the method for p-value adjustment in multiple testing. Default is \code{'BH'} (Benjamini-Hochberg). Other methods available in \code{\link{p.adjust}}, such as \code{'holm'} or \code{'bonferroni'}, can also be used.
#' @param method Character string specifying the statistical test method to use. Options are:
#' \describe{
#'   \item{\code{'Wilcox'}}{Mann-Whitney U test (Wilcoxon rank sum test) without data transformation.}
#'   \item{\code{'ASR'}}{Arcsine Square Root Transformation, followed by a normality test and t-test.}
#'   \item{\code{'Logit'}}{Logit Transformation, followed by a normality test and t-test.}
#' }
#' @param verbose Logical indicating whether to print progress messages and outputs. Default is \code{TRUE}.
#' @return A \code{TockyPrepData} object containing the statistical outputs for Tocky Locus Analysis, stored in \code{x@Tocky$TockyLocusStats}.
#' @export
#' @examples
#' \dontrun{
#' x <- GetStatsTockyLocus(x, method = 'ASR')
#' }
#' @importFrom stats aggregate p.adjust shapiro.test t.test wilcox.test sd
#' @importFrom methods show
#' @importClassesFrom TockyPrep TockyPrepData

GetStatsTockyLocus <- function(x, percentTimer = FALSE, p_adjust_method = 'BH', method = "ASR", verbose = TRUE){
    
    input <- list(percentTimer = percentTimer, p_adjust_method = p_adjust_method, method = method, verbose = verbose)
    
    if(percentTimer){
        df <- x@Tocky$TockyLocusStats$percentTimer$data
        
    }else{
        df <- x@Tocky$TockyLocusStats$percentParent$data
    }
    sampledef <- x@sampledef$sampledef
    group <- as.factor(sampledef$group)
    num_groups <- length(levels(group))
    
    if(!any(c("Wilcox", "ASR", "Logit") %in% method)){
        stop("Use one of the available methods, Wilcox, ASR, or Logit. \n")
    }
    
    if(num_groups !=2){
        stop("Wilcox test is applied to two group comparisons only. \n")
    }
    
    locus_levels <- c('New', 'NPt', 'Persistent', 'PAt', 'Arrested')
    df$Locus <- factor(df$Locus, levels = locus_levels)
    
    temp <- aggregate(Percent ~ Group + Locus, data = df, function(x) c(mean = mean(x), sd = sd(x)))
    
    summary_data <- data.frame(
    Locus = temp$Locus,
    Group = temp$Group,
    Mean_Percentage = temp$Percent[, 'mean'],
    SD_Percentage = temp$Percent[, 'sd']
    )
    summary_data$Locus <- factor(summary_data$Locus, levels = locus_levels)
    summary_data <- summary_data[order(summary_data$Locus), ]
    
    rownames(summary_data) <- paste(summary_data$Locus, summary_data$Group, sep = "_")

    num_groups <- levels(group)
    
    
    if(method == "Wilcox"){
        split_data <- split(df, df$Locus)
        
        p_values_list <- lapply(split_data, function(data) {
            wilcox_test <- wilcox.test(Percent ~ Group, data = data)
            data.frame(Locus = unique(data$Locus), p_value = wilcox_test$p.value)
        })
        
        p_values <- do.call(rbind, p_values_list)
        
        p_values$p_adjusted <- p.adjust(p_values$p_value, method = p_adjust_method)
        p_values$significance <- p_values$p_adjusted < 0.05
        
        if(verbose){
            show(show(p_values))
        }
        
        out <- list(Wilcox = p_values, input = input)
    }
    
    if(method == "ASR"|method == "Logit"){
        
        if(method == "ASR"){
            df[['Transformed_Percent']] <- asin(sqrt(df$Percent / 100))
            
        }else{
            df[['Transformed_Percent']] <- log(df$Percent +0.001 / (100 - df$Percent +0.001 ))
        }
        temp <- aggregate(Transformed_Percent ~ Group + Locus, data = df, function(x) c(mean = mean(x), sd = sd(x)))
        
        summary_data <- data.frame(
        Locus = temp$Locus,
        Group = temp$Group,
        Mean_Transformed_Data = temp$Transformed_Percent[, 'mean'],
        SD_Transformed_Data = temp$Transformed_Percent[, 'sd']
        )
        summary_data$Locus <- factor(summary_data$Locus, levels = locus_levels)
        summary_data <- summary_data[order(summary_data$Locus), ]
        
        shapiro_output <- shapiro.test(df$Transformed_Percent)
        
        if(verbose){
            show(shapiro_output)
        }
        
        
        if(shapiro_output$p.value < 0.05){
            stop("Shapiro test rejected normality. \n")
            
        }else{
            split_data <- split(df, df$Locus)
            p_values_list <- lapply(split_data, function(data) {
                t_test <- t.test(Transformed_Percent ~ Group, data = data)
                data.frame(Locus = unique(data$Locus), p_value = t_test$p.value)
            })
            
            p_values <- do.call(rbind, p_values_list)
            
            p_values$p_adjusted <- p.adjust(p_values$p_value, method = p_adjust_method)
            p_values$significance <- p_values$p_adjusted < 0.05
            
            
            if(verbose){
                show(show(p_values))
            }
            if(method == "ASR"){
                out <- list(ASR = list(transformed_data = df, transformed_summary_data = summary_data,
                T_test = p_values,
                input = input)
                )
                
            }else{
                out <- list(Logit = list(transformed_data = df, transformed_summary_data = summary_data,
                T_test = p_values,
                input = input)
                )
            }
            
            
            
        }
    }
    

    
    if(percentTimer){
        if(length(x@Tocky$TockyLocusStats$percentTimer[['Stats']]) ==0){
            x@Tocky$TockyLocusStats$percentTimer[['Stats']] <- out[[method]]
        }else{
            x@Tocky$TockyLocusStats$percentTimer$Stats[[method]] <- out
        }
        
        
    }else{
        x@Tocky$TockyLocusStats$percentParent[['Stats']] <- out
        
    }
    
    return(x)
}

#' Plot Density of Angles by Group Using Ridge Plots
#'
#' This function takes a TockyPrepData object, which should have been previously
#' processed using the `timer_transform` function, and creates a ridge plot showing
#' the density distribution of angles for each group defined in the dataset.
#'
#' @param x A TockyPrepData object that has been processed with the `timer_transform` function.
#' @param alpha A number between 0 and 1 to be usedby ggridges.
#' @param scale A scaling factor to scale the height of the ridgelines. Used by ggridges.
#' @param group_order Optional. A character vector to define the order of group
#' @param legend Logical. If TRUE, legend is included.
#'
#' @return A ggplot object showing the density distribution of angles by group.
#'
#' @import ggplot2
#' @import ggridges
#' @import RColorBrewer
#' @examples
#' \dontrun{
#' plotAngleDensity(x)
#'}
#'
#' @export
#' @importFrom stats setNames
#' @importFrom rlang sym
#' @importClassesFrom TockyPrep TockyPrepData
plotAngleDensity <- function(x, alpha = 0.3, group_order = NULL, scale = 2, legend = FALSE){
    if (!inherits(x, "TockyPrepData")) {
        stop("Input must be a TockyPrepData object processed with the timer_transform function.", call. = FALSE)
    }
    
    df <- x@Data
    if ("file" %in% names(x@sampledef$sampledef) && "group" %in% names(x@sampledef$sampledef)) {
        df <- merge(df, x@sampledef$sampledef, by = 'file')
    } else {
        stop("The 'sampledef' slot must contain both 'file' and 'group' columns.", call. = FALSE)
    }
    
    df <- df[!is.na(df$Angle),]
    unique_groups <- unique(df$group)
    num_groups <- length(unique_groups)+1
    
    if(is.null(group_order)){
        df$group <- as.factor(df$group)
    }else{
        df$group <- factor(df$group, levels = group_order)
    }
    
    
    if (num_groups > 1) {
        if (num_groups <= 8) {
            colors <- brewer.pal(num_groups, "Set1")
        } else {
            colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_groups)
        }
    } else {
        colors <- "#56B4E9"
    }
    colors <- colors[1:length(unique_groups)]
    colors_vector <- setNames(colors, unique_groups)
    
    p <- ggplot(df, aes(x = !!sym('Angle'), y = !!sym('group'), fill = !!sym('group'))) +
        geom_density_ridges(alpha = alpha, scale = scale) +
        scale_fill_manual(values = colors_vector) +
        labs(title = "Density Plots of Angle by Group",
             x = "Angle",
             y = "Group",
             fill = "Group") +
        theme_ridges(grid = TRUE)
        
    if(legend){
        p <- p +   theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank())
    } else {
        p <- p + theme(legend.position = "none",
                       axis.title.y = element_blank(),
                       axis.text.y = element_text(size = 12),
                       axis.ticks.y = element_blank())
    }


    return(p)
}
