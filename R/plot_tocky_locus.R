# Copyright 2024 Masahiro Ono
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Generate basic QC plots for Tocky data (Timer-Blue vs Timer-Red 2d plots)
#'
#' This function creates quick control plots for the TockyPrepData object
#' analyzing fluorescence changes over time in cellular activities.
#'
#' @param x A TockyPrepData object produced by the function `prep_tocky`.
#' @param file The name of the output file.
#' @param n The number of plots per row and column in the output grid.
#' @param max_cell_number The maximum number of cells to be displayed per panel.
#' @param viridis (Optional). If TRUE, a colour-blind friendly colour set is used.
#' @param interactive (Optional). If TRUE, an interactive session is used to trim plot area.
#'
#' @return An unchanged TockyPrepData object, primarily for consistency in pipeline usage.
#' @export
#' @examples
#' \dontrun{
#' plot_tocky_locus(data)
#'}
#' @importFrom methods show
#' @importFrom utils askYesNo
#' @importFrom graphics locator rect par
#' @importFrom grDevices dev.off
#' @importFrom stats quantile
#' @importClassesFrom TockyPrep TockyPrepData
plot_tocky_locus <- function(x, file='PlotTockyLocus', n = 3, max_cell_number = 20000, viridis = FALSE, interactive = FALSE){

    if(!inherits(x, "TockyPrepData")){
        stop("Use a TockyPrepData for x. \n")
    }
    
    sampledef <- x@sampledef$sampledef
    z <- x@Data[, c('Red_log','Blue_log',"Angle","Intensity", 'file')]
    z <- merge(z, sampledef, by = 'file')
    x.var <- 'Red_log'
    y.var <- 'Blue_log'
    z$group <-as.factor(z$group)
    cat('.')
    X <- split(z, z$group)
    
    df <- matrix(rep(1, length(levels(z$file))*4), nrow=4)
    df <- as.data.frame(df); rownames(df)<-paste("Q",1:4, sep='');colnames(df) <- levels(z$file)
    cat('.')
    
    tmpdata <- z[,c('Red_log','Blue_log')]
    
    if(nrow(tmpdata) > max_cell_number ){
        tmpdata  <-  tmpdata[sample(1:nrow(tmpdata), max_cell_number),]
    }
    
    var.x.name <- c('Red_log','Blue_log')[1]; var.y.name <- c('Red_log','Blue_log')[2]
    xl <- c(quantile(tmpdata$Red_log, 0.05), max(tmpdata$Red_log))
    yl <- c(quantile(tmpdata$Blue_log, 0.05), max(tmpdata$Blue_log))
    

    
    if(interactive){
        plot(tmpdata, pch='.', col = rgb(0,0,0,alpha=0.2), xlim = xl, ylim = yl)
        x.scalegate <- c()
        y.scalegate <- c()
        
        scaleanswer <- askYesNo("Need to change the limits of X and Y Axes?")
        if(scaleanswer){
            repeat{
                scalegate <- locator(type = 'p', col=2)
                x.scalegate <- c(scalegate[[1]])
                xlim <- sort(x.scalegate)
                y.scalegate <- c(scalegate[[2]])
                ylim <- sort(y.scalegate)
                if(length(x.scalegate)==2){
                    cat("Got it! \n")
                    rect(xleft = xlim[1], ybottom = ylim[1], xright = xlim[2], ytop = ylim[2], border = 2)
                    answer3 <- askYesNo("Happy with your scaling? \n")
                    if(answer3){
                         break
                    }else{
                        dev.off()
                        plot(tmpdata, pch='.', col=rgb(0,0,0,alpha=0.2), xlim = xl, ylim = yl)
                    }
                }else{
                    paste("choose exactly 2 points!")
                }
            }
        }
        
        if(!scaleanswer){
            xlim <- c(quantile(tmpdata[,1], 0.05), max(tmpdata[,1]))
            ylim <- c(quantile(tmpdata[,2], 0.05), max(tmpdata[,2]))

            
        }
    }else{
        
        xlim <- c(quantile(tmpdata[,1], 0.05), max(tmpdata[,1]))
        ylim <- c(quantile(tmpdata[,2], 0.05), max(tmpdata[,2]))
        }


    group <- x@sampledef$sampledef[,x@sampledef$Group]
    

        par(mfrow=c(n,n))
        par(mar = c(4, 4, 2.5, 2.5))
        for(i in 1:length(X)){
            
            tpdata <- X[[i]][,c("Red_log", "Blue_log")]
            colnames(tpdata) <- c('x','y')
          #
           if(nrow(tpdata) > max_cell_number ){
               logic <- sample(1:nrow(tpdata), max_cell_number)
                tpdata  <-  tpdata[logic,]
                angles <- X[[i]]$Angle[logic]
            }
#
              density  <-  Locus_to_colour(angles, viridis = viridis)
              title <- names(X)[i]

                plot(tpdata, xlab="Timer Red", ylab="Timer Blue", pch=19, cex=0.2, col=density,xlim = xlim, ylim = ylim, cex.main=1.4, cex.lab = 1.2, main = title)



        }
        TockyLocusLegend(viridis = viridis)

}



#' Generate Tocky Locus Legend For a Plot
#' @param legend A logical arguement.
#' @param cex A numeric value for the text size.
#' @param viridis (Optional). If TRUE, a colour-blind friendly colour set is used.
#' @examples
#' \dontrun{
#' TockyLocusLegend()
#'}
#' @export

TockyLocusLegend <- function(legend = TRUE, cex = 1, viridis = FALSE){
    
    if(legend){
        if(!viridis){
            df <- data.frame(
                colour = c(rgb(0,0,1,alpha=0.5),
                rgb(0.25,0,0.75,alpha=0.3),
                rgb(0.5,0,0.5,alpha=0.5),
                rgb(0.75,0,0.25,alpha=0.3),
                rgb(1,0,0,alpha=0.3),
                rgb(0,0,0,alpha=0.1)
                ),
                population <- c("New", "NPt", "Persistent", "PAt", "Arrested", "Negative"))
        }else {
            palette <- viridis::viridis(5, alpha = c(0.5, 0.3, 0.5, 0.3, 0.3), begin = 0, end = 1, option = "D")
            palette <- c(palette, 'grey')
            df <- data.frame(
                colour = palette,
                population <- c("New", "NPt", "Persistent", "PAt", "Arrested", "Negative"))
             
        }

            
        
        par(mar = c(2, 2, 2,2))
        plot(x = c(-1,1), y = c(-1,1),xaxt="n",yaxt="n",bty="n", type = 'n', ann=F)
        
        text(-0.5, 1, labels = 'Tocky Locus:', cex = cex)
        
        legend(-1,0.9, col = as.vector(df$colour),  pch=19, legend=df$population, ncol=1, cex = cex)
    }

}




    
#' Convert Timer Angle Data into color code
#'
#' This function assigns colors to different ranges of angle values, with an option
#' to use colorblind-friendly colors from the viridis palette.
#'
#' @param x Angle numeric vector.
#' @param viridis Logical, whether to use the viridis color palette.
#' @return a character vector for color code.
#' @examples
#' \dontrun{
#' col <- Locus_to_colour(x = c(0, 25, 45, 65, 90), viridis = TRUE)
#'}
#'
#' @importFrom grDevices rgb
#' @importFrom viridis viridis
#' @export

Locus_to_colour <- function(x, viridis = FALSE){
    nalogic <- is.na(x)
    colour <- x
    x[is.na(x)] <- 900
    new_lg <- x == 0
    persistent_lg <- x >= 30 & x < 60
    np_lg <- x > 0 & x < 30
    pa_lg <- x > 60 & x < 90
    arrested_lg <- x == 90
    
    if (!viridis) {
        colour <- ifelse(new_lg, rgb(0, 0, 1, alpha = 0.5),
                         ifelse(arrested_lg, rgb(1, 0, 0, alpha = 0.3),
                                ifelse(persistent_lg, rgb(0.5, 0, 0.5, alpha = 0.5),
                                       ifelse(np_lg, rgb(0.25, 0, 0.75, alpha = 0.3),
                                              ifelse(pa_lg, rgb(0.75, 0, 0.25, alpha = 0.3),
                                                     rgb(0, 0, 0, alpha = 0.2)))))
                                                     )
        
        colour[nalogic] <- rgb(0, 0, 0, alpha = 0.1)
    } else {
        palette <- viridis::viridis(5, alpha = c(0.5, 0.3, 0.5, 0.3, 0.3), begin = 0, end = 1, option = "D")
        palette <- c(palette, 'grey')
        names(palette) <- c("new", "arrested", "persistent", "np", "pa")
        
        colour <- ifelse(new_lg, palette[1],
                         ifelse(arrested_lg, palette[5],
                                ifelse(persistent_lg, palette[3],
                                       ifelse(np_lg, palette[2],
                                              ifelse(pa_lg, palette[4], palette[6])))))
                                                     
        
    }
    
    return(colour)
}



#' Plot Coloured Tocky Locus Legend
#' @param mar_par parameters for the function mar. The default is c(4, 4, 10, 4)
#' @return A plot with colored rectangles and labels.
#' @examples
#' \dontrun{
#'   plotTockyLocusLegend()
#' }
#' @export
#' @importFrom graphics text axis


plotTockyLocusLegend <- function(mar_par = c(4, 4, 10, 4)) {
  categories <- c("New", "Persistent", "Arrested")
  colors <- c("#ADD8E6", "#9370DB", "#FF69B4")
  labels <- c("NPt", "Persistent", "PAt")
  positions <- c(0, 30, 60, 90)
  rectangle_width <- 2
  par(mar = mar_par)
  # Plotting code as provided
  plot(0, type="n", xlim=c(-5, 98), ylim=c(0, 1.4), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", main="")
  rect(positions[1]-rectangle_width, 0, positions[1]+rectangle_width, 1, col="blue", border="black")
  rect(positions[4]-rectangle_width, 0, positions[4]+rectangle_width, 1, col="red", border="black")
  rect(positions[1]+rectangle_width, 0, positions[2], 1, col=colors[1], border="black")
  rect(positions[2], 0, positions[3], 1, col=colors[2], border="black")
  rect(positions[3], 0, positions[4]-rectangle_width, 1, col=colors[3], border="black")
  text((positions[1]+rectangle_width + positions[2]) / 2, 0.5, labels[1], cex=1.5)
  text((positions[2] + positions[3]) / 2, 0.5, labels[2], cex=1.5)
  text((positions[3] + positions[4]-rectangle_width) / 2, 0.5, labels[3], cex=1.5)
  axis(1, at=positions, labels=c("0", "30", "60", "90"), cex.axis=1.5)
  text(positions[1], 1.1, categories[1], cex=1.5, srt=0, adj=c(0.5, 0))
  text(positions[4], 1.1, categories[3], cex=1.5, srt=0, adj=c(0.5, 0))
}


