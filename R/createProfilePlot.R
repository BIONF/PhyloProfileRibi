#' Create profile heatmap plot
#' @export
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param parm plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2+3) names of 2 variables var1ID and var2ID -
#' default = "var1" & "var2"; (4+5) mid value and color for mid value of var1 -
#' default is 0.5 and #FFFFFF; (6) color for lowest var1 - default = "#FF8C00";
#' (7) color for highest var1 - default = "#4682B4"; (8+9) mid value and color 
#' for mid value of var2 - default is 1 and #FFFFFF;(10) color for lowest var2 -
#' default = "#FFFFFF", (11) color for highest var2 - default = "#F0E68C", (12)
#' color of co-orthologs - default = "#07D000"; (13+14+15) text sizes for x, y
#' axis and legend - default = 9 for each; (16) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (17) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (18) angle of x-axis from 0 to
#' 90 - default = 60; (19) show/hide separate line for reference taxon 1/0 -
#' default = 0; (20) enable/disable coloring gene categories TRUE/FALSE -
#' default = FALSE). NOTE: Leave blank or NULL to use default values.
#' @param rank working taxonomy rank
#' @return A profile heatmap plot as a ggplot object.
#' @import ggplot2
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("superTaxonProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(superTaxonProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE
#' )
#'
#' heatmapPlottingCr(plotDf, plotParameter, "species")

heatmapPlottingCr <- function(data = NULL, parm = NULL, rank = "species"){
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (is.null(parm))
        parm <- list(
            "xAxis" = "taxa", "var1ID" = "var1", "var2ID"  = "var2",
            "midVar1" = 0.5, "midColorVar1" =  "#FFFFFF",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "midVar2" = 1, "midColorVar2" =  "#FFFFFF",
            "lowColorVar2" = "#CB4C4E", "highColorVar2" = "#3E436F",
            "paraColor" = "#07D000", "xSize" = 8, "ySize" = 8, "legendSize" = 8,
            "mainLegend" = "top", "dotZoom" = 0, "xAngle" = 60, "guideline" = 0,
            "colorByGroup" = FALSE)
    geneID <- supertaxon <- category <-var1<-var2 <- presSpec <- paralog <- NULL
    xmin <- xmax <- ymin <- ymax <- NULL
    # get font and dot size
    dotZoom <- parm$dotZoom
    xSize <- parm$xSize
    ySize <- parm$ySize
    if (rank == "species") {
        if (dotZoom > -0.2) dotZoom = -0.2
        if (xSize > 12) xSize = 12
        if (ySize > 12) ySize = 12
    }
    
    # create heatmap plot with geom_point & scale_color_gradient for present
    # ortho & var1, geom_tile & scale_fill_gradient for var2
    if (parm$xAxis == "genes") p <- ggplot(data,aes(x = geneID, y = supertaxon))
    else p <- ggplot(data, aes(y = geneID, x = supertaxon))
    if (parm$colorByGroup == TRUE) {
        p <- p + geom_tile(aes(fill = factor(category)), alpha = 0.3)
    } else {
        if (length(unique(stats::na.omit(data$var2))) != 1)
            p <- p + scale_fill_gradient2(
                low = parm$lowColorVar2, high = parm$highColorVar2,
                mid = parm$midColorVar2, midpoint = parm$midVar2,
                na.value = "gray95", limits = c(0, 1)) +
                geom_tile(aes(fill = var2))
    }
    if (length(unique(stats::na.omit(data$presSpec))) < 3) {
        if (length(unique(stats::na.omit(data$var1))) == 1) {
            p <- p + geom_point(
                aes(colour = var1), na.rm = TRUE,
                size = data$presSpec*5*(1+dotZoom), show.legend = FALSE
            )
        } else
            p <- p +
                geom_point(
                    aes(colour = var1), na.rm = TRUE,
                    size = data$presSpec * 5 * (1 + dotZoom)) +
                scale_color_gradient2(
                    low = parm$lowColorVar1, high = parm$highColorVar1,
                    mid = parm$midColorVar1, midpoint = parm$midVar1,
                    limits = c(0,1)
                )
    } else {
        if (length(unique(stats::na.omit(data$var1))) == 1) {
            p <- p + geom_point(aes(size=presSpec),color="#336a98",na.rm = TRUE)
        } else
            p <- p + geom_point(aes(colour=var1, size = presSpec),na.rm = TRUE)+
                scale_color_gradient2(
                    low = parm$lowColorVar1, high = parm$highColorVar1,
                    mid = parm$midColorVar1, midpoint = parm$midVar1,
                    limits = c(0,1)
                )
    }
    # plot inparalogs (if available)
    if (length(unique(stats::na.omit(data$paralog))) > 0) {
        p <- p +
            geom_point(
                data = data, aes(size = paralog), color = parm$paraColor,
                na.rm = TRUE, show.legend = TRUE) +
            guides(size = guide_legend(title = "# of co-orthologs")) +
            scale_size_continuous(range = c(
                min(stats::na.omit(data$paralogSize)) * (1 + dotZoom),
                max(stats::na.omit(data$paralogSize)) * (1 + dotZoom)))
    } else {
        # remain the scale of point while filtering
        presentVl <- data$presSpec[!is.na(data$presSpec)]
        p <- p + scale_size_continuous(range = c(
            (floor(min(presentVl) * 10) / 10 * 5) * (1 + dotZoom),
            (floor(max(presentVl) * 10) / 10 * 5) * (1 + dotZoom)))
    }
    # color gene categories
    if (parm$colorByGroup == FALSE) {
        p <- p + guides(fill = guide_colourbar(title = parm$var2ID),
                        color = guide_colourbar(title = parm$var1ID))
    } else
        p <- p + guides(fill = guide_legend("Category"),
                        color = guide_colourbar(title = parm$var1ID))
    # guideline for separating ref species
    if (parm$guideline == 1) {
        if (parm$xAxis == "genes") {
            p <- p + labs(y = "Taxon") +
                geom_hline(yintercept = 0.5, colour = "dodgerblue4") +
                geom_hline(yintercept = 1.5, colour = "dodgerblue4")
        } else
            p <- p + labs(x = "Taxon") +
                geom_vline(xintercept = 0.5, colour = "dodgerblue4") +
                geom_vline(xintercept = 1.5, colour = "dodgerblue4")
    }
    p <- p + theme_minimal(base_size = 9)
    vjustValue <- 1
    if (parm$xAngle == 90) vjustValue <- 0.5
    p <- p + theme(
        axis.text.x = element_text(
            angle = parm$xAngle, hjust = 1, size = xSize, 
            vjust = vjustValue
        ),
        axis.text.y = element_text(size = ySize),
        axis.title.x = element_text(size = xSize),
        axis.title.y = element_text(size = ySize),
        legend.title = element_text(size = parm$legendSize),
        legend.text = element_text(size = parm$legendSize),
        legend.position = parm$mainLegend)
    return(p)
}

#' Highlight gene and/or taxon of interest on the phylogenetic profile plot
#' @export
#' @usage highlightProfilePlotCr(data, plotParameter = NULL, taxonHighlight =
#'     "none", rankName = "none", geneHighlight = "none")
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param plotParameter plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2+3) names of 2 variables var1ID and var2ID -
#' default = "var1" & "var2"; (4+5) mid value and color for mid value of var1 -
#' default is 0.5 and #FFFFFF; (6) color for lowest var1 - default = "#FF8C00";
#' (7) color for highest var1 - default = "#4682B4"; (8+9) mid value and color
#' for mid value of var2 - default is 1 and #FFFFFF;(10) color for lowest var2 -
#' default = "#FFFFFF", (11) color for highest var2 - default = "#F0E68C", (12)
#' color of co-orthologs - default = "#07D000"; (13+14+15) text sizes for x, y
#' axis and legend - default = 9 for each; (16) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (17) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (18) angle of x-axis from 0 to
#' 90 - default = 60; (19) show/hide separate line for reference taxon 1/0 -
#' default = 0; (20) enable/disable coloring gene categories TRUE/FALSE -
#' default = FALSE). NOTE: Leave blank or NULL to use default values.
#' @param taxonHighlight taxon of interst. Default = "none".
#' @param rankName working taxonomy rank (needed only for highlight taxon).
#' @param geneHighlight gene of interest. Default = "none".
#' @import ggplot2
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}

highlightProfilePlotCr <- function(
    data = NULL, plotParameter = NULL, taxonHighlight = "none",
    rankName = "species", geneHighlight = "none"
){
    if (is.null(data)) stop("Input data cannot be NULL!")
    xmin <- xmax <- ymin <- ymax <- NULL
    p <- heatmapPlottingCr(data, plotParameter, rankName)
    # highlight taxon
    if (taxonHighlight != "none") {
        # get selected highlight taxon ID
        nameReducedFile <- paste(
            system.file(package = "PhyloRBF"),
            "PhyloProfile/data/taxonNamesReduced.txt", sep="/")
        if (!file.exists(nameReducedFile)) {
            taxonNamesReduced <- NULL
            delayedAssign("taxName", taxonNamesReduced)
        } else {
            taxName <- utils::read.table(
                nameReducedFile, sep="\t", header=TRUE, comment.char = ""
            )
        }
        taxonHighlightID <- taxName$ncbiID[
            taxName$fullName == taxonHighlight & taxName$rank == rankName]
        if (length(taxonHighlightID) == 0L)
            taxonHighlightID <- taxName$ncbiID[taxName$fullName==taxonHighlight]
        # get taxonID together with it sorted index
        selTaxon <- toString(data[data$supertaxonID == taxonHighlightID, 2][1])
        selIndex <- grep(selTaxon, levels(as.factor(data$supertaxon)))
        if (plotParameter$xAxis == "taxa") {
            rect <- data.frame(
                xmin=selIndex-0.5, xmax = selIndex+0.5, ymin = -Inf, ymax = Inf)
        } else
            rect <- data.frame(
                ymin=selIndex-0.5, ymax = selIndex+0.5, xmin = -Inf, xmax = Inf)
        p <- heatmapPlottingCr(data, plotParameter, rankName) + geom_rect(
            data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
    }
    # highlight gene
    if (geneHighlight != "none") {
        selIndex <- match(geneHighlight, levels(as.factor(data$geneID)))
        if (plotParameter$xAxis == "taxa") {
            rect <- data.frame(
                ymin=selIndex-0.5, ymax = selIndex+0.5, xmin = -Inf, xmax = Inf)
        } else
            rect <- data.frame(
                xmin=selIndex-0.5, xmax = selIndex+0.5, ymin = -Inf, ymax = Inf)
        p <- heatmapPlottingCr(data, plotParameter, rankName) + geom_rect(
            data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
    }
    return(p)
}
