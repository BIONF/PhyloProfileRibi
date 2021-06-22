#' Highlight gene and/or taxon of interest on the phylogenetic profile plot
#' @export
#' @usage highlightProfilePlotCr(data, plotParameter = NULL, taxonHighlight =
#'     "none", rankName = "none", geneHighlight = "none")
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param plotParameter plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2+3) names of 2 variables var1ID and var2ID -
#' default = "var1" & "var2"; (4) color for lowest var1 - default = "#FF8C00";
#' (5) color for highest var1 - default = "#4682B4"; (6) color for lowest var2 -
#' default = "#FFFFFF", (7) color for highest var2 - default = "#F0E68C", (8)
#' color of co-orthologs - default = "#07D000"; (9+10+11) text sizes for x, y
#' axis and legend - default = 9 for each; (12) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (13) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (14) angle of x-axis from 0 to
#' 90 - default = 60; (14) show/hide separate line for reference taxon 1/0 -
#' default = 0; (15) enable/disable coloring gene categories TRUE/FALSE -
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
    rankName = "none", geneHighlight = "none"
){
    if (is.null(data)) stop("Input data cannot be NULL!")
    xmin <- xmax <- ymin <- ymax <- NULL
    p <- PhyloProfile::heatmapPlotting(data, plotParameter)
    # highlight taxon
    if (taxonHighlight != "none") {
        # get selected highlight taxon ID
        nameReducedFile <- paste(
            system.file(package = "PhyloProfileRibi"),
            "PhyloProfile/data/taxonNamesReduced.txt", sep="/")
        if (!file.exists(nameReducedFile)) {
            taxonNamesReduced <- NULL
            delayedAssign("taxName", taxonNamesReduced)
        } else
            taxName <- utils::read.table(nameReducedFile, sep="\t", header=TRUE)
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
        p <- PhyloProfile::heatmapPlotting(data, plotParameter) + geom_rect(
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
        p <- PhyloProfile::heatmapPlotting(data, plotParameter) + geom_rect(
            data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
    }
    return(p)
}
