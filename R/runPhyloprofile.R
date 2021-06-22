#' Run PhyloProfile app
#' @export
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import DT
#' @importFrom colourpicker colourInput
#' @import energy
#' @import shinyBS
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(shinyjs, except = colourInput)

runPhyloProfileRibi <- function(){
    appDir <- system.file("PhyloProfile", package = "PhyloProfileRibi")
    if (appDir == "") {
        stop(
            "Could not find apps directory. Try re-installing `PhyloProfileRibi`.",
            call = FALSE
        )
    }

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
