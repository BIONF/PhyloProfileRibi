#' Get list of pre-installed NCBI taxon names
#' @description Get all NCBI taxon names from
#' "PhyloProfile/data/taxonNamesReduced.txt"
#' @export
#' @return List of taxon IDs, their full names, taxonomy ranks and parent IDs
#' obtained from "PhyloProfile/data/taxonNamesReduced.txt"
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getNameListCr <- function() {
    nameReducedFile <- paste(
        system.file(package="PhyloRBF"),
        "PhyloProfile/data/taxonNamesReduced.txt",
        sep="/"
    )

    if (!file.exists(nameReducedFile)) {
        utils::data(taxonNamesReduced)
    } else {
        taxonNamesReduced <- utils::read.table(
            nameReducedFile, sep = "\t", header = TRUE, fill = TRUE, 
            comment.char = ""
        )
    }

    taxonNamesReduced$fullName <- as.character(taxonNamesReduced$fullName)
    taxonNamesReduced$rank <- as.character(taxonNamesReduced$rank)
    taxonNamesReduced <- taxonNamesReduced[!duplicated(taxonNamesReduced), ]

    return(taxonNamesReduced)
}

#' Get taxonomy matrix
#' @description Get the (full or subset) taxonomy matrix from
#' "data/taxonomyMatrix.txt" based on an input taxon list
#' @export
#' @param subsetTaxaCheck TRUE/FALSE subset taxonomy matrix based on input taxon
#' IDs. Default = FALSE.
#' @param taxonIDs list of input taxon IDs (e.g. ncbi1234). Default = NULL.
#' @return Data frame contains the (subset of) taxonomy matrix for list of
#' input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getTaxonomyMatrixCr <- function(subsetTaxaCheck = FALSE, taxonIDs = NULL){
    taxonomyMatrixFile <- paste(
        system.file(package="PhyloRBF"),
        "PhyloProfile/data/taxonomyMatrix.txt",
        sep="/"
    )

    if (!file.exists(taxonomyMatrixFile)) {
        utils::data(taxonomyMatrix)
    } else {
        taxonomyMatrix <- utils::read.table(
            taxonomyMatrixFile, sep = "\t", header = TRUE,
            stringsAsFactors = TRUE
        )
    }

    if (subsetTaxaCheck) {
        if (missing(taxonIDs)) return(taxonomyMatrix)
        taxonomyMatrix <- taxonomyMatrix[
            taxonomyMatrix$abbrName  %in% taxonIDs, ]
    }
    return(taxonomyMatrix)
}

#' Get NCBI taxon names for a selected list of taxa
#' @description Get NCBI taxon names from
#' "PhyloProfile/data/taxonNamesReduced.txt" for a list of input taxa
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param taxonIDs list of taxon IDs (e.g. ncbi1234). Default = NULL
#' @return Data frame contains a list of full names, taxonomy ranks and parent
#' IDs for the input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export

getInputTaxaNameCr <- function(rankName, taxonIDs = NULL){
    # check input parameters
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- PhyloProfile::getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    # load list of unsorted taxa
    Dt <- getTaxonomyMatrixCr(TRUE, taxonIDs)
    # load list of taxon name
    nameList <- getNameListCr()
    # return
    choice <- data.frame(
        "ncbiID" = unlist(Dt[rankName]), stringsAsFactors = FALSE
    )
    choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
    return(choice)
}

#' Get a subset of input taxa based on a selected taxonomy rank
#' @description Get a subset of taxon ncbi IDs and names from an input list of
#' taxa based on a selected supertaxon (identified by its taxonomy rank and
#' supertaxon name or supertaxon ID).
#' @usage getSelectedTaxonNamesCr(inputTaxonIDs, rank, higherRank, higherID,
#'     higherName)
#' @param inputTaxonIDs list of input taxon IDs (e.g. c("10116", "122586"))
#' @param rank taxonomy rank of input taxa (e.g. "species")
#' @param higherRank selected taxonomy rank (e.g. "phylum")
#' @param higherID supertaxon ID (e.g. 7711). NOTE: either supertaxon ID or
#' name is required, not neccessary to give both.
#' @param higherName supertaxon name (e.g. "Chordata"). NOTE: either
#' supertaxon ID or name is required, not neccessary to give both.
#' @export
#' @return A data frame contains ncbi IDs and names of taxa from the input taxon
#' list that belong to the selected supertaxon.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getSelectedTaxonNamesCr <- function(
    inputTaxonIDs = NULL, rank = NULL,
    higherRank = NULL, higherID = NULL, higherName = NULL
) {
    rankName <- NULL
    if (is.null(inputTaxonIDs) | is.null(rank))
        stop("Input taxa and taxonomy rank cannot be NULL!")
    taxDf <- getTaxonomyMatrixCr(TRUE, paste0("ncbi", inputTaxonIDs))
    if (is.null(higherID) & is.null(higherName))
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE))
    if (is.null(higherRank)) {
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE))
    } else {
        if (!is.null(higherName) & is.null(higherID)) {
            taxaList <- getNameListCr()
            superID <- taxaList$ncbiID[
                taxaList$fullName == higherName
                & taxaList$rank %in% c(higherRank, "norank")]
            customizedtaxaID <- levels(
                as.factor(taxDf[rank][taxDf[higherRank] == superID, ]))
            return(
                data.frame(
                    ncbiID = taxaList$ncbiID[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    name = taxaList$fullName[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    stringsAsFactors = FALSE))
        } else if (!is.null(higherID)) {
            return(
                data.frame(
                    ncbiID = taxDf$ncbiID[taxDf[,higherRank] == higherID],
                    name = taxDf$fullName[taxDf[,higherRank] == higherID],
                    stringsAsFactors = FALSE))
        }
    }
}

#' Sort list of (super)taxa based on a selected reference (super)taxon
#' @usage sortInputTaxaCr(taxonIDs = NULL, rankName, refTaxon = NULL,
#'     taxaTree = NULL)
#' @param taxonIDs list of taxon IDs (e.g.: ncbi1234, ncbi9999, ...). Default =
#' NULL.
#' @param rankName working taxonomy rank (e.g. "species", "phylum",...)
#' @param refTaxon selected reference taxon. Default = NULL.
#' @param taxaTree taxonomy tree for the input taxa (optional). Default = NULL.
#' @return A taxonomy matrix for the input taxa ordered by the selected
#' reference taxon. This matrix is sorted either based on the NCBI taxonomy
#' info, or based on an user-defined taxonomy tree (if provided).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export

sortInputTaxaCr <- function(
    taxonIDs = NULL, rankName, refTaxon = NULL, taxaTree = NULL
){
    ncbiID <- fullName <- abbrName <- NULL
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- PhyloProfile::getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    if (is.null(refTaxon))  refTaxon <- taxonNames$fullName[1]
    # get list of taxon names
    fullnameList <- getNameListCr()
    taxonNames <- getInputTaxaNameCr(rankName, taxonIDs)
    # get selected supertaxon ID(s)
    rankNameTMP <- taxonNames$rank[taxonNames$fullName == refTaxon]
    if (rankName == "strain") {
        superID <- fullnameList$ncbiID[fullnameList$fullName == refTaxon]
    } else
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon
            & fullnameList$rank == rankNameTMP[1]]
    # get full taxonomy data & representative taxon
    Dt <- getTaxonomyMatrixCr()
    repTaxon <- Dt[Dt[, rankName] == superID, ][1, ]
    # THEN, SORT TAXON LIST BASED ON TAXONOMY TREE
    if (is.null(taxaTree)) {
        distDf <- subset(Dt, select = -c(ncbiID, fullName))
        row.names(distDf) <- distDf$abbrName
        distDf <- distDf[, -1]
        taxaTree <- PhyloProfile::createRootedTree(
            distDf, as.character(repTaxon$abbrName)
        )
    } else
        taxaTree <- ape::root(
            taxaTree,outgroup=as.character(repTaxon$abbrName),resolve.root=TRUE)
    taxonList <- PhyloProfile::sortTaxaFromTree(taxaTree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]
    # subset to get list of input taxa only
    sortedDt <- subset(sortedDt, abbrName %in% taxonIDs)
    # get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(
        sortedDt,
        select = c("abbrName", "ncbiID", "fullName", as.character(rankName)))
    colnames(sortedOut) <- c("abbrName", "species", "fullName", "ncbiID")
    # add name of supertaxa into sortedOut list
    sortedOut <- merge(
        sortedOut, fullnameList, by = "ncbiID", all.x = TRUE, sort = FALSE)
    sortedOut$species <- paste0("ncbi", sortedOut$species)
    ## create new column for sorted supertaxon
    indexSpec <- unlist(lapply(
        seq_len(nlevels(as.factor(sortedOut$fullName.y))),function (x) 1000+x))
    indexSpecDf <- data.frame(
        fullName.y = unique(as.character(sortedOut$fullName.y)),
        sortedSupertaxon = paste0(
            indexSpec, "_", unique(as.character(sortedOut$fullName.y))
        ), stringsAsFactors = FALSE)
    sortedOut <- merge(indexSpecDf, sortedOut, by = "fullName.y")
    # final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[, c(
        "abbrName", "taxonID", "fullName.x", "species", "ncbiID",
        "sortedSupertaxon", "rank", "category")]
    colnames(sortedOut) <- c(
        "abbrName", "taxonID", "fullName", "ncbiID", "supertaxonID",
        "supertaxon", "rank", "category")
    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)
    return(sortedOut)
}
