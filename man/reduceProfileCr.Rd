% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsePhyloProfile.R
\name{reduceProfileCr}
\alias{reduceProfileCr}
\title{Reduce the filtered profile data into supertaxon level}
\usage{
reduceProfileCr(filteredProfile)
}
\arguments{
\item{filteredProfile}{dataframe contains the filtered profiles (see
?parseInfoProfile, ?filterProfileData and ?filteredProfile)}
}
\value{
A reduced dataframe contains only profile data for the selected
supertaxon rank. This dataframe contains only supertaxa and their value
(mVar1 & mVar2) for each gene.
}
\description{
Reduce data of the processed phylogenetic profiles from input
taxonomy rank into supertaxon level (e.g. from species to phylum)
}
\author{
Vinh Tran {tran@bio.uni-frankfurt.de}
}
