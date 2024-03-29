% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetDE_DIvA.R
\name{GetDE_DIvA}
\alias{GetDE_DIvA}
\title{Get RNA-Seq table DE results for DIvA.}
\usage{

  GetDE_DIvA(path = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Legube/PROCESSED_strand-spe-reverse_10102018/mapping/counts",
  pattern = "_counts.tsv", exclude = "STX", Conditions = c("C1",
  "C1_OHT"), Labels = paste(Conditions, rep(c("1", "2"), each = 2), sep =
  "_"), p.cutoff = 0.1, LogFC.cutoff = 0.5, to.test = "CvsC_OHT")
}
\arguments{
\item{path}{The path directory of table counts matrix (from HTSeq Count)}

\item{pattern}{The pattern of the extension (default "_counts.tsv")}

\item{exclude}{An regex expression to exclude some conditions, default "STX"}

\item{Conditions}{A vector of conditions names : c("C1","C1_OHT")}

\item{Labels}{The names of columns, default : paste(Conditions,rep(c("1","2"),each=2),sep="_")}

\item{p.cutoff}{The adjusted p.value cutoff used to determine wich gene is differentially expressed}

\item{LogFC.cutoff}{The log2(foldchange) cutoff used to determine wich gene is differentially expressed}

\item{to.test}{Which condition to test ? Default "CvsC_OHT"}
}
\value{
A data.frame with all expressed genes in this Condition
}
\description{
Get RNA-Seq table DE results for DIvA.
}
