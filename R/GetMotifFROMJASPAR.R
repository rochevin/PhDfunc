#' GetMotifFROMJASPAR : Search motifs within DNA promoter sequences
#'
#' @param bed_path A path which contains one or more bedfiles
#' @param pattern extension, default ".bed"
#' @param species The species ID from NCBI, default 9606
#' @param collection collection=c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH")
#' @param myJASPAR A JASPAR connection databases, default JASPAR2018::JASPAR2018
#' @param genome A BSgenome object
#' @param w upstream/downstream for promoter, default 3000
#'
#' @return A data.frame of motif
#' @export
#'
#' @examples
GetMotifFROMJASPAR <- function(bed_path = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/results/COMMON_GDR_RESPONSE/BED",
                               pattern =".bed",
                               species = 9606,
                               collection = "CORE",
                               myJASPAR = JASPAR2018::JASPAR2018,
                               genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                               w = 3000){
    ##GDR/DDR BED
    require(magrittr)
    filesbed <- list.files(bed_path,full.names = T,pattern=pattern)

    #MOTIFS ANALYSIS
    #Search motif in seq from JASPAR
    opts <- list()
    opts[["species"]] <- species
    opts[["collection"]] <- collection

    PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)

    beds <- lapply(filesbed,function(x){
        mybed <- rtracklayer::import.bed(x) %>% GenomicRanges::promoters(upstream = w,downstream = w)
        Bed <- basename(x) %>% str_remove(".bed")
        Groups = Bed %>% str_split_fixed("_",n=Inf)
        colnames(Groups) <- stringr::str_c("Group",1:length(Groups))

        seq = BSgenome::getSeq(genome, mybed)

        motif_ix <- motifmatchr::matchMotifs(PFMatrixList, seq,out = c("scores"))
        motif_id <- colnames(motif_ix)
        res <- motifmatchr::motifMatches(motif_ix)
        counts <- motifmatchr::motifCounts(motif_ix)


        tibble(occ=apply(counts,2,sum),id = motif_id) %>% cbind(Groups)

    }) %>% bind_rows()
}
