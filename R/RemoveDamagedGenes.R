#' Remove Damaged Genes from table
#'
#' @param my.vals A data.frame or tibble wich contains Gene id as rowname
#' @param DSBpath Path of a bedfile wich contains genomic position to exclude
#' @param GenePath Path of Gene annotation (ensembl format)
#'
#' @return my.vals without Genes that overlap with DSBpath
#' @export
#'
#' @examples
RemoveDamagedGenes <- function(my.vals,DSBpath = "/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed",
                               GenePath = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed"){
    DSB174 <- rtracklayer::import.bed(DSBpath)
    ens.genes <- read.table(GenePath,sep="\t",h=T) %>% GRanges()
    ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
    seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))


    require(plyranges)
    DSB.extend <- DSB174 %>% anchor_center() %>% mutate(width = 1000000)
    gene.in.DiVA <- ens.genes %>% filter_by_overlaps(DSB.extend)

    my.vals %>% filter(!rowname %in% gene.in.DiVA$gene_id)

}
