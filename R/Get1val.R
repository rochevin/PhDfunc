#' Get1val : Get the coverage of a RleList object for a list of genomic positions
#'
#' @param Name The name of the bigwig file
#' @param one.w A bigwig file imported in RleList format (as = "RleList")
#' @param x A GRanges object
#'
#' @return A tibble with the coverage of each element of x
#' @export
#'
#' @examples
#' to.plot.4C <- sres <- lapply(c(100000,500000,1000000),function(bin){
#' values <-  bless80 %>% anchor_center() %>% mutate(width=bin)
#' lapply(wigs.4C,function(wig){
#'     one.w <- import.bw(str_c(path4C,wig,sep="/"),as="RleList")
#'     File <- replaceName[[str_remove(wig,"_chr.*.bw")]]
#'     vp <- str_extract(wig,"chr[0-9A-Z]+_[0-9]+-[0-9]+") %>% str_replace("_",":")
#'     x <- values %>% filter_by_non_overlaps(GRanges(vp))
#'     Get1val(File,one.w,x) %>% mutate(binsize = bin) %>% mutate(viewpoint = vp)
#' }) %>% bind_rows()
#' }) %>% bind_rows()
Get1val <- function(Name,one.w,x){
    lapply(split(x,droplevels(seqnames(x))),function(zz){
        message(unique(as.character(seqnames(zz))))
        cov <- one.w[[unique(as.character(seqnames(zz)))]]
        score <- rtracklayer::Views( cov, start = start(zz), end = end(zz) ) %>% sum()
        tidyr::tibble(wig = Name,value = score,rowname = zz$name)
    }) %>% bind_rows()
}
