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
Get1val <- function(Name,one.w,x){
    require(magrittr)
    lapply(split(x,droplevels(seqnames(x))),function(zz){
        message(unique(as.character(seqnames(zz))))
        cov <- one.w[[unique(as.character(seqnames(zz)))]]
        score <- IRanges::Views( cov, start = start(zz), end = end(zz) ) %>% sum()
        tidyr::tibble(wig = Name,value = score,rowname = zz$name)
    }) %>% bind_rows()
}
