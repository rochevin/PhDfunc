#' Title
#'
#' @param my.wigs
#' @param one.w
#' @param x
#'
#' @return
#' @export
#'
#' @examples
Get1val <- function(my.wigs,one.w,x){
    lapply(split(x,droplevels(seqnames(x))),function(zz){
        message(unique(as.character(seqnames(zz))))
        cov <- one.w[[unique(as.character(seqnames(zz)))]]
        score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
        tibble(wig = my.wigs,value = score,rowname = zz$name)
    }) %>% bind_rows()
}
