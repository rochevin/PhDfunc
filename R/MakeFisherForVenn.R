#' Fisher Test for a list of name
#'
#' @param venn A list of names
#' @param universe Universe of names
#'
#' @return A vector of pval by fisher.test
#' @export
#'
#' @examples
MakeFisherForVenn <- function(venn,universe){
    possibility <- combn(names(venn),2)
    pval <- NULL
    for(i in 1:ncol(possibility)){
        A.n <- possibility[1,i]
        B.n <- possibility[2,i]
        A <- venn[[A.n]]
        B <- venn[[B.n]]
        common <- (A %in% B) %>% .[. == TRUE] %>% length()
        AvsB <- matrix(c(common,length(B) - common,length(A) - common,length(universe) - length(A) - length(B) + common),ncol=2,byrow = T)
        res <- fisher.test(AvsB)$p.value
        pval <- c(pval,res)
    }
    names(pval) <- apply(possibility,2,function(x){paste(x[1],x[2],sep="_")})
    pval
}
