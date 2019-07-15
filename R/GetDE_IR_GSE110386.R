GetDE_IR_GSE110386 <- function(path = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_GSE110386/counts",
                               pattern ="_counts.tsv",
                               exclude = NULL,
                               Conditions = paste(rep("WT",9),rep(c("ctr","4Gy4h","10Gy4h"),each=3),sep = "_"),
                               Labels = paste(Conditions,rep(c("1","2","3"),3),sep="_"),
                               p.cutoff = 0.1,
                               LogFC.cutoff = 1,
                               to.test = "WT_Cvs4Gy"){
    require(magrittr)
    files <- list.files(path,pattern=pattern,full.names = T)
    if(!is.null(exclude)){
        files <- files[-grep(exclude,files)]
    }
    SampleGroup <- factor(Conditions)

    DG <- edgeR::readDGE(files,header= FALSE,group = SampleGroup,labels =Labels)
    filter <- HTSFilter::HTSFilter(DG$counts,DG$samples$group,plot=F)
    DG$counts <- filter$filteredData
    DG <- edgeR::calcNormFactors(DG, method="TMM")
    designMat <- model.matrix(~0+group,data=DG$samples)
    colnames(designMat) <- levels(DG$samples$group)
    DG <- edgeR::estimateGLMCommonDisp(DG, design=designMat)
    DG <- edgeR::estimateGLMTrendedDisp(DG, design=designMat)
    DG <- edgeR::estimateGLMTagwiseDisp(DG, design=designMat)
    fit <- edgeR::glmFit(DG, designMat)
    my.contrasts <- makeContrasts(
        WT_Cvs4Gy = WT_4Gy4h-WT_ctr,
        WT_Cvs10Gy = WT_10Gy4h-WT_ctr,
        WT_10Gyvs4Gy = WT_10Gy4h-WT_4Gy4h,
        levels=designMat)
    lrt <- edgeR::glmLRT(fit, contrast=my.contrasts[,to.test])

    table.counts <- lrt$table %>%
        tibble::rownames_to_column() %>%
        dplyr::mutate(p.adj=p.adjust(PValue,method="BH")) %>%
        dplyr::mutate(FILTER.FC = ifelse(abs(logFC)>LogFC.cutoff,1,0)) %>%
        dplyr::mutate(FILTER.P = ifelse(p.adj<p.cutoff,1,0))
    return(table.counts)
}
