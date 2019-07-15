GetDE_ETO_GSE120612 <- function(path = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/ETO_GSE120612/counts",
                                pattern ="_counts.tsv",
                                exclude = "STX",
                                Conditions = c("C1","C1_OHT"),
                                Labels = paste(Conditions,rep(c("1","2"),each=2),sep="_"),
                                p.cutoff = 0.1,
                                LogFC.cutoff = 0.5,
                                gene_id_name = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/biomart_export_gene_id_name.txt",
                                README = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/ETO_GSE120612/counts/README"){
    require(tidyverse)
    require(edgeR)
    README <- README %>% read_tsv(col_names = F)
    names(README) <- c("SRA","GSM","Sample")
    README <- README %>% mutate(Cond = Sample) %>% separate(Cond,into = c("Cell","Condition","Replicate")) %>% mutate(group = paste(Condition,Cell,sep="_")) %>% mutate_at(c("Cell","Condition","Replicate","group"),as_factor)
    files <- list.files(path,pattern="_counts.tsv",full.names = T)
    SampleGroup <- factor(README$group)
    DG <- readDGE(files,header= FALSE,group = SampleGroup,labels =README$Sample)
    filter <- HTSFilter::HTSFilter(DG$counts,DG$samples$group,plot=T)
    DG$counts <- filter$filteredData
    DG <- calcNormFactors(DG, method="TMM")
    designMat <- model.matrix(~0+group,data=DG$samples)
    colnames(designMat) <- levels(DG$samples$group)

    DG <- estimateGLMCommonDisp(DG, design=designMat)
    DG <- estimateGLMTrendedDisp(DG, design=designMat)
    DG <- estimateGLMTagwiseDisp(DG, design=designMat)

    fit <- glmFit(DG, designMat)

    my.contrasts <- makeContrasts(
        NALM6 = Etoposide_NALM6-Control_NALM6,
        C697 = Etoposide_697-Control_697,
        RAMOS = Etoposide_RAMOS-Control_RAMOS,
        JURKAT = Etoposide_JURKAT-Control_JURKAT,
        THP1 = Etoposide_THP1-Control_THP1,
        U937 = Etoposide_U937-Control_U937,
        levels=designMat)

    genes.att <- read.delim(gene_id_name,sep="\t",h=T)


    my.Eto.vals <- lapply(colnames(my.contrasts),function(to.test){
        lrt <- glmLRT(fit, contrast=my.contrasts[,to.test])
        lrt$table %>%
            tibble::rownames_to_column() %>%
            mutate(p.adj=p.adjust(PValue,method="BH")) %>%
            mutate(FILTER.FC = ifelse(abs(logFC)>LogFC.cutoff,1,0)) %>%
            mutate(FILTER.P = ifelse(p.adj<p.cutoff,1,0)) %>%
            mutate(group = ifelse(logFC>0,"Upregulated","Downregulated")) %>%
            dplyr::left_join(genes.att,by = c("rowname" = "Gene.stable.ID"))
    })
    names(my.Eto.vals) <- colnames(my.contrasts)
    my.Eto.vals
}
