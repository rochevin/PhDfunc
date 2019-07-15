GenerateDE <- function(out = "",
                       my.ex="DiVA",
                       my.pval = 0.1,
                       my.FC = 0.5,
                       output_file = paste0(out,"AutoReport_DE_edgeR-", my.ex, "_", my.pval,"_",my.FC, ".html")){
    my.experiments <- list(
        "DiVA" = list(
            my.metadata.file = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/results/DE/RNA-SEQ_DIVA_LEGUBE/autoReport/metadata.tsv"
            ,my.files = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Legube/PROCESSED_strand-spe-reverse_10102018/mapping/counts"
            ,DE.1 = "pOHT"
            ,DE.2 = "mOHT"
            ,experiment = "RNA-Seq genes counts on DiVA cells before and after breaks"
            ,gene.list = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/biomart_export_gene_id_name.txt"
        ),
        "IR_4Gy" = list(
            my.metadata.file = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/results/DE/RNA-SEQ_IR_GSE110386/autoReport/metadata.tsv"
            ,my.files = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_GSE110386/counts"
            ,DE.1 = "WT_4Gy4h"
            ,DE.2 = "WT_ctr"
            ,experiment = "RNA-Seq genes counts on U2OS cells before and after breaks 4H after 4Gy of ionizing radiations"
            ,gene.list = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/biomart_export_gene_id_name.txt"
        ),
        "IR_10Gy" = list(
            my.metadata.file = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/results/DE/RNA-SEQ_IR_GSE110386/autoReport/metadata.tsv"
            ,my.files = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_GSE110386/counts"
            ,DE.1 = "WT_10Gy4h"
            ,DE.2 = "WT_ctr"
            ,experiment = "RNA-Seq genes counts on U2OS cells before and after breaks 4H after 10Gy of ionizing radiations"
            ,gene.list = "/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/biomart_export_gene_id_name.txt"
        )
    )
    rmarkdown::render(
        "templates/template_DE_edgeR.Rmd", params = list(
            pval.cutoff = my.pval,
            FC.cutoff = my.FC,
            exp=my.ex,
            metadata =  my.experiments[[my.ex]][["my.metadata.file"]],
            DE.1 = my.experiments[[my.ex]][["DE.1"]],
            DE.2 = my.experiments[[my.ex]][["DE.2"]],
            experiment = my.experiments[[my.ex]][["experiment"]],
            gene.list = my.experiments[[my.ex]][["gene.list"]],
            out = out
        ),
        output_file = output_file
    )
    output_tsv <- paste0(out,paste("DE_result_edgeR",my.ex,LogFC.cutoff,p.cutoff,my.experiments[[my.ex]][["DE.1"]],my.experiments[[my.ex]][["DE.2"]],"table.tsv",sep="|"))
    data <- read_tsv(output_tsv)
    data
}
