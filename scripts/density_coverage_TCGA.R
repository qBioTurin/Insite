library(TCGAbiolinks)
library(dplyr)
library(tidyr)

Tumor.purity$purity<-as.numeric(stringr::str_replace_all(pattern = ",",replacement = ".",string = Tumor.purity$CPE))
Tumor.purity%>%
  filter(purity>0.7)%>%
  pull(Sample.ID)%>%as.vector()->High_pur_sample_ID

query_maf <- GDCquery(
  project = c("TCGA-COAD",
              "TCGA-DLBC",
              "TCGA-ESCA",
              "TCGA-GBM",
              "TCGA-HNSC",
              "TCGA-KICH",
              "TCGA-KIRC",
              "TCGA-KIRP",
              "TCGA-LAML",
              "TCGA-LGG",
              "TCGA-LIHC",
              "TCGA-LUAD",
              "TCGA-LUSC",
              "TCGA-MESO",
              "TCGA-OV",
              "TCGA-PAAD",
              "TCGA-PCPG",
              "TCGA-PRAD",
              "TCGA-READ",
              "TCGA-SARC",
              "TCGA-SKCM",
              "TCGA-STAD",
              "TCGA-TGCT",
              "TCGA-THCA",
              "TCGA-THYM",
              "TCGA-UCEC",
              "TCGA-UCS",
              "TCGA-UVM"),
  data.category = "Simple Nucleotide Variation",
  access = "open",
  barcode = High_pur_sample_ID,
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_maf,files.per.chunk=10)
maf <- GDCprepare(query_maf)

dens_coverage<-density(maf$t_depth)
save(dens_coverage,file = "Data/dens.RData")
