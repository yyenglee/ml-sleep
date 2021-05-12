#pdir = "E:/Lab/sleep/sleep_gene/workDIR/feature_selection/"
sGene = "E:/Lab/sleep/sleep_gene/data/sleep_regulating_gene.xlsx"

geneAliasDict <- list("mouse" = paste0(pdir,"/REFERENCE/Mus_musculus.gene_info.gz"),
                      "drosophila" = paste0(pdir,"/REFERENCE/Drosophila_melanogaster.gene_info.gz"),
                      "human" = paste0(pdir,"/REFERENCE/Homo_sapiens.gene_info.gz"))

taxid <- list("mouse"= "10090", 
              "drosophila" ="7227",
              "human"="9606")
