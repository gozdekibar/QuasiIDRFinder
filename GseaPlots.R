library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library("readxl")


my_data_metapredict <- read_excel("./Dataset/MasterTable_canonicalHuman_qIDR.xlsx")

overlapped_names=my_data_metapredict[my_data_metapredict$overlap_metapredict == "yes", ]
x=data.frame(overlapped_names$gene_symbol,overlapped_names$`min(p-val)`)
x=x[!duplicated(x$overlapped_names.gene_symbol), ]

colnames(x)=c("Gene", "metric")
newdata <- x[order(x$metric),] 
newdata$metric=-log10(newdata$metric)


## Extract the foldchanges
foldchanges <- newdata$metric

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- newdata$Gene

foldchanges <- sort(foldchanges, decreasing = TRUE)
head(foldchanges)



# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = foldchanges, 
                OrgDb = org.Hs.eg.db, 
                ont = 'BP', 
                nPerm = 1000, 
                minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE,keyType = 'SYMBOL') 


gseaGO_results <- gseaGO@result


# GSEA analysis for Rna Binding Proteins 
c2 <- read.gmt("./Dataset/gmt_files/rna_binding_all.gmt.txt")

msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)

gseaGO_results <- msig@result

#pdf("plot_gsea_rna_binding_proteins.pdf", width =5 , height = 5.5)
gseaplot(msig, geneSetID = 'rna binding proteins',title = msig$Description[1])
#dev.off()

#pdf("plot_gsea_rbp_2.pdf", width =5 , height = 5)
gseaplot2(msig, geneSetID = 1, title = msig$Description[1])
#dev.off()


###GO analysis
c2 <- read.gmt("./Dataset/gmt_files/PLDs.gmt.txt")

msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)

#pdf("plot_gsea_plds_1.pdf", width =5 , height = 5.5)
gseaplot(msig, geneSetID = 'PLDs',title = msig$Description[1])
#dev.off()


#pdf("plot_gsea_plds_2.pdf", width =5 , height = 5)
gseaplot2(msig, geneSetID = 1, title = msig$Description[1])
#dev.off()

####

c2 <- read.gmt("./Dataset/gmt_files/tfs.gmt.txt")

msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE,maxGSSize = 700, pvalueCutoff = 1)
gseaGO_results <- msig@result

#pdf("plot_tfs_tfs_1.pdf", width =5 , height = 5.5)
gseaplot(msig, geneSetID = 'tfs',title = msig$Description[1])
#dev.off()

#pdf("plot_gsea_tfs_2.pdf", width =5 , height = 5)
gseaplot2(msig, geneSetID = 1, title = msig$Description[1])
#dev.off()

