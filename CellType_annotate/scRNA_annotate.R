
dir <- '/project/lidong/scRNA/text/GSE166555'
object <- 'All' # All & other

source('/nas/bin/annotate')
data <- prepare_data(dir,object);seob <- data[[1]];annotate <- data[[2]]
#########################################################################################################
## 调整注释内容
gene <- list(T_cell = c('CD3D','CD3E','CD2','CD3G'),
             NK_cell = c('NKG7','CTSW','PRF1','GNLY'),
             B_cell = c('CD79A','CD79B','IGHM','CD19'),
             Macrophages = c('CD68','CD163','CD14','MS4A7'),
             DC_cell = c('HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1'),
             #DC_cell = c('CD83','CD86','HCK','AIF1'),
             Epithelial = c('EPCAM','ELF3','TSTD1','TSPAN8'),
             Endothelial = c('CDH5','CLDN5','EGFL7','PECAM1'),
             Fibroblasts = c('DCN','COL1A1','COL1A2','COL3A1'),
             Mast_cell = c('CPA3','TPSAB1','TPSB2','IL1RL1')
)

gene <- geneshow(seob,gene)
dev.off()

#有0组存在，每个位置对应星状图加1
annotate$labels <- as.numeric(annotate$seurat_clusters)-1
annotate$labels[c(1,5)+1] <- "T_cell"
annotate$labels[c(16,19,21)+1] <- "NK_cell"
annotate$labels[c(8,18,20)+1] <- "B_cell"
annotate$labels[c(0)+1] <- "Plasma_B"
annotate$labels[c(10)+1] <- "Macrophages"
annotate$labels[c(17)+1] <- "DC_cell"
annotate$labels[c(2,3,4,6,7,9,12)+1] <- "Epithelial"
annotate$labels[c(14)+1] <- "Endothelial"
annotate$labels[c(11,13)+1] <- "Fibroblasts"
annotate$labels[c(15)+1] <- "Mast_cell"


##############################################################################################
main(seob,object,col2,gene,annotate,dir)

