##################################################################################
### The scripts used for generating every figures in the Gianluca's manuscript ###
##################################################################################

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(Matrix)
library(gplots)

celltype_col_plate = c("Paraxial mesoderm" = "#A65628",
                       "Neuroectoderm" = "#D9D9D9",
                       "Amniotic mesoderm" = "#E41A1C",
                       "Endothelium" = "#FB9A99",
                       "Parietal endoderm" = "#FDBF6F",
                       "Early development" = "#7570B3",
                       "Cardiac mesoderm" = "#F1E2CC",
                       "Surface ectoderm" = "#FFFF99",
                       "Neuromesodermal progenitors" = "#BF5B17",
                       "Extraembryonic ectoderm" = "#A6CEE3",
                       "Allantois" = "#FFFF99",
                       "Visceral endoderm" = "#B3CDE3",
                       "Primitive erythroid cells" = "#4DAF4A",
                       "Gut" = "#F4CAE4",
                       "Hematoendothelial progenitors" = "#F781BF",
                       "Neural crest" = "#CCEBC5",
                       "Megakaryocytes" = "#8DD3C7",
                       "Heart field" = "#FDB462",
                       "White blood cells" = "#F0027F")

day_col_plate = c("E7.5" = "#440154",
                  "E8" = "#3b528b",
                  "E8.5" = "#21918c",
                  "E8.75" = "#5ec962",
                  "E9.5" = "#0072B2",
                  "Day6" = "#D55E00",
                  "Day8" = "#F0027F",
                  "other"= "#E5E5E5")

##############
### Fig.1f ###
##############

pd = read.csv("global_UMAP.csv", header=T, row.names=1, as.is=T)

p = pd[sample(1:nrow(pd),100000),] %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    scale_color_manual(values=celltype_col_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Fig.1f.png", width = 8, height = 8, dpi = 300)


##############
### Fig.1g ###
##############

pd = read.csv("Fig1g_data.csv", header=T, as.is=T)

pd = pd[!pd$sample_ID %in% c("ETiX-Day6-failed-1", "ETiX-Day6-failed-2",
                             "ETiX-Day8-failed-1", "ETiX-Day8-failed-2", "ETiX-Day8-failed-3", "ETiX-Day8-failed-4",
                             "ETiX-Day8-Pax6neg-1", "ETiX-Day8-Pax6neg-2"),]

pd$sample_ID = factor(pd$sample_ID, levels = rev(c("NE-E7.5-1","NE-E7.5-2","NE-E7.5-3",
                                                   "ETiX-Day6-1","ETiX-Day6-2","ETiX-Day6-3",
                                                   "NE-E8","NE-E8.5","NE-E8.75","NE-E9.5-1","NE-E9.5-2",
                                                   "ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5")))

pd$celltype = factor(pd$celltype, levels = rev(c("Early development",
                                           "Hematoendothelial progenitors",
                                           "Extraembryonic ectoderm",
                                           "Allantois",
                                           "Amniotic mesoderm",
                                           "Visceral endoderm",
                                           "Parietal endoderm",
                                           "Gut",
                                           "Surface ectoderm",
                                           "Neuroectoderm",
                                           "Neural crest",
                                           "Neuromesodermal progenitors",
                                           "Paraxial mesoderm",
                                           "Cardiac mesoderm",
                                           "Heart field",
                                           "Primitive erythroid cells",
                                           "Megakaryocytes",
                                           "White blood cells",
                                           "Endothelium")))

p = pd %>%
    ggplot(aes(fill=celltype, y=cell_number, x=sample_ID)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=celltype_col_plate) +
    labs(x="Individual embryos", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))  +
    ggsave("Fig.1g.pdf",
           dpi = 300,
           height  = 10, 
           width = 10)

##############
### Fig.1h ###
##############

pd = read.csv("global_UMAP.csv", header=T, row.names=1, as.is=T)
df = pd[sample(1:nrow(pd),100000),]

for(i in names(day_col_plate)){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$stage != i | df$condition != "WT"] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Fig.1h.", i,".png"),
                   dpi = 300,
                   height  = 6, 
                   width = 6), silent = TRUE)
}

##############
### Fig.1i ###
##############

dat = read.csv("Fig1i_data.csv", header=T, row.names=1, as.is=T)
### in this Matrix, each row refers to a cell type of natural embryo, 
### each column refers to a cell type of ETiX embryo

pdf("Fig.1i.pdf",12,12)
heatmap.2(as.matrix(t(dat)), 
          col=viridis, 
          scale="col", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.2, 
          cexCol = 0.2,
          margins = c(15,15))
dev.off()

##################
### Ext.Fig.3a ###
##################

pd = read.csv("Ext_Fig3a_data.csv", header=T, row.names=1, as.is=T)

p1 = ggplot(pd, aes(EXON_pct)) + 
    geom_histogram(binwidth = 0.1) + 
    geom_vline(xintercept = 85) +
    theme_classic(base_size = 20) +
    labs(x = "% of reads mapping to exon", y = "Cell number") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.3a_1.png",
           dpi = 300,
           height  = 6, 
           width = 6)

pd_sub = pd[pd$EXON_pct <= 85,]

x1 = mean(pd_sub$log2_umi) - sd(pd_sub$log2_umi)
x2 = mean(pd_sub$log2_umi) + 2*sd(pd_sub$log2_umi)
hist(pd_sub$log2_umi, 500); abline(v = x1); abline(v = x2)

pd_sub = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <=85,]

p2 = ggplot(pd_sub, aes(log2_umi)) + 
    geom_histogram(binwidth = 0.05) + 
    geom_vline(xintercept = x1, color = "red") +
    geom_vline(xintercept = x2, color = "red") +
    theme_classic(base_size = 20) +
    labs(x = "Log2(UMI count per cell)", y = "Cell number") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.3a_2.png",
           dpi = 300,
           height  = 6, 
           width = 6)



##################
### Ext.Fig.3b ###
##################

pd = read.csv("Ext_Fig3b_data.csv", header=T, row.names=1, as.is=T)

p1 = ggplot(pd, aes(EXON_pct)) + 
    geom_histogram(binwidth = 0.1) + 
    geom_vline(xintercept = 85) +
    theme_classic(base_size = 20) +
    labs(x = "% of reads mapping to exon", y = "Cell number") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.3b_1.png",
           dpi = 300,
           height  = 6, 
           width = 6)

pd_sub = pd[pd$EXON_pct <= 85,]

x1 = mean(pd_sub$log2_umi) - sd(pd_sub$log2_umi)
x2 = mean(pd_sub$log2_umi) + 2*sd(pd_sub$log2_umi)
hist(pd_sub$log2_umi, 500); abline(v = x1); abline(v = x2)

pd_sub = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & pd$EXON_pct <=85,]

p2 = ggplot(pd_sub, aes(log2_umi)) + 
    geom_histogram(binwidth = 0.05) + 
    geom_vline(xintercept = x1, color = "red") +
    geom_vline(xintercept = x2, color = "red") +
    theme_classic(base_size = 20) +
    labs(x = "Log2(UMI count per cell)", y = "Cell number") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.3b_2.png",
           dpi = 300,
           height  = 6, 
           width = 6)


##################
### Ext.Fig.3c ###
##################

pd = read.csv("global_UMAP.csv", header=T, row.names=1, as.is=T)

col_plate = c("nextseq_1" = "#E69F00",
              "nextseq_2" = "#0072B2",
              "other"= "#E5E5E5")

for(i in c("nextseq_1","nextseq_2")){
    
    df = pd[sample(1:nrow(pd),100000),]
    
    df$Anno1 = as.vector(df$batch)
    df$Anno1[df$batch != i] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=col_plate) +
            ggsave(paste0("Ext.Fig.3c.", i,".png"),
                   dpi = 300,
                   height  = 6, 
                   width = 6), silent = TRUE)
}


##################
### Ext.Fig.3d ###
##################

pd = read.csv("global_UMAP.csv", header=T, row.names=1, as.is=T)

col_plate = c("nextseq_1" = "#E69F00",
              "nextseq_2" = "#0072B2")

pd$sample_ID = factor(pd$sample_ID, levels = rev(c("NE-E7.5-1","NE-E7.5-2","NE-E7.5-3","NE-E8","NE-E8.5","NE-E8.75","NE-E9.5-1","NE-E9.5-2",
                                                   "ETiX-Day6-1","ETiX-Day6-2","ETiX-Day6-3","ETiX-Day6-failed-1","ETiX-Day6-failed-2",
                                                   "ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5",
                                                   "ETiX-Day8-failed-1","ETiX-Day8-failed-2","ETiX-Day8-failed-3","ETiX-Day8-failed-4","ETiX-Day8-Pax6neg-1","ETiX-Day8-Pax6neg-2")))

pd %>% group_by(sample_ID, batch) %>%
    tally() %>%
    ggplot(aes(x = sample_ID, y = n, fill = batch)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values=col_plate) +
    theme_classic(base_size = 20) +
    labs(x = "Individual natural or synthetic embryos", y = "Cell number") +
    coord_flip() +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.3d.pdf",
           dpi = 300,
           height  = 10, 
           width = 10)

##################
### Ext.Fig.3e ###
##################

pd = read.csv("global_UMAP_NE_only.csv", header=T, row.names=1, as.is=T)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    scale_color_manual(values=celltype_col_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Ext.Fig.3e.png", width = 8, height = 8, dpi = 300)



##################
### Ext.Fig.3f ###
##################

pd = read.csv("global_UMAP_ETiX_only.csv", header=T, row.names=1, as.is=T)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    scale_color_manual(values=celltype_col_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Ext.Fig.3f.png", width = 8, height = 8, dpi = 300)


##################
### Ext.Fig.3g ###
##################

df = read.csv("global_UMAP.csv", header=T, row.names=1, as.is=T)

sample_list = c("ETiX-Day6-1","ETiX-Day6-2","ETiX-Day6-3",
                "ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5")

for(i in sample_list){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$sample_ID != i] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Ext.Fig.3g", i,".png"),
                   dpi = 300,
                   height  = 6, 
                   width = 6), silent = TRUE)
}



##################
### Ext.Fig.3h ###
##################

df = read.csv("Ext_Fig3h_data.csv", header=T, row.names=1, as.is=T)

df$sample_resource = factor(df$sample_resource, levels = c("NE","ETiX"))

p = ggplot(df, aes(PC_1,PC_2, color=stage, shape = sample_resource, label=sample_ID)) + 
    geom_point(size=2) + 
    geom_text(size=3, hjust = -0.1) +
    scale_color_manual(values=day_col_plate) +
    theme_classic(base_size = 12) +
    labs(x = "PC_1 (26.4%)", y = "PC_2 (13.6%)") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.3h.pdf",
           dpi = 300,
           height  = 6, 
           width = 8)

####################
### Ext.Fig.4a-c ###
####################

pd = read.csv("Integration_three_datasets.csv", header=T, row.names=1, as.is=T)
df = pd[sample(1:nrow(pd),100000),]

sample_list = names(table(df$group))

for(i in sample_list){
    
    df$Anno1 = as.vector(df$celltype)
    df$Anno1[df$group != i] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.1) +
            theme_void() +
            theme(legend.position="none") +
            ggsave(paste0("Ext.Fig.4", i,".png"),
                   dpi = 300,
                   height  = 6, 
                   width = 6), silent = TRUE)
}

####################
### Ext.Fig.4d-e ###
####################

pd = read.csv("Fig1g_data.csv", header=T, as.is=T)

pd = pd[pd$sample_ID %in% c("ETiX-Day6-1","ETiX-Day6-2","ETiX-Day6-3",
                            "ETiX-Day6-failed-1", "ETiX-Day6-failed-2",
                            "ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5",
                            "ETiX-Day8-failed-1", "ETiX-Day8-failed-2", "ETiX-Day8-failed-3", "ETiX-Day8-failed-4"),]

pd$sample_ID = factor(pd$sample_ID, levels = rev(c("ETiX-Day6-1","ETiX-Day6-2","ETiX-Day6-3",
                                                   "ETiX-Day6-failed-1", "ETiX-Day6-failed-2",
                                                   "ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5",
                                                   "ETiX-Day8-failed-1", "ETiX-Day8-failed-2", "ETiX-Day8-failed-3", "ETiX-Day8-failed-4")))

pd$celltype = factor(pd$celltype, levels = rev(c("Early development",
                                                 "Hematoendothelial progenitors",
                                                 "Extraembryonic ectoderm",
                                                 "Allantois",
                                                 "Amniotic mesoderm",
                                                 "Visceral endoderm",
                                                 "Parietal endoderm",
                                                 "Gut",
                                                 "Surface ectoderm",
                                                 "Neuroectoderm",
                                                 "Neural crest",
                                                 "Neuromesodermal progenitors",
                                                 "Paraxial mesoderm",
                                                 "Cardiac mesoderm",
                                                 "Heart field",
                                                 "Primitive erythroid cells",
                                                 "Megakaryocytes",
                                                 "White blood cells",
                                                 "Endothelium")))

p = pd %>%
    ggplot(aes(fill=celltype, y=cell_number, x=sample_ID)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=celltype_col_plate) +
    labs(x="Individual embryos", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))  +
    ggsave("Ext.Fig.4d.pdf",
           dpi = 300,
           height  = 10, 
           width = 10)

pd$condition = if_else(pd$sample_ID %in% c("ETiX-Day6-1","ETiX-Day6-2","ETiX-Day6-3",
                                           "ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5"), "WT", "Failed")
pd$cell_fraction_per_embryo  = pd$cell_fraction_per_embryo * 100

col_plate = c("WT" = "#D55E00",
              "Failed" = "#009E73")
p = ggplot(pd[pd$sample_ID %in% c("ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5",
                                   "ETiX-Day8-failed-1", "ETiX-Day8-failed-2", "ETiX-Day8-failed-3", "ETiX-Day8-failed-4"),], aes(celltype, cell_fraction_per_embryo, fill = condition)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1)) +
    scale_fill_manual(values = col_plate) +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black")) +
    ggsave("Ext.Fig.4e.pdf",
           dpi = 300,
           height  = 8, 
           width = 10) 




################
### Fig.2g,h ###
################

pd = read.csv("Neuroectoderm_UMAP.csv", header=T, row.names=1, as.is=T)

neuroectoderm_celltype_col_plate = c("Hindbrain & Spinal cord" = "#E69F00",
                       "Prosencephalon" = "#56B4E9",
                       "Mesencephalon & MHB" = "#009E73",
                       "Floor plate" = "#0072B2",
                       "Roof plate" = "#CC79A7",
                       "Early neurons" = "#F0027F")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    scale_color_manual(values=neuroectoderm_celltype_col_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Fig.2g.png", width = 8, height = 8, dpi = 300)


for(i in names(day_col_plate)){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$stage != i | df$condition != "WT"] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Fig.2h.", i,".png"),
                   dpi = 300,
                   height  = 3, 
                   width = 3.5), silent = TRUE)
}


##############
### Fig.2i ###
##############

pd = read.csv("Fig2i_data.csv", header=T, as.is=T)

pd = pd[pd$sample_ID %in% c("ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5"),]

pd$sample_ID = factor(pd$sample_ID, levels = rev(c("ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5")))

p = pd %>%
    ggplot(aes(fill=celltype, y=cell_number, x=sample_ID)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=celltype_col_plate) +
    labs(x="Individual embryos", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))  +
    ggsave("Fig.2i.pdf",
           dpi = 300,
           height  = 10, 
           width = 10)


####################
### Ext.Fig.7a-e ###
####################

cds = readRDS("cds_Neuroectoderm.rds")

### for example, Pax6
gene_list = c("Pax6")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)

##################
### Ext.Fig.7f ###
##################

cds = readRDS("cds.rds")

### for example, Sox10
gene_list = c("Sox10")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)


#################
### Ext.Fig.8 ###
#################

pd = read.csv("Fig1g_data.csv", header=T, as.is=T)

pd = pd[pd$sample_ID %in% c("ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5",
                            "ETiX-Day8-Pax6neg-1","ETiX-Day8-Pax6neg-2"),]

pd$sample_ID = factor(pd$sample_ID, levels = rev(c("ETiX-Day8-1","ETiX-Day8-2","ETiX-Day8-3","ETiX-Day8-4","ETiX-Day8-5",
                                                   "ETiX-Day8-Pax6neg-1","ETiX-Day8-Pax6neg-2")))

pd$celltype = factor(pd$celltype, levels = rev(c("Early development",
                                                 "Hematoendothelial progenitors",
                                                 "Extraembryonic ectoderm",
                                                 "Allantois",
                                                 "Amniotic mesoderm",
                                                 "Visceral endoderm",
                                                 "Parietal endoderm",
                                                 "Gut",
                                                 "Surface ectoderm",
                                                 "Neuroectoderm",
                                                 "Neural crest",
                                                 "Neuromesodermal progenitors",
                                                 "Paraxial mesoderm",
                                                 "Cardiac mesoderm",
                                                 "Heart field",
                                                 "Primitive erythroid cells",
                                                 "Megakaryocytes",
                                                 "White blood cells",
                                                 "Endothelium")))

p = pd %>%
    ggplot(aes(fill=celltype, y=cell_number, x=sample_ID)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=celltype_col_plate) +
    labs(x="Individual embryos", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))  +
    ggsave("Fig.2i.pdf",
           dpi = 300,
           height  = 10, 
           width = 10)



################
### Fig.3n,o ###
################

pd = read.csv("Heart_UMAP.csv", header=T, row.names=1, as.is=T)

heart_celltype_col_plate = c("Cardiac mesoderm" = "#E69F00",
                       "First heart field" = "#56B4E9",
                       "Second heart field" = "#F0027F")

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    scale_color_manual(values=heart_celltype_col_plate) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Fig.3n.png", width = 8, height = 8, dpi = 300)

for(i in names(day_col_plate)){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$stage != i | df$condition != "WT"] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Fig.3o.", i,".png"),
                   dpi = 300,
                   height  = 3, 
                   width = 3.5), silent = TRUE)
}


###################
### Ext.Fig.10c ###
###################

cds = readRDS("cds.rds")

### for example, Meox2
gene_list = c("Meox2")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)

###################
### Ext.Fig.10d ###
###################

cds = readRDS("cds_Heart.rds")

### for example, Hand1
gene_list = c("Hand1")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)


################
### Fig.4d,e ###
################

pd = read.csv("Endoderm_UMAP.csv", header=T, row.names=1, as.is=T)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Fig.4d.png", width = 8, height = 8, dpi = 300)

for(i in names(day_col_plate)){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$stage != i | df$condition != "WT"] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Fig.4e.", i,".png"),
                   dpi = 300,
                   height  = 3, 
                   width = 3.5), silent = TRUE)
}


##################
### Ext.Fig.12 ###
##################

cds = readRDS("cds_Gut.rds")

### for example, Hand1
gene_list = c("Apela")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)



##############
### Fig.5d ###
##############

pd = read.csv("ExE_endoderm_UMAP.csv", header=T, row.names=1, as.is=T)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Fig.5d.png", width = 8, height = 8, dpi = 300)


################
### Fig.5g,h ###
################

pd = read.csv("ExE_ectoderm_UMAP.csv", header=T, row.names=1, as.is=T)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Fig.5g.png", width = 8, height = 8, dpi = 300)

for(i in names(day_col_plate)){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$stage != i | df$condition != "WT"] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Fig.5h.", i,".png"),
                   dpi = 300,
                   height  = 3, 
                   width = 3.5), silent = TRUE)
}



##############
### Fig.5i ###
##############

cds = readRDS("cds_ExE_ectoderm.rds")

### for example, Prl4a1
gene_list = c("Prl4a1")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)






###################
### Ext.Fig.14c ###
###################

pd = read.csv("ExE_endoderm_UMAP.csv", header=T, row.names=1, as.is=T)

p = pd %>%
    ggplot() +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.6) +
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
    theme_void() +
    theme(legend.position="none") + 
    ggsave("Ext.Fig.14c.png", width = 8, height = 8, dpi = 300)

for(i in names(day_col_plate)){
    
    df$Anno1 = as.vector(df$stage)
    df$Anno1[df$stage != i | df$condition != "WT"] = "other"
    
    try(ggplot(df) +
            geom_point(aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            geom_point(data = subset(df, Anno1 != 'other'),
                       aes(x = UMAP_1, y = UMAP_2, color = Anno1), size=0.3) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=day_col_plate) +
            ggsave(paste0("Ext.Fig.14c.", i,".png"),
                   dpi = 300,
                   height  = 3, 
                   width = 3.5), silent = TRUE)
}



###################
### Ext.Fig.14d ###
###################

cds = readRDS("cds_ExE_endooderm.rds")

### for example, Prl4a1
gene_list = c("Prl4a1")

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "NE"], genes = gene_list)

plot_cells(cds[,cds$condition == "WT" & cds$sample_resource == "ETiX"], genes = gene_list)



