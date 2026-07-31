
# 
# 1 - Setup, grouping & data import
# 2 - PCA plots
# 3 - get biomart data
# 4 - limma exploration
# 5 - Sleuth across genes initial

########## 1 Setup, grouping & data import ###################
#if (!require("BiocManager", quietly = T)) { 
#    install.packages("BiocManager") }

#install.packages("aggregation") # for genes from tx
library(aggregation)
#install.packages("devtools")
library(devtools)
#BiocManager::install("rhdf5")
library(rhdf5)
#devtools::install_github("pachterlab/sleuth")
library(sleuth) #vignette('intro', package = 'sleuth')
library(ggplot2) # v3.5.2
library(ggrepel)
library(dplyr)
library(tidyr)
library(dbplyr)
library(ggpubr)
library(grid) # v4.4.1
library(gridExtra)
#BiocManager::install("limma")
library(limma) # v3.60.6
#BiocManager::install("edgeR")
library(edgeR)
#BiocManager::install("tximport")
library(tximport)
library(readr)
library(ggvenn) #
library(ggVennDiagram) # v1.5.2
library(VennDiagram)
# BiocManager::install(c("biomaRt", "BiocFileCache"), update = TRUE)
# BiocManager::install("biomaRt")
library(biomaRt)

# [1] Load a table describing our sample, conditions
# and the source directories such that the 1st column
# contains the sample names, the middle column(s)
# contain the conditions, and the last column has the
# folder containing the kallisto output

s2c <- read.csv("LIMMA/table.csv", header=T, sep="\t")
str(s2c) # check table 

# reallocate table into groups: 
s2c_c <- subset(s2c, condition=="0") # 
s2c_3 <- subset(s2c, condition=="3") # 
s2c_7 <- subset(s2c, condition=="7") # 
s2c_10 <- subset(s2c, condition=="10") # 
s2c_14 <- subset(s2c, condition=="14") # 
s2c_21 <- subset(s2c, condition=="21") # 
s2c_28 <- subset(s2c, condition=="28") # 
s2c_35 <- subset(s2c, condition=="35") # 

# Now, we will use sleuth_prep() to make an object 
# with our experiment info, model and groups.
so <- sleuth_prep(s2c, extra_bootstrap_summary=T, read_bootstrap_tpm=T)
# XX targets pass filter

# all sample
#pdf("pca.pdf", width=5, height=8)
#plot_pca(so, color_by = 'condition', text_labels=T) # do PCA
#dev.off()

######## 2 PCA plots ##################################

# Get TPM values, filter, make PCA plots

tpm_data <- so$obs_norm %>% # extract TPM values
  dplyr::select(target_id, sample, tpm) %>%
  tidyr::pivot_wider(names_from = sample, values_from = tpm)
# Convert to matrix format for PCA
tpm_matrix <- as.matrix(tpm_data[,-1])
rownames(tpm_matrix) <- tpm_data$target_id
var_per_gene <- apply(tpm_matrix, 1, var) # filter for non-zero values
tpm_matrix_filtered <- tpm_matrix[var_per_gene > 0, ]
pca_result <- prcomp(t(tpm_matrix_filtered), scale = T)
pca_df <- as.data.frame(pca_result$x)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, so$sample_to_covariates, by = "sample")
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
unique_combinations <- unique(pca_df$condition)

library(ggforce) # Required for geom_mark_ellipse
# 1. Define Color Palette
viridis_8 <- hcl.colors(n = 8, palette = "viridis")
base_colours <- c("0" = viridis_8[1], "3" = viridis_8[2], "7" = viridis_8[3],
                 "10" = viridis_8[4], "14" = viridis_8[5], "21" = viridis_8[6],
                 "28" = viridis_8[7], "35" = viridis_8[8])
base_colors <- c("0" = viridis_8[1], "3" = viridis_8[2], "7" = viridis_8[3],
                 "10" = viridis_8[4], "14" = viridis_8[5], "21" = viridis_8[6],
                 "28" = viridis_8[7], "35" = viridis_8[8])
color_palette <- base_colors[as.character(unique_combinations)]
oval_stats <- pca_df %>% group_by(condition) %>%
  summarise(  # PC1 & PC2
    mean_PC1 = mean(PC1), mean_PC2 = mean(PC2),
    sd_PC1   = ifelse(is.na(sd(PC1)) || sd(PC1) == 0, 0.5, sd(PC1) * 1.5),
    sd_PC2   = ifelse(is.na(sd(PC2)) || sd(PC2) == 0, 0.5, sd(PC2) * 1.5),
    mean_PC3 = mean(PC3), mean_PC4 = mean(PC4),
    sd_PC3   = ifelse(is.na(sd(PC3)) || sd(PC3) == 0, 0.5, sd(PC3) * 1.5),
    sd_PC4   = ifelse(is.na(sd(PC4)) || sd(PC4) == 0, 0.5, sd(PC4) * 1.5) )

# Define Color Palette
viridis_8 <- hcl.colors(n = 8, palette = "viridis")
base_colors <- c("0" = viridis_8[1], "3" = viridis_8[2], "7" = viridis_8[3],
                 "10" = viridis_8[4], "14" = viridis_8[5], "21" = viridis_8[6],
                 "28" = viridis_8[7], "35" = viridis_8[8])
color_palette <- base_colors[as.character(unique_combinations)]

# --- PLOT 1: PC1 vs PC2 ---
pdf("LIMMA/PC1.PC2.pdf", width=8, height=5)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(condition), fill = factor(condition))) +
  geom_point(size = 5, alpha = 0.5) +
  geom_text(aes(label = condition), color = "black", size = 4) +
  scale_color_manual(values = color_palette, name = "Condition") +
  scale_fill_manual(values = color_palette, name = "Condition") +
  labs(x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 21))
print(p)
dev.off()
ggsave("LIMMA/PC1.PC2.png", plot = p, width = 8, height = 5, dpi = 300)

pdf("LIMMA/PC3.PC4.pdf", width=8, height=5)
p2 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = factor(condition), fill = factor(condition))) +
  geom_point(size = 5, alpha = 0.5) +
  geom_text(aes(label = condition), color = "black", size = 4) +
  scale_color_manual(values = color_palette, name = "Condition") +
  scale_fill_manual(values = color_palette, name = "Condition") +
  labs(x = paste0("PC3 (", var_explained[3], "%)"), y = paste0("PC4 (", var_explained[4], "%)")) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 21))
print(p2)
dev.off()
ggsave("LIMMA/PC3.PC4.png", plot = p2, width = 8, height = 5, dpi = 300)

pdf("LIMMA/PC1.PC4.pdf", width=8, height=5)
p14 <- ggplot(pca_df, aes(x = PC1, y = PC4, color = factor(condition), fill = factor(condition))) +
  geom_point(size = 5, alpha = 0.5) +
  geom_text(aes(label = condition), color = "black", size = 4) +
  scale_color_manual(values = color_palette, name = "Condition") +
  scale_fill_manual(values = color_palette, name = "Condition") +
  labs(x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC4 (", var_explained[4], "%)")) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 21))
print(p14)
dev.off()
ggsave("LIMMA/PC1.PC4.png", plot = p14, width = 8, height = 5, dpi = 300)


# --- PLOT 4: PC2 vs PC4 ---
pdf("LIMMA/PC2.PC4.pdf", width=8, height=5)
p24 <- ggplot(pca_df, aes(x = PC2, y = PC4, color = factor(condition), fill = factor(condition))) +
   geom_point(size = 5, alpha = 0.5) +
  geom_text(aes(label = condition), color = "black", size = 4) +
  scale_color_manual(values = color_palette, name = "Condition") +
  scale_fill_manual(values = color_palette, name = "Condition") +
  labs(x = paste0("PC2 (", var_explained[2], "%)"),
       y = paste0("PC4 (", var_explained[4], "%)")) +
  theme_bw() + 
  theme(axis.text =element_text(size = 16), axis.title=element_text(size = 20),
        legend.text =element_text(size=20), legend.title=element_text(size=21))
print(p24)
dev.off()
ggsave("LIMMA/PC2.PC4.png", plot = p24, width = 8, height = 5, dpi = 300)

pdf("plot_sample_heatmap.pdf", width=11, height=11)
plot_sample_heatmap(so ,use_filtered = T, color_high = "white",
  color_low = "dodgerblue", x_axis_angle = 50,
  annotation_cols=setdiff(colnames(so$sample_to_covariates), "sample"),
  cluster_bool = T)
dev.off()

####### 3 get biomart data ##########################

# Let's test across genes  
# We need to map the isoforms to genes with BiomaRt
# Note we need to have annotation matching our reference cDNAs

# BiocManager::install("BiocFileCache", update = TRUE, ask = FALSE)
# devtools::install_version("dbplyr", version = "2.3.4")

library(biomaRt)

mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "sscrofa_gene_ensembl",
  host    = "https://sep2025.archive.ensembl.org/"
)

# Add useCache = FALSE to bypass the broken dbplyr/BiocFileCache step
t2g <- getBM(
  attributes = c(
    "ensembl_transcript_id", 
    "transcript_version", 
    "ensembl_gene_id", 
    "external_gene_name",
    "description", 
    "transcript_biotype"
  ), 
  mart = mart,
  useCache = FALSE
)

str(t2g) # 60273 rows of:
        # ensembl_transcript_id
        # ensembl_gene_id
        # external_gene_name
#  rename the transcripts
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
t2g$target_id <- paste(t2g$target_id, t2g$transcript_version,sep=".") # Re-name
str(t2g) # 60,440

############ 4 limma exploration ########################

# transcripts first
#  Examine the TPM data using tximport  
files = paste(s2c$path, "abundance.h5", sep="/")

# Import Kallisto abundance.h5 files with tximport
txi.kallisto <- tximport(files, type = "kallisto", txOut = T)
str(txi.kallisto)
head(txi.kallisto$counts)
y <- DGEList(txi.kallisto$counts)
dim(y) # check
full <- read.csv("LIMMA/table.csv", header=T, sep="\t") # metadata
full$condition <- factor(full$condition)
full$day <- as.numeric(as.character(full$condition))
condition <- factor(full$condition,
                    levels = c(0,3,7,10,14,21,35))

# Treat day as a factor
design <- model.matrix(~0 + condition, data=full)
cont <- makeContrasts(
  D0_vs_D3 = condition0-condition3,
  D0_vs_D7 = condition0-condition7,
  D0_vs_D10 = condition0-condition10,
  D0_vs_D14 = condition0-condition14,
  D0_vs_D21 = condition0-condition21,
  D0_vs_D28 = condition0-condition28,
  D0_vs_D35 = condition0-condition35,
  levels=design ) 
# filtering using the design information:
keep <- filterByExpr(y, design)
y <- y[keep, ]
str(y) # check -> 29,486 transcripts 
y <- calcNormFactors(y) # normalize and run voom transformation
v <- voom(y, design) # v is now ready for lmFit()  

fit <- lmFit(v, design) # eBayes stands for empirical Bayes
fitm <- eBayes(fit, trend=T) # test comparison
str(fitm)
pdf("plotSA.limma.pdf")
plotSA(fitm) # we have a variance trend, so keep trend=T in eBayes()
dev.off()
results <- topTable(fitm, n=dim(fitm$coefficients)[1])

fit2 <- contrasts.fit(fit, cont) # main comparison
fit2 <- eBayes(fit2, trend=T) # get p/w comparison
res_D3 <- topTable(fit2, coef= "D0_vs_D3", number=Inf)
res_D7 <- topTable(fit2, coef= "D0_vs_D7", number=Inf)
res_D10 <- topTable(fit2, coef="D0_vs_D10", number=Inf)
res_D14 <- topTable(fit2, coef="D0_vs_D14", number=Inf)
res_D21 <- topTable(fit2, coef="D0_vs_D21", number=Inf)
res_D28 <- topTable(fit2, coef="D0_vs_D28", number=Inf)
res_D35 <- topTable(fit2, coef="D0_vs_D35", number=Inf) 

write.csv(res_D3, "limma.res_D3.csv") # XXX transcripts 
write.csv(res_D7, "limma.res_D7.csv") # XXX transcripts 
write.csv(res_D10, "limma.res_D10.csv") # XXX transcripts 
write.csv(res_D14, "limma.res_D14.csv") # XXX transcripts 
write.csv(res_D21, "limma.res_D21.csv") # XXX transcripts 
write.csv(res_D28, "limma.res_D28.csv") # XXX transcripts 
write.csv(res_D35, "limma.res_D35.csv") # XXX transcripts 

dim(subset(res_D3, adj.P.Val<=0.05))
dim(subset(res_D7, adj.P.Val<=0.05))
dim(subset(res_D10, adj.P.Val<=0.05))
dim(subset(res_D14, adj.P.Val<=0.05))
dim(subset(res_D21, adj.P.Val<=0.05))
dim(subset(res_D28, adj.P.Val<=0.05))
dim(subset(res_D35, adj.P.Val<=0.05))

summary(res_D3$adj.P.Val)
summary(res_D7$adj.P.Val)
summary(res_D10$adj.P.Val)
summary(res_D14$adj.P.Val)
summary(res_D21$adj.P.Val)
summary(res_D28$adj.P.Val)
summary(res_D35$adj.P.Val)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- abs(cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {  cex.cor <- 1 + 0.4 / strwidth(txt)  }
    text(0.5, 0.5, txt, cex = 1 + cex.cor * Cor) 
    } # Resize the text by level of correlation

pdf("pairs.limma.pdf", width=6, height=6)
pairs( results[,1:4], cex=0.1, upper.panel = panel.cor,
       lower.panel = panel.smooth, cex.labels = 1.5)
dev.off()

# make gene-level volcano plots 

make_volcano <- function(results,  comparison,  t2g,
                         fdr_cutoff = 0.05, logfc_cutoff = 1.39) {
  transcript_ids <- rownames(results)
  if (grepl("\\.", transcript_ids[1])) {
    results$target_id <- rownames(results)
    results_with_genes <- left_join(
      data.frame(  target_id = rownames(results),
        results,  row.names = NULL ),
      t2g %>% dplyr::select(
        target_id,  ens_gene,   ext_gene  ), by = "target_id"  )
  } else { results$ensembl_transcript_id <- rownames(results)
    results_with_genes <- left_join(
      data.frame( ensembl_transcript_id = rownames(results),
        results,  row.names = NULL ),
      t2g %>% dplyr::select(
        ensembl_transcript_id = target_id, ens_gene,
        ext_gene  ), by = "ensembl_transcript_id"  )  }
  results_with_genes$gene_label <- ifelse(
    !is.na(results_with_genes$ext_gene) &
      results_with_genes$ext_gene != "",
    results_with_genes$ext_gene,
    ifelse( !is.na(results_with_genes$ens_gene) &
        results_with_genes$ens_gene != "",
      results_with_genes$ens_gene,
      ifelse( grepl("\\.", transcript_ids[1]),
        results_with_genes$target_id,
        results_with_genes$ensembl_transcript_id)  ))
  # Make infection-induced genes positive
  results_with_genes$logFC <- -results_with_genes$logFC
  results_with_genes$significant <- ifelse(
    results_with_genes$adj.P.Val < fdr_cutoff &
      abs(results_with_genes$logFC) > logfc_cutoff,
    "Significant",   "Not"  )
  top_genes <- results_with_genes %>%
    filter(significant == "Significant")
  cat(  comparison, ":",  nrow(top_genes),  "DE transcripts\n" )
  write.csv(  results_with_genes,
    paste0(comparison, ".all.csv"),   row.names = FALSE  )
  write.csv(  top_genes,
    paste0(comparison, ".DE.csv"), row.names = FALSE )
top_genes2 <- top_genes %>% arrange(adj.P.Val) %>% slice_head(n=10)

  pdf(paste0(comparison, ".volcano.pdf"), width = 6,height =4)
  p <- ggplot(  results_with_genes,
    aes( x = logFC, y = -log10(adj.P.Val),color = significant )) +
    geom_point(  shape = 1, size = 1.2, alpha = 0.2) +
    geom_hline( yintercept = -log10(fdr_cutoff),
      linetype = "dashed", color = "darkgray" ) +
    geom_vline( xintercept = c( -logfc_cutoff, logfc_cutoff),
      linetype = "dashed", color = "darkgray" ) +
    scale_color_manual(values = c(Significant="red",Not ="black" ) ) +
    geom_text_repel(  data = top_genes2,
      aes(label=gene_label),box.padding=0.4,min.segment.length=0,
      max.overlaps =220, force=10, size = 6 ) +
    labs(x = "log2(FC)", y = "-log10(FDR)")  +theme_bw()+ #+ xlim(-7, 7)
    theme(  legend.position = "none",
      plot.title = element_text( size = 16, face = "bold" ),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12) )
  print(p)
  dev.off()
  return(results_with_genes) }

# do plots and get data
vol_D3  <- make_volcano(res_D3,  "D3",  t2g)
vol_D7  <- make_volcano(res_D7,  "D7",  t2g)
vol_D10  <- make_volcano(res_D10,  "D10",  t2g)
vol_D14 <- make_volcano(res_D14, "D14", t2g)
vol_D21 <- make_volcano(res_D21, "D21", t2g)
vol_D28 <- make_volcano(res_D28, "D28", t2g)
vol_D35 <- make_volcano(res_D35, "D35", t2g)

 d3_de = subset( res_D3, adj.P.Val<=0.05 & (logFC > 1.39| logFC < -1.39) )
 d7_de = subset( res_D7, adj.P.Val<=0.05 & (logFC > 1.39 | logFC < -1.39) )
d10_de = subset(res_D10, adj.P.Val<=0.05 & (logFC > 1.39 | logFC < -1.39) )
d14_de = subset(res_D14, adj.P.Val<=0.05 & (logFC > 1.39 | logFC < -1.39) )
d21_de = subset(res_D21, adj.P.Val<=0.05 & (logFC > 1.39 | logFC < -1.39) )
d28_de = subset(res_D28, adj.P.Val<=0.05 & (logFC > 1.39 | logFC < -1.39) )
d35_de = subset(res_D35, adj.P.Val<=0.05 & (logFC > 1.39 | logFC < -1.39) )

# get all unique DE get, min logFC
# d7_de, d10_de, d14_de, d28_de = duplicate
set1 <- rownames(rbind(d3_de, d7_de, d10_de, d14_de, d28_de, d35_de)) %in% rownames(d21_de)
str(set1) # 1103 transcripts
# only the top on in d3_de

all_de <- rbind(d3_de, d7_de, d10_de, d14_de, d28_de, d35_de)
all_de <- all_de[order(rownames(all_de)), ]
dim(all_de)
head(all_de) # 1103 unique genes

vol_all <- make_volcano(all_de, "all", t2g) # not good

# ---------------------------------------------------------
# Plot expression trajectories for all DE genes
#
# Assumes:
#   all_de = DE results table
#   rownames(all_de) = gene IDs
#   v = voom object
#   full = sample metadata
#
# Produces:
#   DE_gene_expression_panels.pdf
#   DE_gene_expression_panels.png
#
# Layout:
#   3 columns
#   ~157 rows (471 genes)
#   free y-axis per gene
#   coloured by day
# ---------------------------------------------------------

colnames(v$E) <- full$sample
viridis_8 <- hcl.colors(8,"viridis")
base_colors <- c( "D0"=viridis_8[1], "D3"=viridis_8[2],
  "D7"=viridis_8[3],  "D10"=viridis_8[4], "D14"=viridis_8[5],
  "D21"=viridis_8[6],  "D28"=viridis_8[7],  "D35"=viridis_8[8] )

gene_lookup <- t2g %>%  dplyr::select(target_id,ens_gene,ext_gene) %>%
  distinct()

gene_lookup$symbol <- ifelse(
  !is.na(gene_lookup$ext_gene) & gene_lookup$ext_gene!="",
  gene_lookup$ext_gene, gene_lookup$ens_gene )

expr_mat <- as.data.frame(v$E)
expr_mat$gene <- rownames(expr_mat)
expr_mat <- expr_mat %>%   filter(gene %in% rownames(all_de))

if(grepl("^ENSSSCT",rownames(all_de)[1])){
  expr_mat <- expr_mat %>%    left_join(
      gene_lookup %>%     dplyr::select(target_id,symbol),
      by=c("gene"="target_id")  ) }else{
  expr_mat <- expr_mat %>%    left_join(
      gene_lookup %>%        dplyr::select(ens_gene,symbol),
      by=c("gene"="ens_gene") ) }

expr_mat$plot_name <- ifelse(is.na(expr_mat$symbol) | expr_mat$symbol=="",
  expr_mat$gene,  paste0(expr_mat$symbol," | ",expr_mat$gene))

expr_long <- expr_mat %>%   pivot_longer(
    cols=-c(gene,symbol,plot_name),names_to="sample", values_to="expression"
  ) %>%  left_join( full[,c("sample","condition")],  by="sample"  )

expr_long$day <- factor( paste0("D",expr_long$condition),
  levels=c( "D0","D3","D7","D10","D14","D21","D28","D35")) 

expr_cluster <- expr_long %>%  group_by(gene,day) %>%
  summarise(    mean_expr=mean(expression,na.rm=TRUE),
    .groups="drop" ) %>%
  pivot_wider(    names_from=day, values_from=mean_expr )
gene_ids <- expr_cluster$gene
cluster_matrix <- expr_cluster %>%  dplyr::select(-gene) %>%
  as.matrix()
rownames(cluster_matrix) <- gene_ids

gene_cor <- cor(  t(cluster_matrix), use="pairwise.complete.obs" )
hc <- hclust(  as.dist(1-gene_cor),  method="average" )
gene_order <- gene_ids[hc$order]
plot_lookup <- expr_long %>% dplyr::select(gene,plot_name) %>%
  distinct()

plot_order <- plot_lookup$plot_name[match(gene_order,plot_lookup$gene)]
expr_long$plot_name <- factor(expr_long$plot_name,levels=plot_order )
str(expr_long) 
all_de <- read.csv("LIMMA/all.DE.csv", header=T, sep=",") # Read CSV file
str(all_de)
# get  expr_long from all_de based on target_id or ens_gene
# thus making a new smaller expr_long for plotting

expr_long3 <- expr_long %>%
  dplyr::filter(gene %in% all_de$target_id)

expr_long3$day <- factor(
  expr_long3$day,
  levels = c("D0","D3","D7","D10","D14","D21","D28","D35")
)

expr_long3$plot_name <- factor(
  expr_long3$plot_name,
  levels = plot_order
)

# use gene symbol as strip text but keep transcript-level faceting
lab_df <- expr_long3 %>%
  dplyr::select(plot_name, symbol) %>%
  dplyr::distinct()

labeller_gene <- function(x) {
  lab_df$symbol[match(x, lab_df$plot_name)]
}

p <- ggplot(expr_long3,
            aes(x = day, y = expression, colour = day)) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0),
    size = 3,
    alpha = 0.6
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    colour = "black",
    linewidth = 1.2,
    alpha = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    colour = "black",
    size = 1.1
  ) +
  scale_colour_manual(values = base_colors) +
  facet_wrap(
    ~plot_name,
    ncol = 15,
    scales = "free_y",
    labeller = labeller(plot_name = labeller_gene)
  ) +
  labs(x = NULL, y = "Normalised expression") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 5, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.spacing = unit(0.02, "lines"))
ggsave( "DE_gene_expression_panels.pdf", p,  width =30,
  height = 60,  dpi = 300, limitsize = F)
ggsave(  "DE_gene_expression_panels.png", p,  width =30,
  height = 60,  dpi = 300, limitsize = FALSE,  bg = "white")

#########

#pdf("ENSSSCT000000973541.pdf", height=6, width=8)
#plot_bootstrap(so2, "ENSSSCT000000973541.1", x_axis_angle = 90,
#               units="est_counts", color_by="condition")
#dev.off()

# Spline model # "Does expression change over time in any non-linear way?"

library(splines)
fitsp <- lmFit(v, model.matrix(  ~ ns(day, df=3),  data=full))
fitsp <- eBayes(fitsp, trend = T)
res <- topTableF(fitsp ,  number=Inf,  p.value = 1) # limma F-test
dim(res) # XX DE
write.csv(res, "spline.toptable.csv") # 29.5k transcripts 
head(res[order(res$adj.P.Val), ])
sig_spline <- subset(res, adj.P.Val < 0.05)  # 25,500 !!! transcripts 

# plot ignificance-versus-effect-size plot 
res$significant <- ifelse( res$adj.P.Val < 0.05,
    "Significant", "Not" )
spline1 <- ggplot( res, aes(  x = log10(F), y = -log10(adj.P.Val),
        colour = significant    )) +
    geom_point(alpha = 0.4) +    theme_bw()
pdf("spline.pdf", width=7, height=7)
print(spline1)
dev.off()

####### 5 Sleuth across genes ###########

so2 <- sleuth_prep(s2c, target_mapping=t2g, aggregation_column='ens_gene',
                  extra_bootstrap_summary=T, read_bootstrap_tpm=T)

# Next, we will smooth the tpm per sample using a parameter
# based on our model - so here we estimate parameters for
# response error measurement (full) model
# this is our alternative model with DE
so2 <- sleuth_fit(  so2, ~ ns(condition, df=3),  'full')

# get our null model r where the isoform levels are equal
so2 <- sleuth_fit(  so2, ~1, 'reduced' )
so2 <- sleuth_lrt(  so2,  'reduced', 'full')
# we compare our null and alternative models
models(so2) # check model details

# differential analysis using a likelihood ratio test (LRT)
sr <- sleuth_results(so2, 'reduced:full', 'lrt', show_all=F)
sr2 <- dplyr::filter(sr, qval <= 0.05)  
str(sr2)
write.csv(sr2, "sleuth_DE.csv") # 4793 DE genes
##########

library(readxl)
# list genes with multiple transcripts #

# all_de <- read_xlsx("TABLES/All_DE.xlsx") # Read Excel file
all_de <- read.csv("LIMMA/all.DE.csv", header=T, sep=",") # Read CSV file
str(all_de)
length(unique(all_de$ens_gene))
gene_counts <- all_de %>% count(ens_gene, name = "n_transcripts")
multi_transcript_genes <- gene_counts %>% dplyr::filter(n_transcripts > 1)
multi_transcript_de <- all_de %>% # Extract all rows for those genes
  semi_join(multi_transcript_genes, by = "ens_gene") %>% arrange(ens_gene)

# View(multi_transcript_de)
dim(multi_transcript_de) # 221+ genes
write.csv(  multi_transcript_de,  "all.DE.multiple_transcripts.csv",row.names = F )

write.csv( unique(multi_transcript_de$ens_gene), "all.DE.multiple_txn.csv",row.names = F )

### PRRSV data

library(readr) 

expr <- read_csv("MAPPING/gene_expression_log.csv") 
sample_info <- data.frame(  sample = c(
    "P22-4042","P22-4050","P22-4032",
    "P22-4044","P22-4045","P22-4040",
    "P22-4043","P22-4048","P22-4052",
    "P22-4038","P22-4051","P22-4039",
    "P22-4035","P22-4047","P22-4041",
    "P22-4037","P22-4036","P22-4046",
    "P22-4033","P22-4034","P22-4049"  ),
  day = c( rep(3,3),   rep(7,3),  rep(10,3),  rep(14,3),
    rep(21,3), rep(28,3),  rep(35,3) ) ) 
expr <- left_join( expr,  sample_info,  by = "sample" ) 

sample_means <- expr %>% group_by(sample, day) %>%
  summarise(  mean_expression = mean(expr,na.rm=T),.groups = "drop" )

# mean per timepoint
time_means <- sample_means %>% group_by(day) %>%
  summarise( mean_expression = mean(mean_expression),.groups = "drop"  )

# add D0 baseline
time_means <- bind_rows(  data.frame(day=0,mean_expression=0), time_means)

p <- ggplot() + geom_point( data = sample_means,
    aes( x = factor(day, levels = c(0,3,7,10,14,21,28,35)),
      y = mean_expression, colour = as.character(day)  ),
    size = 5, alpha = 0.8,  position = position_jitter(
      width = 0.05, height = 0 )) +
  geom_line( data = time_means,
    aes(  x = factor(day, levels = c(0,3,7,10,14,21,28,35)),
      y = mean_expression,
      group = 1 ),  colour = "black", alpha=0.6, linewidth = 1.5 ) +
  geom_point( data = time_means,
    aes( x = factor(day, levels = c(0,3,7,10,14,21,28,35)),
      y = mean_expression ),
    colour = "black", size =6, alpha=0.7 ) +
  scale_x_discrete(
    limits = c("0","3","7","10","14","21","28","35")) +
  scale_colour_manual(    values = base_colours ) +
  labs(   x = "Day post infection",
    y = "Log10-scaled mean expression" ) +
  theme_bw(base_size = 16) +
  theme(    legend.position = "none"  )
print(p)
View(sample_means)
ggsave( "Mean_gene_expression_by_timepoint.pdf", p,
  width =6,  height =5 )
ggsave(  "Mean_gene_expression_by_timepoint.png",
  p, width =6, height =5,  dpi = 300)

######

# model DE genes vs viral profile

# DE genes
all_de <- read.csv("LIMMA/all.DE.csv", stringsAsFactors = FALSE)

expr_long2 <- expr_long %>%  dplyr::filter(symbol %in% all_de$gene_label)

# Viral trajectory
viral_profile <- sample_means %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(
    viral_expr = mean(mean_expression, na.rm = TRUE),
    .groups = "drop"
  )

viral_profile$day <- paste0("D", viral_profile$day)

# Mean gene expression per day
pig_profile <- expr_long2 %>%
  dplyr::group_by(gene, plot_name, day) %>%
  dplyr::summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

# Correlation with viral trajectory
cor_table <- dplyr::bind_rows(
  lapply(unique(pig_profile$gene), function(g) {

    x <- merge(
      pig_profile[pig_profile$gene == g, ],
      viral_profile,
      by = "day"
    )

    if (nrow(x) < 3) return(NULL)

    data.frame(
      gene = g,
      viral_cor = cor(
        x$mean_expr,
        x$viral_expr,
        use = "pairwise.complete.obs"
      )
    )

  })
)

# Gene labels
plot_lookup <- expr_long2 %>%
  dplyr::select(gene, plot_name) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Rank by viral correlation
cor_table2 <- cor_table %>%
  dplyr::left_join(plot_lookup, by = "gene") %>%
  dplyr::arrange(dplyr::desc(viral_cor)) %>%
  dplyr::mutate(
    facet_order = dplyr::row_number(),
    facet_id = paste0(
      sprintf("%04d", facet_order),
      "__",
      plot_name,
      " | r=",
      sprintf("%.3f", viral_cor)
    )
  )

# Save correlation table
write.csv(
  cor_table2,
  "PigGenes_correlated_with_viral_trajectory.csv",
  row.names = FALSE
)

# Final plotting dataframe
expr_long2 <- expr_long2 %>%
  dplyr::left_join(
    cor_table2 %>%
      dplyr::select(gene, facet_id),
    by = "gene" )
expr_long2$day <- factor(
  as.character(expr_long2$day),
  levels = c("D0","D3","D7","D10","D14","D21","D28","D35") ) 
str(expr_long2)

expr_long2$day <- factor(
  as.character(expr_long2$day),
  levels = c("D0","D3","D7","D10","D14","D21","D28","D35")
)

p <- ggplot(
  expr_long2,
  aes(x = day, y = expression)
) +
  geom_point(
    aes(colour = day),
    position = position_jitter(width = 0.15),
    size = 0.5,
    alpha = 0.5
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    colour = "black",
    linewidth = 0.5
  ) +
  facet_wrap(
    ~ facet_id,
    ncol = 20,
    scales = "free_y" ) +
  scale_colour_manual(values = base_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 3),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(  "PigGenes_ranked_by_ViralCorrelation.pdf",
  p,width = 30, height = 80, dpi = 300, limitsize = F )
ggsave(  "PigGenes_ranked_by_ViralCorrelation.png",
  p,width = 30, height = 80, dpi = 300,  bg = "white", limitsize = F )
 
#### ppins #### 

library(dplyr)
library(STRINGdb)
library(igraph)
library(ggraph)

# DE genes
genes_for_string <- data.frame(
  gene = unique(all_de$gene_label)
)

# STRING pig
string_db <- STRINGdb$new(
  version = "12",
  species = 9823,
  score_threshold = 700
)

# map genes
mapped <- string_db$map(
  genes_for_string,
  "gene"
)

cat(
  "Input genes:", nrow(genes_for_string), "\n",
  "Mapped genes:", nrow(mapped), "\n"
)

# retrieve interactions
ppi <- string_db$get_interactions(mapped$STRING_id)

ids <- unique(mapped$STRING_id)

ppi <- ppi %>%
  dplyr::filter(
    from %in% ids,
    to %in% ids
  )

# lookup table
lookup <- mapped %>%
  dplyr::select(STRING_id, gene)

edges <- ppi %>%
  left_join(
    lookup,
    by = c("from" = "STRING_id")
  ) %>%
  rename(gene1 = gene) %>%
  left_join(
    lookup,
    by = c("to" = "STRING_id")
  ) %>%
  rename(gene2 = gene)

# graph
g <- graph_from_data_frame(
  edges[, c("gene1","gene2")],
  directed = FALSE )
# degree
V(g)$degree <- degree(g)

g2 <- subset(g, V(g)$degree > 4) # save tables
write.csv(
  edges,
  "STRING_edges.csv",
  row.names = FALSE  )

write.csv(  data.frame(    gene = V(g)$name,
    degree = V(g)$degree ),
  "STRING_nodes.csv",  row.names = FALSE)
# plot
p <- ggraph(g, layout = "fr") +
  geom_edge_link(
    colour = "grey80",
    alpha = 0.3) +
  geom_node_point(
    aes(size = degree),
    colour = "firebrick") +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3) +
  theme_void()

ggsave(
  "STRING_PPIN.pdf",
  p,
  width = 12,
  height = 10)

ggsave(
  "STRING_PPIN.png",
  p,
  width = 12,
  height = 10,
  dpi = 300)

cat(  "Nodes:", vcount(g), "\n",
  "Edges:", ecount(g), "\n")

# largest connected component
comp <- components(g)
g2 <- induced_subgraph( g,
  V(g)[comp$membership == which.max(comp$csize)] )
# keep top 50 hubs
deg <- degree(g2)
g3 <- induced_subgraph(
  g2,  names(sort(deg, decreasing = TRUE))[1:50])

#V(g2)$degree <- degree(g2)
p2 <- ggraph(g2, layout = "fr") +
  geom_edge_link( colour = "grey80",alpha = 2) +
  geom_node_point( aes(size = degree*1.2),
    colour = "firebrick", alpha=0.8) +
  geom_node_text( aes(label = name), repel = TRUE,
    size = 4, max.overlaps = Inf) +  theme_void()

ggsave(  "STRING_PPIN2.pdf",
  p2,  width = 16, height = 10)

ggsave(  "STRING_PPIN2.png",
  p2,  width = 16, height = 10, dpi = 300)

str(g2)
V(g2)$name
write.csv(  data.frame(    gene = V(g2)$name,
  degree = V(g2)$degree ),
  "STRING_nodes_largest_component.csv",  row.names = FALSE)

## GSEA

BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("org.Ss.eg.db")
library(org.Ss.eg.db)

library(clusterProfiler)
library(dplyr)
library(biomaRt)

# Combine all limma contrasts
all_res <- bind_rows(
  res_D3,
  res_D7,
  res_D10,
  res_D14,
  res_D21,
  res_D28,
  res_D35,
  .id = "comparison"
)

all_res$target_id <- rownames(all_res)

# Add gene IDs
all_res <- left_join(
  all_res,
  t2g[, c("target_id", "ens_gene")],
  by = "target_id"
)

# One score per gene across all timepoints
gene_rank <- all_res %>%
  filter(!is.na(ens_gene), ens_gene != "") %>%
  group_by(ens_gene) %>%
  summarise(
    rank = sqrt(sum(t^2)),
    .groups = "drop"
  )

# Ensembl -> Entrez
entrez_map <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "entrezgene_id"
  ),
  mart = mart
)

gene_rank <- gene_rank %>%
  inner_join(
    entrez_map,
    by = c("ens_gene" = "ensembl_gene_id")
  ) %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct(entrezgene_id, .keep_all = TRUE)

# Ranked vector for GSEA
ranks <- gene_rank$rank
names(ranks) <- gene_rank$entrezgene_id
ranks <- sort(ranks, decreasing = TRUE)

# KEGG GSEA
gsea <- gseKEGG(
  geneList = ranks,
  organism = "ssc",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

write.csv(
  as.data.frame(gsea),
  "KEGG_GSEA_All_Timepoints.csv",
  row.names = FALSE
)

pdf("KEGG_GSEA_All_Timepoints.pdf", width = 8, height = 6)
print(dotplot(gsea, showCategory = 20))
dev.off()

head(as.data.frame(gsea))