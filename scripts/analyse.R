library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(limma)
library(splines)

# -------------------------------
# Settings
# -------------------------------
xlims <- c(0, 15010)
xbreaks <- seq(0, 16000, 1000)
input_folder <- "DEPTH/"

# -------------------------------
# CDS annotation
# -------------------------------
cds <- data.frame(
  start = c(105, 7277, 11679, 11684, 12287, 12829, 13377, 13382, 13970, 14481),
  end   = c(7295,11668,12428,11896,13084,13380,13982,13513,14491,14867),
  name  = c("1a replicase","1b replicase / RdRp","GP2 envelope",
            "GP2b envelope","GP3 envelope","GP4 envelope",
            "GP5 envelope","ORF5a","membrane envelope","nucleocapsid")
)
cds$mid <- (cds$start + cds$end)/2

# -------------------------------
# Sample mapping
# -------------------------------
sample_info <- data.frame(
  sample = c("P22-4042","P22-4050","P22-4032",
             "P22-4044","P22-4045","P22-4040",
             "P22-4043","P22-4048","P22-4052",
             "P22-4038","P22-4051","P22-4039",
             "P22-4035","P22-4047","P22-4041",
             "P22-4037","P22-4036","P22-4046",
             "P22-4033","P22-4034","P22-4049"),
  condition = c(rep(3,3), rep(7,3), rep(10,3),
                rep(14,3), rep(21,3), rep(28,3), rep(35,3))
)

# -------------------------------
# Sliding window
# -------------------------------
compute_sliding <- function(df, window_size=1000, step_size=100) {
  if(nrow(df) < window_size) return(NULL)
  starts <- seq(1, max(df$pos) - window_size, by=step_size)

  means <- sapply(starts, function(s) {
    subset <- df[df$pos >= s & df$pos < (s + window_size), "depth"]
    if(length(subset) > 0) mean(subset) else NA  })

  out <- data.frame(pos = starts + window_size/2,
                    depth = means)
  out <- out[complete.cases(out), ]
  if(nrow(out)==0) return(NULL)

  return(out)
}

files <- list.files(input_folder, pattern="\\.txt$", full.names=TRUE)

all_data <- lapply(files, function(file) {
  df <- read.table(file, sep="\t", header=FALSE,
                   col.names=c("chrom","pos","depth"))
  df <- df[complete.cases(df), ]
  sw <- compute_sliding(df)
  if(is.null(sw)) return(NULL)
  sample_name <- gsub("\\.PV173709\\.fasta\\.depth\\.txt$", "",
                      basename(file))
  sample_name <- gsub("_","-",sample_name)

  sw$sample <- sample_name
  return(sw)
}) %>% bind_rows()

all_data <- merge(all_data, sample_info, by="sample")
all_data <- all_data[all_data$depth > 0, ]
all_data$log_depth <- log10(all_data$depth)

gene_expr <- lapply(1:nrow(cds), function(i) {
  region <- cds[i,]
  sub <- all_data %>%
    filter(pos >= region$start & pos <= region$end)

  sub %>%    group_by(sample) %>%
    summarise(expr = mean(log_depth), .groups="drop") %>%
    mutate(gene = region$name) }) %>% bind_rows()

write.csv(gene_expr, "gene_expression_log.csv", row.names=FALSE)

expr_mat <- gene_expr %>% pivot_wider(names_from=sample, values_from=expr)

expr_mat <- as.data.frame(expr_mat)
rownames(expr_mat) <- expr_mat$gene
expr_mat$gene <- NULL
expr_mat_norm <- normalizeBetweenArrays(expr_mat)
write.csv(expr_mat_norm, "limma_input_matrix.csv")

# -------------------------------
# Spline model
# -------------------------------
sample_table <- sample_info
sample_table <- sample_table[
  match(colnames(expr_mat_norm), sample_table$sample), ]

design <- model.matrix(~ ns(condition, df=3), data=sample_table)

fit <- lmFit(expr_mat_norm, design)
fit <- eBayes(fit)

res <- topTable(fit, number=Inf)
write.csv(res, "limma_results_spline.csv")

anova_res <- topTableF(fit, number=Inf)
write.csv(anova_res, "limma_spline_ANOVA.csv")

plot_condition <- function(df, cond, yvar, ylabel) {
  sub <- df[df$condition == cond, ]
  if(nrow(sub)==0) return(NULL)
  ymax <- max(sub[[yvar]], na.rm=TRUE)
  label_y <- rep(c(ymax*1.1, ymax*0.95), length.out=nrow(cds))

  ggplot(sub, aes(x=pos, y=.data[[yvar]], colour=sample)) +
    geom_rect(data=cds,
              aes(xmin=start, xmax=end,
	      ymin=ymax*0.9, ymax=ymax*1.15),
              inherit.aes=FALSE,  fill="grey70", alpha=0.1) +
    geom_line(linewidth=1.2, alpha=0.5) +
    geom_label(data=cds,  aes(x=mid, y=label_y, label=name),
        inherit.aes=F, size=2, fill="black", colour="white") +
    scale_x_continuous(limits=xlims, breaks=xbreaks,
                       labels=function(x) x/1000) +
    expand_limits(y=ymax*1.2) +
    labs(title=paste0("D",cond),
         x="Genomic position (Kb)",
         y=ylabel,         colour="Sample") +
    theme_minimal(base_size=9) }

conds <- sort(unique(sample_info$condition))

# RAW
plots_raw <- lapply(conds, function(c)
  plot_condition(all_data, c, "log_depth", "Log10-scaled depth"))
panel_raw <- plot_grid(plotlist=plots_raw, ncol=1)
plot_df <- gene_expr %>% left_join(sample_info, by="sample")
top_genes <- rownames(anova_res)[1:6]

p_spline <- ggplot(
  plot_df %>% filter(gene %in% top_genes),
  aes(x=condition, y=expr)) +
  geom_point(aes(color=sample)) +
  geom_smooth(method="lm",
              formula = y ~ ns(x,3),
              se=FALSE,   colour="black",
              linewidth=1) +
  facet_wrap(~gene, scales="free_y") +
  theme_minimal() +
  labs(x="Day", y="Expression (log depth)")

#ggsave("spline_gene_trajectories.png",       p_spline,
#       width=10, height=6, dpi=400)

# -------------------------------
# Combine output panels
# -------------------------------
final_plot <- panel_raw

ggsave("depth_plot.png",
       final_plot,
       width=6, height=12,
       dpi=400, bg="white")
