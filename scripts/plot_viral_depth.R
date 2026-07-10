# -------------------------------
# Libraries
# -------------------------------
library(ggplot2)
library(cowplot)
library(dplyr)

# -------------------------------
# Settings
# -------------------------------
xlims <- c(0, 15010)
xbreaks <- seq(0, 16000, 1000)
input_folder <- "DEPTH/"

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
                rep(14,3), rep(21,3), rep(28,3), rep(35,3)),
  stringsAsFactors = FALSE
)

# -------------------------------
# Sliding window
# -------------------------------
compute_sliding <- function(df, window_size=1000, step_size=100) {
  
  starts <- seq(1, max(df$pos) - window_size, by=step_size)
  
  means <- sapply(starts, function(s) {
    subset <- df[df$pos >= s & df$pos < (s + window_size), "depth"]
    if (length(subset) > 0) mean(subset) else NA
  })
  
  data.frame(
    pos = starts + window_size/2,
    depth = means
  ) %>% na.omit()
}

# -------------------------------
# Read depth files
# -------------------------------
files <- list.files(input_folder, pattern="\\.txt$", full.names=TRUE)

all_data <- lapply(files, function(file) {
  
  df <- read.table(file, sep="\t", header=FALSE,
                   col.names=c("chrom","pos","depth"))
  
  df <- df[complete.cases(df), ]
  
  sw <- compute_sliding(df)
  
  sample_name <- gsub("\\.PV173709\\.fasta\\.depth\\.txt$", "",
                      basename(file))
  sample_name <- gsub("_","-",sample_name)
  
  sw$sample <- sample_name
  return(sw)
}) %>% bind_rows()

# merge condition
all_data <- merge(all_data, sample_info, by="sample")

# -------------------------------
# ORIGINAL SCALE
# -------------------------------
all_data <- all_data[all_data$depth > 0, ]
all_data$log_depth <- log10(all_data$depth)

# -------------------------------
# NORMALISED (Z-score per sample)
# -------------------------------
all_data <- all_data %>%
  group_by(sample) %>%
  mutate(z_depth = (log_depth - mean(log_depth)) / sd(log_depth)) %>%
  ungroup()

# -------------------------------
# Plot function
# -------------------------------
plot_condition <- function(df, cond, yvar, ylabel) {
  
  sub <- df[df$condition == cond, ]
  
  ggplot(sub, aes_string(x="pos", y=yvar, colour="sample")) +
    geom_line(
      linewidth = 1.2,   # ✅ thicker lines
      alpha = 0.5        # ✅ transparency
    ) +
    scale_x_continuous(
      limits=xlims,
      breaks=xbreaks
    ) +
    labs(
      title=paste("Timepoint", cond),
      x="Position",
      y=ylabel,
      colour="Sample"
    ) +
    theme_minimal(base_size=9) +
    theme(
      legend.position="right"
    )
}

conditions <- sort(unique(sample_info$condition))

# -------------------------------
# RAW panels
# -------------------------------
plots_raw <- lapply(conditions, function(cond) {
  plot_condition(all_data, cond,
                 "log_depth",
                 "Log10-scaled depth")
})

panel_raw <- plot_grid(plotlist=plots_raw, ncol=1)

# -------------------------------
# NORMALISED panels
# -------------------------------
plots_norm <- lapply(conditions, function(cond) {
  plot_condition(all_data, cond,
                 "z_depth",
                 "Z-scored depth")
})

panel_norm <- plot_grid(plotlist=plots_norm, ncol=1)

# -------------------------------
# Combine panels
# -------------------------------
final_plot <- plot_grid(  panel_raw,
  panel_norm,  ncol=2,
  labels=c("Raw depth", "Normalised (Z-score)") )

ggsave("depth_comparison_panel.png", final_plot, width=12,
       height=15,      dpi=400,       bg="white")
ggsave("depth_comparison_panel.pdf", final_plot, width=12,
       height=15,      dpi=400,       bg="white")