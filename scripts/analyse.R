library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

# -------------------------------
# Settings (1 to 14,999 bp)
# -------------------------------
xlims <- c(1, 14999)
xbreaks <- seq(0, 15000, 1000)
input_folder <- "DEPTH/"

# -------------------------------
# CDS Annotation & Alternating Tracks
# -------------------------------
cds <- data.frame(
  start = c(105, 1274,  7277, 11679, 11684, 12287, 12829, 13377, 13382, 13970, 14481, 14974),
  end   = c(7295, 4648, 11668, 12428, 11896, 13084, 13380, 13982, 13513, 14491, 14867, 14975),
  name  = c("1a replicase", "nsp2", "RdRp", "GP2",
            "GP2b", "GP3", "GP4",            "GP5", "", "M", "N", "")
)

# Remove unnamed annotations and calculate midpoints
cds <- cds[cds$name != "", ]
cds$mid <- (cds$start + cds$end) / 2

# Assign alternating tracks (1 & 2) so overlapping genes stack cleanly
cds$track <- rep(c(1, 2), length.out = nrow(cds))

# -------------------------------
# Sample Mapping
# -------------------------------
sample_info <- data.frame(
  sample = c("P22-4042", "P22-4050", "P22-4032",
             "P22-4044", "P22-4045", "P22-4040",
             "P22-4043", "P22-4048", "P22-4052",
             "P22-4038", "P22-4051", "P22-4039",
             "P22-4035", "P22-4047", "P22-4041",
             "P22-4037", "P22-4036", "P22-4046",
             "P22-4033", "P22-4034", "P22-4049"),
  condition = c(rep(3, 3), rep(7, 3), rep(10, 3),
                rep(14, 3), rep(21, 3), rep(28, 3), rep(35, 3))
)

# -------------------------------
# Sliding Window (Spans 1 to 14,999 bp)
# -------------------------------
compute_sliding <- function(df, window_size=500, step_size=100, min_pos=1, max_pos=14999) {
  if(nrow(df) < 10) return(NULL)
  
  starts <- seq(min_pos, max_pos, by=step_size)
  
  sw <- lapply(starts, function(s) {
    w_start <- max(min_pos, s - window_size/2)
    w_end   <- min(max_pos, s + window_size/2)
    
    subset <- df[df$pos >= w_start & df$pos <= w_end, "depth"]
    
    if(length(subset) > 0) {
      data.frame(pos = s, depth = mean(subset, na.rm=TRUE))
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  return(sw)
}

# -------------------------------
# Data Import & Processing
# -------------------------------
files <- list.files(input_folder, pattern="\\.txt$", full.names=TRUE)

all_data <- lapply(files, function(file) {
  df <- read.table(file, sep="\t", header=FALSE,
                   col.names=c("chrom", "pos", "depth"))
  df <- df[complete.cases(df), ]
  
  # Pad missing positions from 1 to 14,999 bp
  full_pos <- data.frame(pos = 1:14999)
  df <- merge(full_pos, df[, c("pos", "depth")], by="pos", all.x=TRUE)
  df$depth[is.na(df$depth)] <- 0
  
  sw <- compute_sliding(df, window_size=500, step_size=100, min_pos=1, max_pos=14999)
  if(is.null(sw)) return(NULL)
  
  sample_name <- gsub("\\.PV173709\\.fasta\\.depth\\.txt$", "", basename(file))
  sample_name <- gsub("_", "-", sample_name)
  sw$sample <- sample_name
  return(sw)
}) %>% bind_rows()

all_data <- merge(all_data, sample_info, by="sample")
all_data <- all_data[all_data$depth > 0, ]
all_data$log_depth <- log10(all_data$depth+1)

# -------------------------------
# Global Y-Axis Parameters
# -------------------------------
global_ymax <- max(all_data$log_depth, na.rm=TRUE)
# Calculate uniform round breaks for the data portion of the axis
ybreaks <- seq(0, floor(global_ymax), by=1)

# -------------------------------
# Plotting Function (Uniform Y Scale)
# -------------------------------
plot_condition <- function(df, cond, yvar, ylabel, ymax_data) {
  sub <- df[df$condition == cond, ]
  if(nrow(sub) == 0) return(NULL)

  # Dynamic height placement for CDS tracks based on uniform global_ymax
  cds_plot <- cds
  cds_plot$ymin <- ifelse(cds_plot$track == 1, ymax_data * 1.05, ymax_data * 1.20)
  cds_plot$ymax <- ifelse(cds_plot$track == 1, ymax_data * 1.15, ymax_data * 1.30)
  cds_plot$mid_y <- (cds_plot$ymin + cds_plot$ymax) / 2

  ggplot(sub, aes(x=pos, y=.data[[yvar]], colour=sample)) +
    # Staggered grey boxes for gene bounds
    geom_rect(data=cds_plot, aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax),
              inherit.aes=FALSE, fill="grey80", color="grey40", alpha=0.8) +
    # Sliding window depth line plot
    geom_line(linewidth=1.2, alpha=0.5) +
    # Gene names centered inside grey boxes
    geom_text(data=cds_plot, aes(x=mid, y=mid_y, label=name),
              inherit.aes=FALSE, size=2, fontface="bold", colour="black") +
    scale_x_continuous(limits=xlims, breaks=xbreaks, labels=function(x) x/1000) +
    # Fixed y-limits and uniform tick breaks across every panel
    scale_y_continuous(limits=c(0, ymax_data * 1.35), breaks=ybreaks) +
    labs(title=paste0("D", cond),
         x="Genomic position (Kb)",
         y=ylabel, colour="Sample") +
    theme_minimal(base_size=9) +
    theme(panel.grid.minor = element_blank())
}

# -------------------------------
# Generate Output
# -------------------------------
conds <- sort(unique(sample_info$condition))

# Generate plots with uniform ymax
plots_raw <- lapply(conds, function(c)
  plot_condition(all_data, c, "log_depth", "Log10-scaled depth", ymax_data = global_ymax))

panel_raw <- plot_grid(plotlist=plots_raw, ncol=1)

# Save
final_plot <- panel_raw
ggsave("depth_plot.png", final_plot, width=6, height=12,
       dpi=400, bg="white")