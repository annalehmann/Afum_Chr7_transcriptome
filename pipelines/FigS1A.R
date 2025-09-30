library(tidyverse)
library(ggplot2)
library(patchwork)

# DNA seq coverage from each RNA-seq condition (Fig S1A)
setwd("--/bedfiles")

bednames <- list.files(pattern = ".bed")
samples <- sub(".regions.bed", "", bednames)
bedfiles <- list()

# chromosome start/end coords from CEA10 ref assembly
chr_position <- c(0, 4910137, 9697127, 13967564, 17787680, 21587582, 25488251, 27373204)
chr_name <- c("CP097563.1", "CP097564.1", "CP097565.1", "CP097566.1", "CP097567.1", "CP097568.1", "CP097569.1", "CP097570.1")

# read in files
for (f in seq_along(bednames)) {
  d <- read.delim(bednames[f], 
                  sep = "\t", 
                  header = FALSE, 
                  col.names = c("Chromosome", "Coordinate", "end", "average_coverage"))
  global_coords <- numeric(nrow(d)) 
  for (i in chr_name) {
    global_coord_index <- d$Chromosome == i
    global_coords[global_coord_index] <- d$Coordinate[global_coord_index] + chr_position[match(i, chr_name)]
  }
  
  d$global_coord <- global_coords
  
  global_avg <- mean(d$average_coverage)
  
  d$Normalized_Coverage <- d$average_coverage / global_avg
  bedfiles[[f]] <- d
}

names(bedfiles) <- samples

# last coord of each chr
breaks <- c(0, 4910137, 9697127, 13967564, 17787680, 21587582, 25488251, 27373204, 29322347)

colors <- c()
# no FK506
for (f in names(bedfiles)) {
  # assign colors to dif samples
  # euploid no drug
  if (grepl("Eu.*AMM", f)) {
    colors[f] <- "#a7adba"  
    # euploid +FK506
    } else if (grepl("Eu.*FK", f)) {
      colors[f] <- "#444545"  
    # aneuploid no drug
    } else if (grepl("An.*AMM", f)) {
      colors[f] <- "#b4a8d2"  
    # aneuploid +FK506
    } else if (grepl("An.*FK", f)) {
      colors[f] <- "#6850a1" 
    }
  # make plots
  p <- ggplot(bedfiles[[f]], aes(x = global_coord, y = Normalized_Coverage)) + 
    geom_point(aes(fill = Chromosome), color = colors[f], size = 0.001) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", 
                                      fill=NA, 
                                      linewidth = 0.25),
          axis.ticks.y = element_line(color = "black",
                                      linewidth = 0.25)) +
    geom_vline(xintercept = breaks,
               linetype = "dashed",
               linewidth = 0.25) +
    xlim(0, 29322347) +
    ylim(0,2.5)
  assign(f, p)
}

# make combined plots
euploid <- Eu01_AMM /
  Eu02_AMM /
  Eu03_AMM /
  Eu01_FK /
  Eu02_FK /
  Eu03_FK

ggsave(file="euploid_cov.svg", plot=euploid, width=6, height=4)

aneuploid <- An01_AMM /
  An02_AMM /
  An03_AMM /
  An01_FK /
  An02_FK /
  An03_FK

ggsave(file="aneuploid_cov.svg", plot=aneuploid, width=6, height=4)
