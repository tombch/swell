library(tidyverse)
swell <- read.table(file("stdin"), head=TRUE, sep='\t', na.strings=c("NaN", "-"))

positive_coverage <- c(">=1", ">=5", ">=10", ">=20", ">=50", ">=100", ">=200")
percentage_of_tiles <- c(swell$pc_pos_cov_gte1, swell$pc_pos_cov_gte5, swell$pc_pos_cov_gte10, swell$pc_pos_cov_gte20, swell$pc_pos_cov_gte50, swell$pc_pos_cov_gte100, swell$pc_pos_cov_gte200)
pos_cov_df <- data.frame(positive_coverage , percentage_of_tiles)
pos_cov_df$percentage_of_tiles <- as.numeric(as.character(pos_cov_df$percentage_of_tiles))
# Positive coverage graph
ggplot(data = pos_cov_df, mapping = aes(x = positive_coverage, y = percentage_of_tiles)) +
    geom_col() + 
    scale_x_discrete(limits=pos_cov_df$positive_coverage) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_bw()

median_coverage <- c(">=1", ">=5", ">=10", ">=20", ">=50", ">=100", ">=200")
percentage_of_tiles <- c(swell$pc_tiles_medcov_gte1, swell$pc_tiles_medcov_gte5, swell$pc_tiles_medcov_gte10, swell$pc_tiles_medcov_gte20, swell$pc_tiles_medcov_gte50, swell$pc_tiles_medcov_gte100, swell$pc_tiles_medcov_gte200)
med_cov_df <- data.frame(median_coverage , percentage_of_tiles)
med_cov_df$percentage_of_tiles <- as.numeric(as.character(med_cov_df$percentage_of_tiles))
# Median coverage graph
ggplot(data = med_cov_df, mapping = aes(x = median_coverage, y = percentage_of_tiles)) +
    geom_col() + 
    scale_x_discrete(limits=med_cov_df$median_coverage) +
    scale_y_continuous(limits = c(0, 100)) +
    theme_bw()

depth <- unlist(strsplit(swell$tile_vector, split=','))
tile_df <- data.frame(depth)
tile_df$depth <- as.numeric(as.character(tile_df$depth))
# Add position column to data frame
tile_df$position <- seq.int(nrow(tile_df))
# Tile depth graph
ggplot(data = tile_df, mapping = aes(x = position, y = depth)) +
    geom_col() + 
    theme_bw()
# Tile depth histogram
ggplot(data = tile_df, mapping = aes(x = depth)) +
    geom_histogram() + 
    theme_bw()