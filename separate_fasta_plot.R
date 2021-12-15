library(tidyverse)
library(gridExtra)

#Read and merge data
args <- commandArgs(trailingOnly = TRUE)
swell <- read.table(args[1], head=TRUE, sep='\t', na.strings=c("NaN", "-"))
metadata <- read.table(args[2], head=TRUE, sep='\t', fill=TRUE)
start_date <- args[3]
end_date <- args[4]
df <- merge(x=swell, y=metadata, by.x="header", by.y="fasta_header")

# Calculate ISO week for each resulting entry
df$sequencing_submission_week = as.integer(strftime(df$sequencing_submission_date, format="%V"))

# Point opacity parameter
alpha_val = 0.4

# File name prefix
if ((start_date != "None") && (end_date != "None")) {
    prefix = paste("fasta", start_date, end_date, sep="_")
} else if (start_date != "None") {
    prefix = paste("fasta", start_date, "present", sep="_")
} else if (end_date != "None") {
    prefix = paste("fasta", "beginning", end_date, sep="_")
} else {
    prefix = paste("fasta", "all", sep="_")
}

# Generating PDF of graphs similar to/the same as ones made by Nick Loman
pdf(file = paste(prefix, "_summary.pdf", sep=""))
# Percentage acgt scatterplot, sequencing org on x axis
print(ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_point(alpha = alpha_val) + 
    theme_bw())
# Percentage acgt scatterplot, faceted by sequencing org
print(ggplot(data = df, mapping = aes(x = pc_acgt, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank()))
# Percentage acgt histogram, faceted by sequencing org
print(ggplot(data = df, mapping = aes(x = pc_acgt)) +
    geom_histogram() + 
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw())
# Average pc_acgt per ISO week, coloured by sequencing org
print(ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_acgt, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "average pc_acgt"))
# Percentage ambiguous scatterplot, sequencing org on x axis
print(ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_point(alpha = alpha_val) + 
    theme_bw())
# Percentage ambiguous scatterplot, faceted by sequencing org
print(ggplot(data = df, mapping = aes(x = pc_ambiguous, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank()))
# Percentage ambiguous histogram, faceted by sequencing org
print(ggplot(data = df, mapping = aes(x = pc_ambiguous)) +
    geom_histogram() +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw())
# Average pc_ambiguous per ISO week, coloured by sequencing org
print(ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "average pc_ambiguous"))
dev.off()

# PDFs of similar stats, but split by sequencing_org_code and run_name
seq_org_codes <- unique(df$sequencing_org_code)
# Loop through sequencing_org_codes to create separate PDFs for each
for (code in seq_org_codes) {
    pdf(file = paste(prefix, "_", code, ".pdf", sep=""))
    code_df <- filter(df, sequencing_org_code == code)
    runs <- sort(unique(code_df$run_name), decreasing = FALSE)
    for (run in runs) {
        run_df <- filter(code_df, run_name == run)
        # Percentage acgt scatterplot
        p1 = ggplot(data = run_df, mapping = aes(x = pc_acgt, y = 0)) +
            geom_jitter(width = 0, alpha = alpha_val) +
            scale_x_continuous(limits = c(0, 100)) +
            scale_y_continuous(breaks = NULL) +
            theme_bw() +
            theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())
        # Percentage acgt histogram
        p2 = ggplot(data = run_df, mapping = aes(x = pc_acgt)) +
            geom_histogram() + 
            scale_x_continuous(limits = c(0, 100)) +
            theme_bw()
        # Percentage ambiguous scatterplot
        p3 = ggplot(data = run_df, mapping = aes(x = pc_ambiguous, y = 0)) +
            geom_jitter(width = 0, alpha = alpha_val) +
            scale_x_continuous(limits = c(0, 100)) +
            scale_y_continuous(breaks = NULL) +
            theme_bw() +
            theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())
        # Percentage ambiguous histogram
        p4 = ggplot(data = run_df, mapping = aes(x = pc_ambiguous)) +
            geom_histogram() +
            scale_x_continuous(limits = c(0, 100)) +
            theme_bw()
        # Organise into grid
        print(grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2, top = run))
    }
    dev.off()
}