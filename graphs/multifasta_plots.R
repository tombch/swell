library(tidyverse)
library(gridExtra)

# Prevent Rplots.pdf generation
pdf(NULL)

#Read and merge data
args <- commandArgs(trailingOnly = TRUE)
swell <- read.table(args[1], head=TRUE, sep='\t', na.strings=c("NaN", "-"))
metadata <- read.table(args[2], head=TRUE, sep='\t', fill=TRUE)
start_date <- args[3]
end_date <- args[4]
graph_runs <- args[5]
df <- merge(x=swell, y=metadata, by.x="header", by.y="fasta_header")

# Calculate ISO week for each resulting entry
df$sequencing_submission_week = strftime(df$sequencing_submission_date, format="%V")

# Change LOND_BART to BART for readability when sequencing_org_code is on x axis
 df$sequencing_org_code[df$sequencing_org_code=="LOND_BART"] <- "BART"

# Point opacity parameter
alpha_val = 0.4

# File name prefix
if ((start_date != "None") && (end_date != "None")) {
    prefix = paste(start_date, end_date, sep="_")
} else if (start_date != "None") {
    prefix = paste(start_date, "present", sep="_")
} else if (end_date != "None") {
    prefix = paste("beginning", end_date, sep="_")
} else {
    prefix = paste("all", sep="_")
}

# PNG width and height
w = 16
h = 6

# Average pc_acgt per ISO week, coloured by sequencing org
acgt_iso = ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_acgt, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "avg pc_acgt")
# Average pc_ambiguous per ISO week, coloured by sequencing org
ambig_iso = ggplot(data = df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous, colour = sequencing_org_code)) + 
    geom_point(stat = "summary", fun = "mean") + 
    geom_line(stat = "summary", fun = "mean") + 
    labs(y = "avg pc_ambiguous")
# Organise into grid
iso = grid.arrange(acgt_iso, ambig_iso, ncol=2, nrow=1)
ggsave(paste(prefix, "iso_week.png", sep="_"), iso, width=w, height=h)


# Percentage acgt scatterplot, sequencing org on x axis
acgt_scatter_seq = ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_point(alpha = alpha_val) + 
    theme_bw()
# Percentage ambiguous scatterplot, sequencing org on x axis
ambig_scatter_seq = ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_point(alpha = alpha_val) + 
    theme_bw()
# Organise into grid
scatter_seq = grid.arrange(acgt_scatter_seq, ambig_scatter_seq, ncol=2, nrow=1)
ggsave(paste(prefix, "scatter_seq.png", sep="_"), scatter_seq, width=w, height=h)


# Percentage acgt violin, sequencing org on x axis
acgt_violin = ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_violin() + 
    theme_bw()
# Percentage ambiguous violin, sequencing org on x axis
ambig_violin = ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_violin() + 
    theme_bw()
# Organise into grid
violin = grid.arrange(acgt_violin, ambig_violin, ncol=2, nrow=1)
ggsave(paste(prefix, "violin.png", sep="_"), violin, width=w, height=h)


# Percentage acgt boxplot, sequencing org on x axis
acgt_boxplot = ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_acgt)) +
    geom_boxplot() + 
    theme_bw()
# Percentage ambiguous boxplot, sequencing org on x axis
ambig_boxplot = ggplot(data = df, mapping = aes(x = sequencing_org_code, y = pc_ambiguous)) +
    geom_boxplot() + 
    theme_bw()
# Organise into grid
boxplot = grid.arrange(acgt_boxplot, ambig_boxplot, ncol=2, nrow=1)
ggsave(paste(prefix, "boxplot.png", sep="_"), boxplot, width=w, height=h)


# Percentage acgt scatterplot, faceted by sequencing org
acgt_scatter = ggplot(data = df, mapping = aes(x = pc_acgt, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())
# Percentage ambiguous scatterplot, faceted by sequencing org
ambig_scatter = ggplot(data = df, mapping = aes(x = pc_ambiguous, y = 0)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    scale_y_continuous(breaks = NULL) +
    theme_bw() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())
# Organise into grid
scatter = grid.arrange(acgt_scatter, ambig_scatter, ncol=2, nrow=1)
ggsave(paste(prefix, "scatter.png", sep="_"), scatter, width=w, height=h)


# Percentage acgt scatterplot, week on y axis, faceted by sequencing org
acgt_scatter_week = ggplot(data = df, mapping = aes(x = pc_acgt, y = sequencing_submission_week)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()
# Percentage ambiguous scatterplot, week on y axis, faceted by sequencing org
ambig_scatter_week = ggplot(data = df, mapping = aes(x = pc_ambiguous, y = sequencing_submission_week)) +
    geom_jitter(width = 0, alpha = alpha_val) +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()
scatter_week = grid.arrange(acgt_scatter_week, ambig_scatter_week, ncol=2, nrow=1)
ggsave(paste(prefix, "scatter_week.png", sep="_"), scatter_week, width=w, height=h)


# Percentage acgt histogram, faceted by sequencing org
acgt_hist = ggplot(data = df, mapping = aes(x = pc_acgt)) +
    geom_histogram() + 
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()
# Percentage ambiguous histogram, faceted by sequencing org
ambig_hist = ggplot(data = df, mapping = aes(x = pc_ambiguous)) +
    geom_histogram() +
    facet_wrap(. ~ sequencing_org_code, ncol=3, scale="free_y") + 
    theme_bw()
hist = grid.arrange(acgt_hist, ambig_hist, ncol=2, nrow=1)
ggsave(paste(prefix, "hist.png", sep="_"), hist, width=w, height=h)


seq_org_codes <- unique(df$sequencing_org_code)
for (code in seq_org_codes) {
    code_df <- filter(df, sequencing_org_code == code)
    runs <- sort(unique(code_df$run_name), decreasing = FALSE)

    # Average pc_acgt per ISO week, coloured by sequencing org
    acgt_iso_runs = ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_acgt, group = run_name)) + 
        geom_point(stat = "summary", fun = "mean") + 
        geom_line(stat = "summary", fun = "mean") + 
        labs(y = "average pc_acgt")
    # Average pc_ambiguous per ISO week, coloured by run_name
    ambig_iso_runs = ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous, group = run_name)) + 
        geom_point(stat = "summary", fun = "mean") + 
        geom_line(stat = "summary", fun = "mean") + 
        labs(y = "average pc_ambiguous")
    # Organise into grid
    iso_runs = grid.arrange(acgt_iso_runs, ambig_iso_runs, ncol=2, nrow=1, top = code)
    ggsave(paste(prefix, code, "runs_iso_week.png", sep="_"), iso_runs, width=w, height=h)


    # Percentage acgt violin, ISO week on x axis
    acgt_violin_iso = ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_acgt)) +
        geom_violin() + 
        theme_bw()
    # Percentage ambiguous violin, ISO week on x axis
    ambig_violin_iso = ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous)) +
        geom_violin() + 
        theme_bw()
    # Organise into grid
    violin_iso = grid.arrange(acgt_violin_iso, ambig_violin_iso, ncol=2, nrow=1)
    ggsave(paste(prefix, code, "violin_iso_week.png", sep="_"), violin_iso, width=w, height=h)


    # Percentage acgt boxplot, ISO week on x axis
    acgt_boxplot_iso = ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_acgt)) +
        geom_boxplot() + 
        theme_bw()
    # Percentage ambiguous boxplot, ISO week on x axis
    ambig_boxplot_iso = ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous)) +
        geom_boxplot() + 
        theme_bw()
    # Organise into grid
    boxplot_iso = grid.arrange(acgt_boxplot_iso, ambig_boxplot_iso, ncol=2, nrow=1)
    ggsave(paste(prefix, code, "boxplot_iso_week.png", sep="_"), boxplot_iso, width=w, height=h)
}


if (graph_runs == "True") {
    # PDFs of similar stats, but split by sequencing_org_code and run_name
    seq_org_codes <- unique(df$sequencing_org_code)
    # Loop through sequencing_org_codes to create separate PDFs for each
    for (code in seq_org_codes) {
        pdf(file = paste(prefix, "_", code, ".pdf", sep=""))
        code_df <- filter(df, sequencing_org_code == code)
        runs <- sort(unique(code_df$run_name), decreasing = FALSE)

        # Average pc_acgt per ISO week, coloured by sequencing org
        print(ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_acgt, colour = run_name)) + 
            geom_point(stat = "summary", fun = "mean") + 
            geom_line(stat = "summary", fun = "mean") + 
            labs(y = "average pc_acgt"))
        
        # Average pc_ambiguous per ISO week, coloured by run_name
        print(ggplot(data = code_df, mapping = aes(x = sequencing_submission_week, y = pc_ambiguous, colour = run_name)) + 
            geom_point(stat = "summary", fun = "mean") + 
            geom_line(stat = "summary", fun = "mean") + 
            labs(y = "average pc_ambiguous"))
        
        for (run in runs) {
            run_df <- filter(code_df, run_name == run)
            
            # Percentage acgt scatterplot
            p1 = ggplot(data = run_df, mapping = aes(x = pc_acgt, y = 0)) +
                geom_jitter(width = 0, alpha = alpha_val) +
                scale_y_continuous(breaks = NULL) +
                theme_bw() +
                theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())
            
            # Percentage acgt histogram
            p2 = ggplot(data = run_df, mapping = aes(x = pc_acgt)) +
                geom_histogram() + 
                theme_bw()
            
            # Percentage ambiguous scatterplot
            p3 = ggplot(data = run_df, mapping = aes(x = pc_ambiguous, y = 0)) +
                geom_jitter(width = 0, alpha = alpha_val) +
                scale_y_continuous(breaks = NULL) +
                theme_bw() +
                theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.y=element_blank())
            
            # Percentage ambiguous histogram
            p4 = ggplot(data = run_df, mapping = aes(x = pc_ambiguous)) +
                geom_histogram() +
                theme_bw()
            
            # Organise into grid
            print(grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2, top = run))
        }
        dev.off()
    }
}