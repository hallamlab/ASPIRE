#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
qc <- read.delim(args[1], check.names=FALSE)
qc$decontam_class <- factor(qc$decontam_class,
  levels=c("biological", "bio_control", "technical", "positive"))
summary <- do.call(rbind, lapply(levels(qc$decontam_class), function(cls) {
  x <- qc$post_qc_reads[qc$decontam_class == cls]
  data.frame(class=cls, samples=length(x), zero_reads=sum(x == 0),
    min=if(length(x)) min(x) else NA, median=if(length(x)) median(x) else NA,
    max=if(length(x)) max(x) else NA)
}))
write.table(summary, file.path(args[2], "read_depth_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
set.seed(1)
p <- ggplot(qc, aes(decontam_class, log10(post_qc_reads + 1), color=qc_status)) +
  geom_boxplot(aes(group=decontam_class), outlier.shape=NA, color="grey50") +
  geom_jitter(width=0.15, height=0) + theme_bw() +
  scale_x_discrete(labels=c(biological="Biological samples", bio_control="Biological controls",
                           technical="Technical controls", positive="Positive controls")) +
  labs(x="Sample class", y="log10(post-QC reads + 1)", color="QC status",
       caption="The depth inclusion cutoff applies only to biological samples; controls are exempt.")
ggsave(file.path(args[2], "read_depth_by_class.svg"), p, device=grDevices::svg, width=9, height=5)
