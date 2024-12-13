#!/usr/bin/env Rscript

# Function that analyzes the debug output part of demultiplex_hto_and_scRNA.R script.
#
# Example:
# debug_demultiplex_hto_and_scRNA.R \
#    -i /gstore/project/tead/multinome/NGS4467/output/SAM24416357/debug/tables/TABLE_summary_hashedDrops_debug_per_HTO_SAM24416357.tsv \
#
# Roche/Genentech
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com

debug_demultiplex_hto_and_scRNA <- function(infile=infile, infile_per_HTO=infile_per_HTO, ngs_id=NULL, sample_id=NULL, outdir=outdir) {
    
    #########################################
    # parameters. Manually set for debug.   #
    # comment for real run                  @
    #########################################
    #num_threads = 1
    
    #hto_frs_id = "FRS17654"
    #arcseq_frs_id = "FRS14086"
    #outdir = "/gstore/project/tead/multinome/NGS4467_test/OUTPUT/"
    #test_ranks = NULL
    #rank = 20000
    #n_iters = 10000
    #skip_sample = NULL
    #debug_hashedDrops = TRUE
    #infile = "/gstore/project/tead/multinome/NGS4467/output/SAM24416357/debug/tables/TABLE_summary_hashedDrops_debug_SAM24416357.tsv"
    #infile_per_HTO = "/gstore/project/tead/multinome/NGS4467/output/SAM24416357/debug/tables/TABLE_summary_hashedDrops_debug_per_HTO_SAM24416357.tsv"
    #ngs_id = "NGS4467"
    #sample_id = "SAM24416357"
    #outdir = "/gstore/project/tead/multinome/NGS4467/output/SAM24416357/debug/"
    
    library(gridExtra)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(ggpubr)
    library(scales)
    
    hash_stats_debug_df = data.frame(fread(infile, header=T, sep="\t"), check.names=F)
    hash_stats_debug_df2 = melt(hash_stats_debug_df, id.vars=c("hashedDrops_config"))
    hash_stats_debug_df2 = hash_stats_debug_df2[hash_stats_debug_df2$variable %in% c("number_of_doublets", "number_of_confidents", "number_of_non_confidents_non_doublets"),]
    hash_stats_debug_df2$dmin = gsub(".*_dmin:(\\S+)_cnmads.*", "\\1", hash_stats_debug_df2$hashedDrops_config)
    hash_stats_debug_df2$cmin = gsub(".*_cmin:(\\S+).*", "\\1", hash_stats_debug_df2$hashedDrops_config)
    hash_stats_debug_df2$dnmads = gsub("^dnmads:(\\S+)_dmin.*", "\\1", hash_stats_debug_df2$hashedDrops_config)
    hash_stats_debug_df2$cnmads = gsub(".*_cnmads:(\\S+)_cmin.*", "\\1", hash_stats_debug_df2$hashedDrops_config)
    hash_stats_debug_df2$cnmads_cmin = paste0("cnmads:", hash_stats_debug_df2$cnmads, "_cmin:", hash_stats_debug_df2$cmin)
    hash_stats_debug_df2$dnmads_dmin = paste0("dnmads:", hash_stats_debug_df2$dnmads, "_dmin:", hash_stats_debug_df2$dmin)
    
    p_debug1 <- ggplot(data=hash_stats_debug_df2, aes(x=dnmads_dmin, y=value, fill=variable)) +
        facet_grid(cnmads_cmin ~ dnmads_dmin, space="free_x", scales="free_x") + 
        geom_bar(stat="identity", color="NA", position=position_dodge()) +
        scale_y_continuous(label=comma) +
        ylab("counts") +
        ggtitle(paste0(sample_id, " - ", ngs_id)) + 
        theme_minimal() + 
        theme(
            axis.text.x=element_text(angle=90, size=7, hjust=1),
            strip.text.y=element_text(angle=0, hjust=0),
            strip.text.x=element_text(angle=90, hjust=0.5)
        )
    pdf( file=paste0(outdir, "/FIGURE_summary_hashedDrops_debug_", sample_id, ".pdf"), height=10, width=10)
    print(p_debug1)
    dev.off()
    png( file=paste0(outdir, "/FIGURE_summary_hashedDrops_debug_", sample_id, ".png"), height=10, width=10, units="in", res=300)
    print(p_debug1)
    dev.off()
    
    # investigate the difference between high vs low doubets. Is it overrepresented by one HTO in particular?
    hash_stats_debug_per_HTO_df = data.frame(fread(infile_per_HTO, header=T, sep="\t"), check.names=F)
    hash_stats_debug_per_HTO_df = hash_stats_debug_per_HTO_df[hash_stats_debug_per_HTO_df$config %in% c(
                                                                            "dnmads:1_dmin:1_cnmads:1_cmin:1",
                                                                            "dnmads:3_dmin:1_cnmads:1_cmin:1",
                                                                            "dnmads:2_dmin:2_cnmads:2_cmin:1.5",
                                                                            "dnmads:2.5_dmin:1_cnmads:3_cmin:2"),]
    hash_stats_debug_per_HTO_df2 = hash_stats_debug_per_HTO_df[,c("number_of_HTOs_raw_freq","number_of_counts_per_HTO_raw","number_of_HTOs_filt_freq","number_of_counts_per_HTO_filt","config", "HTO_ID")]
    hash_stats_debug_per_HTO_df3 = melt(hash_stats_debug_per_HTO_df2, id.vars=c("HTO_ID", "config"))
    hash_stats_debug_per_HTO_df4 = hash_stats_debug_per_HTO_df3[hash_stats_debug_per_HTO_df3$HTO_ID != "Total",]
    p_debug2 <- ggplot(data=hash_stats_debug_per_HTO_df4, aes(x=HTO_ID, y=value)) +
        facet_grid(variable ~ config, scales="free_y") + 
        geom_bar(stat="identity", color="black", position=position_dodge()) +
        ggtitle(paste0("Counts per HTO; ", sample_id, " - ", ngs_id)) + 
        theme_bw() + theme(axis.text.x=element_text(angle=90, size=7, hjust=1), strip.text.y=element_text(angle=0) ) + scale_y_continuous(labels=comma)
    pdf( file=paste0(outdir, "/FIGURE_summary_hashedDrops_debug_per_HTO_", sample_id, ".pdf"), height=5, width=9)
    print(p_debug2)
    dev.off()
    png( file=paste0(outdir, "/FIGURE_summary_hashedDrops_debug_per_HTO_", sample_id, ".png"), height=5, width=9, units="in", res=300)
    print(p_debug2)
    dev.off()
}

usage=function(errM) {
    cat("\nUsage : Rscript debug_demultiplex_hto_and_scRNA.R [option] <Value>\n")
    cat("       -i        : Global debug output\n")
    cat("       -b        : HTO debug output per HTO file.\n")
    cat("       -n        : NGS ID (i.e. NGS[0-9]{4}\n")
    cat("       -s        : SAM ID (i.e. SAM[0-9]{?}\n")
    cat("       -o        : outdir location to write results\n")
}

ARG = commandArgs(trailingOnly = TRUE)

if(length(ARG) < 5) {
    usage("Potentially missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        infile=ARG[i+1]
    } else if (ARG[i] == "-b") {
        infile_per_HTO=ARG[i+1]
    } else if (ARG[i] == "-n") {
        ngs_id=ARG[i+1]
    } else if (ARG[i] == "-s") {
        sample_id=ARG[i+1]
    } else if (ARG[i] == "-o") {
        outdir=ARG[i+1]
    }
}

debug_demultiplex_hto_and_scRNA(infile=infile, infile_per_HTO=infile_per_HTO, ngs_id=ngs_id, sample_id=sample_id, outdir=outdir)
