#!/usr/bin/env Rscript

# Function that runs gp.sc.archr::runArchR starting from a MultiAssayExperimentObject 
# object previously generated with 'merge_sc_and_archr.R'
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

archr_qc <- function(indir=NULL, mapping_file=NULL) {

    # For debugging:
    #setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    #indir = "../OUTPUT/archr/"
    #mapping_file = "../data/mapping_file.tsv"
    
    if(is.null(indir)){ # TODO check if dir
        stop("an indir has to be included")
    }
    
    if(is.null(mapping_file)){ # TODO check if dir
        stop("a mapping_file has to be included")
    }
    
    dir.create(paste0(indir, "/qc"))
    
    library(ArchR)
    library(MASS)
    library(viridis)
    library(scales)
    library(BSgenome.Hsapiens.Genentech.GRCh38)
    library(ggpubr)
    source("./utils.R")
    
    ggsave <- function(..., bg='white') ggplot2::ggsave(..., bg=bg)
    
    mapping = data.frame(fread(mapping_file, header=T, sep="\t"), check.names=F)
    row.names(mapping) = mapping$SampleID
    mapping$hash_assignment2 = mapping$SampleID
    
    archr_proj = ArchR::loadArchRProject(indir)
    
    # filter doublets. no output, maybe use sink()...?
    message("Filtering out doublets.")
    doublets_file = file("../OUTPUT/archr/qc/doublets.txt", open="wt")
    sink(doublets_file ,type="output")
    sink(doublets_file, type="message")
    archr_proj = filterDoublets(
        ArchRProj = archr_proj,
        cutEnrich = 1,
        cutScore = -Inf,
        filterRatio = 1
    )
    closeAllConnections()
    message("Done filtering out doublets.")
    
    x = readLines("../OUTPUT/archr/qc/doublets.txt")
    x = x[grep("SAM", x)]
    x = gsub("\t", "", x)
    doublets = paste(x, collapse="\n")
    
    # filter cells that do not contain rna. Keep info for table.
    number_of_cells = nrow(getCellColData(archr_proj))
    number_of_cells_with_RNA = nrow(getCellColData(archr_proj[!is.na(archr_proj$Gex_nUMI)]))
    number_of_cells_without_RNA = number_of_cells - number_of_cells_with_RNA
    archr_proj = archr_proj[!is.na(archr_proj$Gex_nUMI)]
    
    tmp_df = data.frame(getCellColData(archr_proj, select=c("hash_assignment2")))
    tmp_df2 = join(tmp_df, mapping, by="hash_assignment2")
    # Make sure hash_assignment2 are in the same order in both vectors.
    if(identical(tmp_df2$hash_assignment2, tmp_df$hash_assignment2) == FALSE){
        stop("something went wrong in joining metadata...")
    }
    
    ########################################
    # QC                                   #
    ########################################
    archr_proj$CNAME = tmp_df2$CNAME
    archr_proj$TEST_ARTICLE = tmp_df2$TEST_ARTICLE
    archr_proj$HTO = tmp_df2$HTO
    archr_proj$Treatment = paste0(archr_proj$CNAME, "_", archr_proj$TEST_ARTICLE)
    archr_proj$Treatment_HTO = paste0(archr_proj$CNAME, "_", archr_proj$TEST_ARTICLE, " (", archr_proj$HTO, ")")
    archr_proj$Treatment_library = paste0(archr_proj$Treatment, "_", archr_proj$library)
    
    df = getCellColData(archr_proj, select = c("log10(nFrags)", "TSSEnrichment", "hash_assignment2", "CNAME", "TEST_ARTICLE", "Treatment", "library", "HTO", "Treatment_HTO", "Treatment_library"))
    df$density = get_density(df[,1], df[,2], n = 100)
    
    # Compute density per sample.
    df2 = NULL
    for(i in 1:length(unique(df$Treatment_library))){
        curr_treatment = unique(df$Treatment_library)[i]
        df_tmp = df[df$Treatment_library == curr_treatment,]
        df_tmp$density = get_density(df_tmp[,1], df_tmp[,2], n = 100)
        
        if(i == 1){
            df2 = df_tmp
        }else{
            df2 = rbind(df2, df_tmp)
        }
    }
    
    message("Generating QC plots")
    p1 = ggplot(data=data.frame(df2), aes(x=df2[["log10(nFrags)"]], y=TSSEnrichment, color=density)) +
        geom_point(size=1) +
        xlim(c(log10(500), quantile(df2[,1], probs = 0.99))) +
        ylim(c(0, quantile(df2[,2], probs = 0.99))) +
        scale_color_viridis() +
        xlab("Log10(nFrags)") +
        facet_grid(library ~ Treatment_HTO) + 
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0)) +
        ggtitle("TSS enrichment score in function of fragment abundance (log10)") + 
        labs(caption="These plots show the ditribution of cells per sample (each point = 1 cell) according to their TSS score vs number of fragments. Density was computed separately for each sample (i.e. panel).\n In general, we should obtain a majority of cell having TSS score = ~12.5 and Log10(nFrags) = ~4.25")
    p1
    ggsave(paste0(indir, "/qc/density_scatter.png"), plot=p1, device="png", height=3.7, width=12.5, units="in")
    ggsave(paste0(indir, "/qc/density_scatter.pdf"), plot=p1, device="pdf", height=3.7, width=12.5, units="in")
    
    p2 = ggplot(data=data.frame(df2), aes(x=TSSEnrichment)) + #, color=HTO)) +
        geom_histogram(binwidth=0.5) +
        facet_grid(library ~ Treatment_HTO) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0)) +
        ggtitle("TSS enrichment score histograms") + 
        xlab("TSS enrichment score") +
        labs(caption="These plots show the distribution count of TSS enrichment scores per sample.")
    p2
    ggsave(paste0(indir, "/qc/TSS_enrichment_score_count_histo.png"), plot=p2, device="png", height=3.7, width=12.5, units="cm")
    ggsave(paste0(indir, "/qc/TSS_enrichment_score_count_histo.pdf"), plot=p2, device="pdf", height=3.7, width=12.5, units="cm")
    
    p3 = ggplot(data=data.frame(df2), aes(x=TSSEnrichment)) + #, color=HTO, fill=HTO)) +
        geom_density(alpha=1) +
        facet_grid(library ~ Treatment_HTO) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0)) +
        ggtitle("TSS enrichment score density plots") + 
        xlab("TSS enrichment score") +
        labs(caption="Plots showing the distribution density of TSS enrichment scores per sample.")
    p3
    ggsave(paste0(indir, "/qc/TSS_enrichment_score_count_density.png"), plot=p3, device="png", height=9.5, width=30.3, units="cm")
    ggsave(paste0(indir, "/qc/TSS_enrichment_score_count_density.pdf"), plot=p3, device="pdf", height=9.5, width=30.3, units="cm")
    
    p4 = ggplot(data=data.frame(df2), aes(y=TSSEnrichment, x=hash_assignment2, color=library)) +
        geom_violin() +
        geom_point(aes(fill=library), alpha=0.2, size=0.5, pch = 21, position=position_jitterdodge(0.2)) +
        facet_grid(. ~ Treatment_HTO, scales="free_x", space="free_x") +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.text.x=element_text(angle=90), strip.text.y=element_text(angle=0)) +
        ggtitle("TSS enrichment score violin plots") + 
        xlab("") +
        labs(caption="Violin plots showing the distribution of TSS enrichment scores per sample. Each point = 1 cell")
    p4
    ggsave(paste0(indir, "/qc/TSS_enrichment_score_count_density.png"), plot=p4, device="png", height=15.6, width=24, units="cm")
    ggsave(paste0(indir, "/qc/TSS_enrichment_score_count_density.pdf"), plot=p4, device="pdf", height=15.6, width=24, units="cm")
    
    p5 = ggplot(data=data.frame(df2), aes(x=df2[["log10(nFrags)"]])) + #, color=HTO)) +
        geom_histogram(binwidth=0.025) +
        facet_grid(library ~ Treatment_HTO) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0)) +
        xlab("log10(nFrags)") + 
        ggtitle("TSS enrichment score histograms") + 
        labs(caption="Plots showing the distribution count of the number of fragments (log10 scale) per sample.")
    p5
    ggsave(paste0(indir, "/qc/log10nfrags_histo.png"), plot=p5, device="png", height=10.55, width=30, units="cm")
    ggsave(paste0(indir, "/qc/log10nfrags_histo.pdf"), plot=p5, device="pdf", height=10.55, width=30, units="cm")
    
    p6 = ggplot(data=data.frame(df2), aes(x=df2[["log10(nFrags)"]])) + # , color=HTO, fill=HTO)) +
        geom_density(alpha=1) +
        facet_grid(library ~ Treatment_HTO) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0)) +
        ggtitle("log10(nFrags) score density plots") + 
        xlab("log10(nFrags)") +
        labs(caption="Plots showing the distribution density of the number of fragments (log10 scale) per sample.")
    p6
    ggsave(paste0(indir, "/qc/log10nfrags_density.png"), plot=p6, device="png", height=10.55, width=30, units="cm")
    ggsave(paste0(indir, "/qc/log10nfrags_density.pdf"), plot=p6, device="pdf", height=10.55, width=30, units="cm")
    
    p7 = ggplot(data=data.frame(df2), aes(y=df2[["log10(nFrags)"]], x=hash_assignment2, color=library)) +
        geom_violin() +
        geom_point(aes(fill=library), alpha=0.2, size=0.5, pch = 21, position=position_jitterdodge(0.2)) +
        facet_grid(. ~ Treatment_HTO, scales="free_x", space="free_x") +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.text.x=element_text(angle=90), strip.text.y=element_text(angle=0)) +
        ggtitle("Log10(nFrags) score violin plots") + 
        xlab("") +
        ylab("log10(nFrags)") +
        labs(caption="Violin plots showing the distribution of the number of fragments in log10 scale per sample. Each point = 1 cell")
    p7    
    ggsave(paste0(indir, "/qc/log10nfrags_violin.png"), plot=p7, device="png", height=10.55, width=30, units="cm")
    ggsave(paste0(indir, "/qc/log10nfrags_violin.pdf"), plot=p7, device="pdf", height=10.55, width=30, units="cm")
    
    # compute fragment size.
    chrom_sizes = getChromSizes(archr_proj)
    chr = paste0(seqnames(chrom_sizes))
    groups = getCellColData(archr_proj, select="hash_assignment2", drop=FALSE)
    uniq_groups = gtools::mixedsort(unique(groups[,1]))
    max_size = 750
    
    message("Computing fragment sizes.")
    dfFS = lapply(seq_along(uniq_groups), function(x){
        cellx = rownames(groups)[which(paste0(groups[,1]) == uniq_groups[x])]
        
        for(i in seq_along(chr)){
            if(i == 1){
                fsi = unlist(suppressMessages(getFragmentsFromProject(
                    ArchRProj = archr_proj,
                    subsetBy = chrom_sizes[paste0(seqnames(chrom_sizes)) %in% chr[i]],
                    cellNames = cellx,
                )), use.names=FALSE) %>% width %>% tabulate(nbins = max_size)
            }else{
                fsi = fsi + unlist(suppressMessages(getFragmentsFromProject(
                    ArchRProj = archr_proj,
                    subsetBy = chrom_sizes[paste0(seqnames(chrom_sizes)) %in% chr[i]],
                    cellNames = cellx,
                )), use.names=FALSE) %>% width %>% tabulate(nbins = max_size)
            }
        }
        
        df = DataFrame(
            group = uniq_groups[x],
            fragmentSize = seq_along(fsi),
            fragmentNumber = fsi,
            fragmentPercent = round(100*fsi/sum(fsi),4)
        )
        df
    }) %>% Reduce("rbind", .)
    
    dfFS2 = data.frame(dfFS)
    dfFS2$SAM_ID = gsub("^(SAM\\d+)_.*", "\\1", dfFS2$group)
    dfFS2$HTO = gsub("^SAM\\d+_*", "", dfFS2$group)
    dfFS2 = merge(dfFS2, mapping[,c("hash_assignment2", "CNAME", "TEST_ARTICLE")], by.x="group", by.y="hash_assignment2")
    dfFS2$Treatment = paste0(dfFS2$CNAME, "_", dfFS2$TEST_ARTICLE, "_", dfFS2$HTO)
    
    p8 = ggplot(dfFS2, aes(fragmentSize, fragmentPercent)) +
        geom_line(size=0.5) +
        facet_grid(SAM_ID ~ Treatment) +
        xlab("ATAC-seq Fragment Size (bp)") +
        ylab("Percentage of Fragments") +
        ggtitle("Percentages of fragments in function of fragment sizes.") +
        scale_y_continuous(limits = c(0, max(dfFS2$fragmentPercent)*1.05), expand = c(0,0)) +
        scale_x_continuous(limits = c(min(dfFS2$fragmentSize), max(dfFS2$fragmentSize)), expand = c(0,0)) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour="black", linewidth=0.5), strip.text.y=element_text(angle=0))
    p8
    ggsave(paste0(indir, "/qc/fragment_size_distrib_perc.png"), plot=p8, device="png", height=10, width=30, units="cm")
    ggsave(paste0(indir, "/qc/fragment_size_distrib_perc.pdf"), plot=p8, device="pdf", height=10, width=30, units="cm")
    
    p9 = ggplot(dfFS2, aes(fragmentSize, fragmentNumber)) +
        geom_line(size=0.5) +
        facet_grid(SAM_ID ~ Treatment) +
        xlab("ATAC-seq Fragment Size (bp)") +
        ylab("Number of Fragments") +
        ggtitle("Number of fragments in function of fragment sizes.") +
        scale_y_continuous(label=comma, limits = c(0, max(dfFS2$fragmentNumber)*1.05), expand = c(0,0)) +
        scale_x_continuous(limits = c(min(dfFS2$fragmentSize), max(dfFS2$fragmentSize)), expand = c(0,0)) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0))
    p9
    ggsave(paste0(indir, "/qc/fragment_size_distrib_count.png"), plot=p9, device="png", height=10, width=30, units="cm")
    ggsave(paste0(indir, "/qc/fragment_size_distrib_count.pdf"), plot=p9, device="pdf", height=10, width=30, units="cm")
    # TODO print the plot with a large width to appreciate the patterns and differences between SAM_IDs.
    
    message("Computing TSS enrichment scores")
    dfTSS = plotTSSEnrichment(ArchRProj = archr_proj,
                           groupBy = "hash_assignment2",
                           returnDF = TRUE)
    head(dfTSS)
    dfTSS2 = data.frame(dfTSS)
    dfTSS2$SAM_ID = gsub("^(SAM\\d+)_.*", "\\1", dfTSS2$group)
    dfTSS2$HTO = gsub("^SAM\\d+_*", "", dfTSS2$group)
    dfTSS2 = merge(dfTSS2, mapping[,c("hash_assignment2", "CNAME", "TEST_ARTICLE")], by.x="group", by.y="hash_assignment2")
    dfTSS2$Treatment = paste0(dfTSS2$CNAME, "_", dfTSS2$TEST_ARTICLE, "_", dfTSS2$HTO)
    
    p10 = ggplot(dfTSS2, aes(x, smoothValue)) +
        geom_line(size = 0.5) +
        facet_grid(SAM_ID ~ Treatment) +
        xlab("Distance to TSS") +
        ylab("Percentage of aggregate TSS enrichment score values") +
        ggtitle("Percentage of TSS enrichment values in function of distance to TSS") +
        scale_y_continuous(limits = c(0, max(dfTSS2$smoothValue)*1.05), expand = c(0,0)) +
        scale_x_continuous(limits = c(min(dfTSS2$x), max(dfTSS2$x)), expand = c(0,0)) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), strip.text.y=element_text(angle=0), axis.text.x=element_text(angle=90))
    p10
    ggsave(paste0(indir, "/qc/TSS_score_distrib_perc.png"), plot=p10, device="png", height=10, width=30, units="cm")
    ggsave(paste0(indir, "/qc/TSS_score_distrib_perc.pdf"), plot=p10, device="pdf", height=10, width=30, units="cm")
    
    p11 = ggplot(dfTSS2, aes(x, value)) +
        geom_line(size = 0.5) +
        facet_grid(SAM_ID ~ Treatment) +
        xlab("Distance to TSS") +
        ylab("Aggregate TSS enrichment score values") +
        ggtitle("Number of TSS enrichment values in function of distance to TSS") +
        scale_y_continuous(label=comma, limits = c(0, max(dfTSS2$value)*1.05), expand = c(0,0)) +
        scale_x_continuous(limits = c(min(dfTSS2$x), max(dfTSS2$x)), expand = c(0,0)) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.text.x=element_text(angle=90), strip.text.y=element_text(angle=0))
    p11
    ggsave(paste0(indir, "/qc/TSS_score_distrib_count.png"), plot=p11, device="png", height=10, width=30, units="cm")
    ggsave(paste0(indir, "/qc/TSS_score_distrib_count.pdf"), plot=p11, device="pdf", height=10, width=30, units="cm")
    
    # Create table of total counts (TSS and fragments).
    dfFS3 = dfFS2 %>% 
        dplyr::group_by(group) %>%
        dplyr::summarize_at(vars(fragmentNumber), list(sum = ~ sum(.))) %>%
        as.data.frame()
    dfFS3$perc = round(dfFS3$sum / sum(dfFS3$sum) * 100, digits=2)
    colnames(dfFS3) = c("group", "fragment_size_count", "fragment_size_perc")
    
    dfTSS3 = dfTSS2 %>% 
        dplyr::group_by(group) %>%
        dplyr::summarize_at(vars(value), list(sum = ~ sum(.))) %>%
        as.data.frame()
    dfTSS3$perc = round(dfTSS3$sum / sum(dfTSS3$sum) * 100, digits=2)
    colnames(dfTSS3) = c("group", "TSS_count", "TSS_perc")
    
    summary_table = merge(dfFS3, dfTSS3, by="group")
    summary_table$fragment_size_count = comma(summary_table$fragment_size_count)
    summary_table$TSS_count = comma(summary_table$TSS_count)
    
    table_theme = ttheme("lBlack", base_size=7, padding = unit(c(1, 1), "mm"),
        tbody.style = tbody_style(
            hjust=1, x=0.9, size=7
        )
    )
    tab = ggtexttable(summary_table, theme=table_theme)
    tab = tab_add_title(tab, str_wrap("Fragment and TSS count summary.", width=30))
    tab = tab_add_footnote(tab, paste0("Doublets that were filtered out (i.e. not present in table):\n ", doublets), size=10)
    tab = tab_add_footnote(tab, paste0("Number of cells without RNA (i.e. not present in table) :\n ", number_of_cells_without_RNA), size=10)
    tab2 = tab %>%
        tab_add_vline(at.column=2, column.side="left", from.row=3, linetype=1, to.row=nrow(summary_table))
    ggsave(paste0(indir, "/qc/summary_table.png"), plot=tab2, device="png", height=15, width=20, units="cm")
    ggsave(paste0(indir, "/qc/summary_table.pdf"), plot=tab2, device="pdf", height=15, width=20, units="cm")
    
    # Save project
    saveArchRProject(archr_proj, overwrite=TRUE, logFile=createLogFile("saveArchRProject"), threads=1)
    
    message("archr_qc.R completed!")
}

usage=function(errM) {
    cat("\nUsage : Rscript my_script_name.R [option] <Value>\n")
    cat("       -i        : indir containing the archR data structure.")
    cat("       -m        : mapping_file containing relevant metadata for the project.")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
    usage("missing arguments")
    stop("missing -i <string> or -m <string> argument.")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        indir=ARG[i+1]
    }else if(ARG[i] == "-m") {
        mapping_file=ARG[i+1]
    }
}
archr_qc(indir=indir, mapping_file=mapping_file)
