#!/usr/bin/env Rscript

# Function that runs DEGs analyses for a scMultiome project.
# object previously generated with 'archr_call_peaks.R'
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

archr_markers_analysis_pairwise <- function(indir=NULL, num_threads=num_threads, my_FDR=0.1, my_log2FC=0.5) {

    setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    indir = "/gstore/project/tead/scMultiome/NGS4467_pilot_JT/OUTPUT/archr/"
    num_threads = 4
    my_FDR = 0.05
    my_log2FC = 1
    
    my_FDR = as.numeric(my_FDR)
    my_log2FC = as.numeric(my_log2FC)
    num_threads = as.numeric(num_threads)
    message("Using num_threads=", num_threads, ". log2FC=", my_log2FC, ". FDR=", my_FDR)
    if(is.null(indir)){ # TODO check if dir
        stop("an indir has to be included")
    }
    
    if(is.null(num_threads)){
        stop("num_threads=<positive_integer> has to be incuded.")
    }else{
        if(!is.numeric(num_threads)){
            stop("num_threads has to be numeric.")
        }
    }
        
    library(ArchR)
    library(parallel)
    library(BSgenome.Hsapiens.Genentech.GRCh38)
    library(gridExtra)
    library(genomitory)
    library(ggpubr)
    library(grid)
    library(ggridges)
    library(ggrepel)
    library(epiregulon)
    source("./tools/R/lib/utils.R")
   
    ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)

    archr_proj = ArchR::loadArchRProject(indir)
    getAvailableMatrices(archr_proj)
    
    curr_outdir = paste0(indir, "/Plots/")
    dir.create(curr_outdir)
    
    # Loading metadata
    mapping = getCellColData(archr_proj, select = c("hash_assignment2", "CNAME", "TEST_ARTICLE", "Treatment", "library", "HTO", "Treatment_HTO", "Treatment_library", "Clusters_TileMatrix"))
    mapping$Treatment_library_HTO = paste0(mapping$Treatment_library, "_", mapping$HTO)
    archr_proj$Treatment_library_HTO = mapping$Treatment_library_HTO
    mapping = unique(mapping)
    mapping = data.frame(mapping)
    
    ############################################################################
    # Check genes of interest regardless of if they are sig                    #
    # Regardless of if signature genes were found to be significant or not.    # 
    # Plot genes of interest.                                                  #
    ############################################################################
    reference_ordinations = readRDS(file="../OUTPUT/reference_ordinations.rds")
   
    my_marker_genes = c("TEAD1", "YAP1", "FOSL1")
    my_marker_tfs =  c("TEAD1", "YAP", "FOSL1")
    # define use and control groups 
    use_groups =   c("Resistant_main", "Resistant_main",    "Sensitive_GNE7883")
    use_controls = c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_DMSO")
    
    differential_peaks = list()
    enrichments = list()
    enrichment_dfs = list()
    ps = list()
    df_enrichment = NULL
    for(i in 1:length(use_groups)){
        curr_use = use_groups[i]
        curr_control = use_controls[i]
        curr_differential_peaks = getMarkerFeatures(archr_proj, useMatrix="PeakMatrix", groupBy="Clusters_TileMatrix_named", bias=c("TSSEnrichment", "log10(nFrags)"), 
                                                    testMethod="wilcoxon", logFile="differentialpeaks", useGroups=curr_use, bgdGroups=curr_control)
        
        tmp_up = getMarkers(curr_differential_peaks, cutOff = paste0("FDR <= ", my_FDR, " & Log2FC >= ", my_log2FC))
        tmp_down = getMarkers(curr_differential_peaks, cutOff =  paste0("FDR <= ", my_FDR, " & Log2FC <= -", my_log2FC))
        
        differential_peaks[["Clusters_TileMatrix_named"]][[paste0(curr_use, "_vs_", curr_control)]] = curr_differential_peaks
        
        message("computing motif enrichment")
        cut_off = c()
        cut_off["Up"] = paste0("FDR <= ", my_FDR, " & Log2FC >= ", my_log2FC)
        cut_off["Down"] = paste0("FDR <= ", my_FDR, " & Log2FC <= -", my_log2FC)
        
        for (peak_annotation in c("Motif", "ENCODE_and_cistromeDB_TF_peaks")){
            for (direction in c("Up", "Down")) {
                curr_enrichment = peakAnnoEnrichment(
                    seMarker = curr_differential_peaks,
                    archr_proj,
                    peakAnnotation=peak_annotation,
                    cutOff=cut_off[direction],
                    logFile = "differentialpeaks"
                )
                
                df = data.frame(TF = rownames(curr_enrichment), mlog10Padj = assay(curr_enrichment)[,1])
                df = df[order(df$mlog10Padj, decreasing = TRUE),]
                df$rank <- seq_len(nrow(df))
                label = paste0("Clusters_TileMatrix_named", "\n", curr_use, "_vs_", curr_control, "\n", peak_annotation, ";", direction, " ;", peak_annotation)
                df$comparison = paste0(curr_use,"_vs_", curr_control)
                df$direction = "Up"
                df$annotation = peak_annotation
                
                if(is.null(df_enrichment)){
                    df_enrichment = df
                }else{
                    df_enrichment = rbind(df_enrichment, df)
                }
                
                ggUp = ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
                    geom_point(size = 1) +
                    ggtitle(label) + 
                    ggrepel::geom_label_repel(
                        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
                        size = 1.8,
                        nudge_x = 0.5, max.overlaps=15, label.r=0.005, label.padding=0.09,
                        color = "black"
                    ) + theme_ArchR() + 
                    ylab("-log10(P-adj) Motif Enrichment") + 
                    xlab("Rank Sorted TFs Enriched") +
                    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
                    theme_ArchR() + theme(axis.title=element_text(size=10), axis.text.x=element_text(size=9), plot.title=element_text(size=6),
                                          legend.title=element_text(size=7), legend.text=element_text(size=8, angle=90, hjust=1), legend.key.size=unit(0.25, 'cm'))
                
                ps[[paste0(curr_use,"_vs_", curr_control)]][[peak_annotation]][[direction]] = ggUp
            }
        }
    }

    head(df_enrichment)
    df_enrichment$mlog10Padj = ifelse(df_enrichment$direction == "Down", df_enrichment$mlog10Padj * -1, df_enrichment$mlog10Padj * 1)
    df_enrichment2_motif_mean = reshape2::dcast(df_enrichment[df_enrichment$annotation == "Motif",], TF ~ comparison, value.var="mlog10Padj", 
                                     fun.aggregate=mean)
    write.table(df_enrichment2, "../OUTPUT/peaks_comparisons_enrichment_motifs.tsv", row.names=F, sep="\t", quote=F)
    
    # Here write tables for end user
    cbind.fill <- function(...){
        nm = list(...)
        nm = lapply(nm, as.matrix)
        n = max(sapply(nm, nrow))
        do.call(cbind, lapply(nm, function (x)
            rbind(x, matrix(, n-nrow(x), ncol(x)))))
    }

    for(annotation in c("Motif", "ENCODE_and_cistromeDB_TF_peaks")){
        my_final_df=NULL
        for(comp in unique(names(enrichment_dfs))){
            header = colnames(enrichment_dfs[[comp]][[annotation]][["Up"]])
            tmp = enrichment_dfs[[comp]][[annotation]][["Up"]]
            colnames(tmp) = rep(comp, 3)
            up = rbind(rep("Up", 3), header, tmp)

            header = colnames(enrichment_dfs[[comp]][[annotation]][["Down"]])
            tmp = enrichment_dfs[[comp]][[annotation]][["Down"]]
            colnames(tmp) = rep(comp, 3)
            down = rbind(rep("Down", 3), header, tmp)

            up_down = cbind(up, down)

            if(is.null(my_final_df)){
                my_final_df = up_down
            }else{
                my_final_df = cbind.fill(my_final_df, up_down)
            }
        }
        write.table(my_final_df, paste0(curr_outdir, "./peaks_enrichment_", annotation, ".csv"), sep=",", row.names=F)
    }
           
    for(comparison in names(ps)){
        p = ggarrange(
            ps[[comparison]][["Motif"]][["Up"]],                          ps[[comparison]][["Motif"]][["Down"]],
            ps[[comparison]][["ENCODE_and_cistromeDB_TF_peaks"]][["Up"]], ps[[comparison]][["ENCODE_and_cistromeDB_TF_peaks"]][["Down"]],
            ncol=2, nrow=2
        )
        #p
        ggsave(paste0(curr_outdir, "/ranks_named_clusters_", comparison, ".pdf"), plot=p, device="pdf", height=9, width=6, units="in", limitsize=FALSE)
    }
    message("Completed archr_differential_peaks.R!")
}

usage=function(errM) {
    cat("\nUsage : Rscript my_script_name.R [option] <Value>\n")
    cat("       -i        : indir containing the archR data structure.")
    cat("       -n        : num_threads.")
    cat("       -s        : FDR cutoff. Default=0.1")
    cat("       -f        : log2FC. Default=0.5")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
    usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        indir=ARG[i+1]
    }else if(ARG[i] == "-n") {
        num_threads=ARG[i+1]
    }else if(ARG[i] == "-s") {
        my_FDR=ARG[i+1]
    }else if(ARG[i] == "-f") {
        my_log2FC=ARG[i+1]
    }
    
}
archr_differential_peaks(indir=indir, num_threads=num_threads, my_FDR, my_log2FC)
