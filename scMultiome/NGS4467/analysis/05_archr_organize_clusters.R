#!/usr/bin/env Rscript

# Workflow template to help organize clusters, merging and naming them.
#
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

archr_organize_clusters <- function(indir=NULL, num_threads=num_threads, my_FDR=0.1, my_log2FC=0.5) {

    #setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    #indir = "/gstore/project/tead/scMultiome/NGS4467_pilot_JT/OUTPUT/archr/"
    #num_threads = 4
    #my_FDR = 0.05
    #my_log2FC = 1
    
    my_FDR = as.numeric(my_FDR)
    my_log2FC = as.numeric(my_log2FC)
    num_threads = as.numeric(num_threads)
    message("Using num_threads=", num_threads, ". log2FC=", my_log2FC, ". FDR=", my_FDR)
    if(is.null(indir)){
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
    source("./tools/R/lib/utils.R")
    
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
    # First, give a relevant name to each clusters of interest                 #
    # To do this, first plot side by side                                      #
    #    (umap colored by clusters) vs umap colored by samples                 #
    ############################################################################
    samples_vColors = vColors[1:length(unique(mapping[["Treatment_HTO"]]))]
    names(samples_vColors) = unique(mapping[["Treatment_HTO"]])
    
    embed = getEmbedding(archr_proj, embedding=paste0("UMAP_TileMatrix"))
    embed = cbind(embed, getCellColData(archr_proj, select=c("Treatment", "HTO", "Treatment_HTO", "Treatment_library", "library", "Clusters_TileMatrix")))
    colnames(embed)[1] = "UMAP_dim_1"
    colnames(embed)[2] = "UMAP_dim_2"
    
    my_labels_df = data.frame(Treatment_HTO=unique(embed$Treatment_HTO),label=paste0("", unique(embed$Treatment_HTO)))
    my_labels = embed %>%
        dplyr::group_by(Treatment_HTO) %>%
        dplyr::summarize(UMAP_dim_1 = median(UMAP_dim_1), UMAP_dim_2 = median(UMAP_dim_2)) %>%
        dplyr::left_join(my_labels_df) %>% as.data.frame()
    
    p_init_samples = NULL; p_init_samples = ggplot(embed, aes(x=UMAP_dim_1, y=UMAP_dim_2, color=Treatment_HTO)) +
        geom_point(size=0.3) +
        scale_color_manual(values=samples_vColors) +
        xlab(paste0("UMAP 1")) +
        ylab(paste0("UMAP 2")) +
        guides(color=guide_legend(ncol=2, override.aes=list(size=1.5))) +
        ggrepel::geom_label_repel(data=my_labels, box.padding=0.0, alpha=0.8, label.size=0.05, fill="gray70", color="black", size=3, aes(label=my_labels[["label"]])) +
        theme_minimal() +
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="bottom",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'), legend.title=element_blank())
    
    my_labels_df2 = data.frame(Clusters_TileMatrix=unique(embed$Clusters_TileMatrix),label=paste0("", unique(embed$Clusters_TileMatrix)))
    my_labels2 = embed %>%
        dplyr::group_by(Clusters_TileMatrix) %>%
        dplyr::summarize(UMAP_dim_1 = median(UMAP_dim_1), UMAP_dim_2 = median(UMAP_dim_2)) %>%
        dplyr::left_join(my_labels_df2) %>% as.data.frame()
    
    clusters_vColors = vColors20[1:length(unique(mapping[["Clusters_TileMatrix"]]))]
    names(clusters_vColors) = unique(mapping[["Clusters_TileMatrix"]])
    
    p_init_clusters = NULL; p_init_clusters = ggplot(embed, aes(x=UMAP_dim_1, y=UMAP_dim_2, color=.data[["Clusters_TileMatrix"]])) +
        geom_point(size=0.3) +
        scale_color_manual(values=clusters_vColors) +
        xlab(paste0("UMAP 1")) +
        ylab(paste0("UMAP 2")) +
        #ggtitle(paste0("UMAP; ", curr_signature)) +
        guides(color = guide_legend(override.aes=list(size=1.5))) + guides(fill=guide_legend(ncol=3)) + #guides(fill = guide_legend(override.aes = list(size=4))) +
        ggrepel::geom_label_repel(data=my_labels2, fill="gray70", box.padding=0.0, alpha=0.8, label.size=0.05, color="black", size=3, aes(label=my_labels2[["label"]])) +
        theme_minimal() +
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="bottom",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'), legend.title=element_blank())
    p1 = ggarrange(p_init_samples, p_init_clusters)

    # Define clusters names:
    # C1 = Sensitive_common_small;   Sensitive_side
    # C2 = Sensitive_no_agent_big    Sensitive_DMSO
    # C3 = Sensitive_agent;          Sensitive_GNE7883
    # C4 = Sensitive_agent;          Sensitive_GNE7883
    # C5 = Resistant_small;          Resistant_side
    # C6 = Resistant_big;            Resistant_main
    # C7 = Resistant_big;            Resistant_main 
    archr_proj$Clusters_TileMatrix_named = plyr::mapvalues(archr_proj$Clusters_TileMatrix,
                                    from = paste0("C_TileMatrix",1:7),
                                    to = c("Sensitive_side", "Sensitive_DMSO", "Sensitive_GNE7883", 
                                           "Sensitive_GNE7883",  "Resistant_side", "Resistant_main", 
                                           "Resistant_main")
    )
    
    # Then integrate it so that we have all 3 panels in the same plot.
    embed$Clusters_TileMatrix_named = getCellColData(archr_proj)$Clusters_TileMatrix_named
    my_labels_df3 = data.frame(Clusters_TileMatrix_named=unique(embed$Clusters_TileMatrix_named),label=paste0("", unique(embed$Clusters_TileMatrix_named)))
    my_labels3 = embed %>%
        dplyr::group_by(Clusters_TileMatrix_named) %>%
        dplyr::summarize(UMAP_dim_1 = median(UMAP_dim_1), UMAP_dim_2 = median(UMAP_dim_2)) %>%
        dplyr::left_join(my_labels_df3) %>% as.data.frame()
    
    clusters_named_vColors = vColors30[1:length(unique(getCellColData(archr_proj)[["Clusters_TileMatrix_named"]]))]
    names(clusters_named_vColors) = unique(mapping[["Clusters_TileMatrix_named"]])
    
    p_init_clusters_named = NULL; p_init_clusters_named = ggplot(embed, aes(x=UMAP_dim_1, y=UMAP_dim_2, color=Clusters_TileMatrix_named)) +
        geom_point(size=0.3) +
        scale_color_manual(values=clusters_named_vColors) +
        xlab(paste0("UMAP 1")) +
        ylab(paste0("UMAP 2")) +
        guides(color=guide_legend(ncol=3, override.aes=list(size=1.5))) +
        ggrepel::geom_label_repel(data=my_labels3, box.padding=0.0, alpha=0.8, label.size=0.05, fill="gray70", color="black", size=3, aes(label=my_labels3[["label"]])) +
        theme_minimal() +
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="bottom",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'), legend.title=element_blank())
 
    p2 = ggarrange(p_init_samples, p_init_clusters, p_init_clusters_named, nrow=1)
    ggsave(paste0(curr_outdir, "/marker_UMAP_TileMatrix_named.pdf"), plot=p1, device="pdf", height=5.5, width=16.5, units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/marker_UMAP_TileMatrix_named.png"), plot=p1, device="png", height=5.5, width=16.5, units="in", limitsize=FALSE)
    
    reference_ordinations = list("init_samples"=p_init_samples, "init_clusters"=p_init_clusters, "init_clusters_named"=p_init_clusters_named, "my_labels"=my_labels, "my_labels2"=my_labels2, "my_labels3"=my_labels3)
    saveRDS(reference_ordinations, file="../OUTPUT/reference_ordinations.rds")
    saveArchRProject(archr_proj, overwrite=TRUE, logFile=createLogFile("saveArchRProject"), threads=1)
    message("Completed archr_organize_clusters.R!")
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
archr_organize_clusters(indir=indir, num_threads=num_threads, my_FDR, my_log2FC)
