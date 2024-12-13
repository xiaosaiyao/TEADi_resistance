#!/usr/bin/env Rscript

# Function that runs DEGs analyses for a scMultiome project.
# object previously generated with 'archr_call_peaks.R'
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

archr_marker_genes <- function(indir=NULL, num_threads=num_threads, my_FDR=0.1, my_log2FC=0.5) {

    #setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    #indir = "/gstore/project/tead/scMultiome/NGS4467_pilot_JT/OUTPUT/archr/"
    #num_threads = 4
    #my_FDR = 0.05
    #my_log2FC = 1
    
    num_threads = as.numeric(num_threads)
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
    
    clusters_named_vColors = c("Resistant_main" = "#000000", 
                               "Resistant_side" ="#696969", 
                               "Sensitive_DMSO" = "#8b4513", 
                               "Sensitive_GNE7883" = "#006400",
                               "Sensitive_side" = "#808000")
    
    ############################################################################
    # Check genes of interest regardless of if they are sig                    #
    # Regardless of if signature genes were found to be significant or not.    # 
    # Plot genes of interest.                                                  #
    # YAP1, TEAD1, FOSL1                                                       #
    ############################################################################
    reference_ordinations = readRDS(file="../OUTPUT/reference_ordinations.rds")
    
    my_marker_genes = c("TEAD1", "YAP1", "FOSL1", "SATB2")
    p_gex1 = plotEmbedding(archr_proj, pal = c("gray", "lightblue", "blue"), colorBy="GeneExpressionMatrix", name=my_marker_genes, embedding="UMAP_TileMatrix", quantCut=c(0.1, 0.90), imputeWeights=NULL, plotAs="points", size=0.7)
    
    # get ChromVAR MotifMatrix
    my_marker_tfs = c("YAP", "TEAD", "FOSL", "SATB2")
    plotVarDev = getVarDeviations(archr_proj, name="MotifMatrix", plot=FALSE)
    
    markersMM = ArchR::getFeatures(archr_proj, select=paste(my_marker_tfs, collapse="|"), useMatrix="MotifMatrix")
    markersMM = grep("z:", markersMM, value=TRUE)
    markersMM = gsub("z:", "", markersMM)
    
    # get ChromVAR TFPeaksDeviationsMatrix
    plotVarDev = getVarDeviations(archr_proj, name="TFPeaksDeviationsMatrix", plot=FALSE)
    plotPDF(plotVarDev, name = "Variable-TFPeaksDeviationsMatrix-Deviation-Scores", width=10, height=10, ArchRProj=archr_proj, addDOC=FALSE)
    markersTF = ArchR::getFeatures(archr_proj, select=paste(my_marker_tfs, collapse="|"), useMatrix = "TFPeaksDeviationsMatrix")
    markersTF = grep("z:", markersTF, value=TRUE)
    markersTF = gsub("z:", "", markersTF)
    
    markersTFBs = c("TEAD1", "YAP1", "FOSL1", "SATB2")
    
    ############################################################################
    #   Check if genes of interest are sig                                     #
    ############################################################################
   
    # Plot signatures for GeneExpressionMatrix, GeneScoreMatrix and Activity(epiregulon)
    seGEX = getMatrixFromProject(archr_proj, useMatrix="GeneExpressionMatrix")
    seMM = getMatrixFromProject(archr_proj, useMatrix="MotifMatrix")
    seTF = getMatrixFromProject(archr_proj, useMatrix="TFPeaksDeviationsMatrix")
    seTFB = getMatrixFromProject(archr_proj, useMatrix="TF_bindingMatrix")
    
    # GeneExpressionMatrix
    idx = which(rowData(seGEX)$name %in% my_marker_genes)
    seGEX2 = seGEX[idx,]
    embed_gex = getEmbedding(archr_proj,  embedding="UMAP_TileMatrix")
    selected_genes_cells = data.frame(as.matrix(t(assays(seGEX2)$GeneExpressionMatrix)))
    colnames(selected_genes_cells) = rowData(seGEX2)$name
    identical(row.names(embed_gex), row.names(selected_genes_cells))
    embed_gex = cbind(embed_gex, selected_genes_cells)
    embed_gex$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
    colnames(embed_gex)[1] = "UMAP_dim_1"
    colnames(embed_gex)[2] = "UMAP_dim_2"
    
    # MotifMatrix
    idx = which(rowData(seMM)$name %in% markersMM)
    seMM2 = seMM[idx,]
    embed_MM = getEmbedding(archr_proj,  embedding="UMAP_TileMatrix")
    selected_genes_cells = data.frame(as.matrix(t(assays(seMM2)$z)))
    colnames(selected_genes_cells) = rowData(seMM2)$name
    identical(row.names(embed_MM), row.names(selected_genes_cells))
    embed_MM = cbind(embed_MM, selected_genes_cells)
    embed_MM$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
    colnames(embed_MM)[1] = "UMAP_dim_1"
    colnames(embed_MM)[2] = "UMAP_dim_2"
    
    # TF
    idx = which(rowData(seTF)$name %in% markersTF)
    seTF2 = seTF[idx,]
    embed_TF = getEmbedding(archr_proj, embedding="UMAP_TileMatrix")
    selected_genes_cells = data.frame(as.matrix(t(assays(seTF2)$z)))
    colnames(selected_genes_cells) = rowData(seTF2)$name
    identical(row.names(embed_TF), row.names(selected_genes_cells))
    embed_TF = cbind(embed_TF, selected_genes_cells)
    embed_TF$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
    colnames(embed_TF)[1] = "UMAP_dim_1"
    colnames(embed_TF)[2] = "UMAP_dim_2"
    
    # TFB
    idx = which(rowData(seTFB)$name %in% markersTFBs)
    seTFB2 = seTFB[idx,]
    embed_TFB = getEmbedding(archr_proj,  embedding="UMAP_TileMatrix")
    selected_genes_cells = data.frame(as.matrix(t(assays(seTFB2)$z)))
    colnames(selected_genes_cells) = rowData(seTFB2)$name
    identical(row.names(embed_TFB), row.names(selected_genes_cells))
    embed_TFB = cbind(embed_TFB, selected_genes_cells)
    embed_TFB$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
    colnames(embed_TFB)[1] = "UMAP_dim_1"
    colnames(embed_TFB)[2] = "UMAP_dim_2"
    
    # Generate plots
    ##########################
    # GEX                    #
    ##########################
    embed_gex2 = NULL
    local(for(k in 1:length(my_marker_genes)){
        curr_marker_gene = my_marker_genes[k]
        tmp_df = embed_gex[,c("UMAP_dim_1","UMAP_dim_2", "Clusters_TileMatrix_named")]
        tmp_df = cbind(tmp_df, embed_gex[[curr_marker_gene]])
        colnames(tmp_df)[ncol(tmp_df)] = "value"
        tmp_df$symbol = curr_marker_gene
        if(is.null(embed_gex2)){
            embed_gex2 <<- tmp_df
        }else{
            embed_gex2 <<- rbind(embed_gex2, tmp_df)
        }
    })
    embed_gex2$Clusters_TileMatrix_named = factor(embed_gex2$Clusters_TileMatrix_named, levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))

    my_comparisons = list(
        c("Sensitive_DMSO", "Sensitive_GNE7883"),
        c("Sensitive_DMSO", "Resistant_main"),
        c("Resistant_main", "Sensitive_GNE7883")
    )
    
    ps_gex = list()  #umap
    ps_gex2 = list() #violin
    local(
    for(gene in my_marker_genes){
        p_gex = NULL; p_gex = ggplot(embed_gex2[embed_gex2$symbol == gene,], aes(x=UMAP_dim_1, y=UMAP_dim_2, color=log2(value + 1)), fill="black") +
            geom_point(size=0.3, alpha=1) +
            scale_color_gradientn(colors=colorRampPalette(c("gray95", "gray40", "blue"))(100), labels=function(x) paste0(format(round(x, 2), nsmall = 2))) +
            ggtitle(paste0(gene)) +
            xlab(paste0("UMAP dim 1")) +
            ylab(paste0("UMAP dim 2")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
                  axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
                  legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        
        ps_gex[[gene]] <<- p_gex
       
        embed_gex2 = embed_gex2[!embed_gex2$Clusters_TileMatrix_named %in% c("Resistant_side", "Sensitive_side"),]
        embed_gex2$Clusters_TileMatrix_named = factor(embed_gex2$Clusters_TileMatrix_named, unique(embed_gex2$Clusters_TileMatrix_named)) 
        p_gex2 = NULL; p_gex2 = ggplot(embed_gex2[embed_gex2$symbol == gene,], aes(x=Clusters_TileMatrix_named, y=log2(value + 1), fill=Clusters_TileMatrix_named), fill="black") +
            geom_violin() +  geom_boxplot(outlier.shape=NA, width=0.1, color="grey70") +
            scale_fill_manual(values=clusters_named_vColors) +
            ggtitle(paste0(gene)) +
            ylab(paste0("gex (log2)")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            stat_compare_means(comparisons=my_comparisons, method="wilcox.test", paired=F, size=3) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
                  axis.text=element_text(size=8, angle=90, hjust=1, vjust=0.5), strip.text.y=element_text(size=8, angle=0), legend.position="right",
                  legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        ps_gex2[[gene]] <<- p_gex2
    }
    )
    saveRDS(file="../OUTPUT/gex_ordinations_YAP_FOSL_TEAD.rds", ps_gex)
    saveRDS(file="../OUTPUT/gex_violin_YAP_FOSL_TEAD.rds", ps_gex2)
    p_gex_all = ggarrange(plotlist=ps_gex)
    p_gex_all2 = ggarrange(plotlist=ps_gex2)
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_UMAP_TileMatrix_gene_expression", ".pdf"), plot=p_gex_all, device="pdf", height=5.48, width=10.12, units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_violin_TileMatrix_gene_expression", ".pdf"), plot=p_gex_all2, device="pdf", height=5.48, width=10.12, units="in", limitsize=FALSE)
    
    #############################
    # chromvar TF_bindingMatrix #
    #############################
    embed_TFB2 = NULL
    local(for(k in 1:length(markersTFBs)){
        curr_marker_gene = markersTFBs[k]
        tmp_df = embed_TFB[,c("UMAP_dim_1","UMAP_dim_2", "Clusters_TileMatrix_named")]
        tmp_df = cbind(tmp_df, embed_TFB[[curr_marker_gene]])
        colnames(tmp_df)[ncol(tmp_df)] = "value"
        tmp_df$symbol = curr_marker_gene
        if(is.null(embed_TFB2)){
            embed_TFB2 <<- tmp_df
        }else{
            embed_TFB2 <<- rbind(embed_TFB2, tmp_df)
        }
    })
    embed_TFB2$Clusters_TileMatrix_named = factor(embed_TFB2$Clusters_TileMatrix_named, levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))
    
    ps_tfb = list()  #umap
    ps_tfb2 = list() #violin
    local(for(tf in markersTFBs){
        p_TFB = NULL; p_TFB = ggplot(embed_TFB2[embed_TFB2$symbol == tf,], aes(x=UMAP_dim_1, y=UMAP_dim_2, color=value), fill="black") +
            geom_point(size=0.3, alpha=1) +
            scale_color_gradientn(colors=colorRampPalette(c("gray95", "red"))(100), limit=c(0,4), oob=scales::squish, labels=function(x) paste0(format(round(x, 2), nsmall = 2))) +
            ggtitle(paste0(tf)) +
            xlab(paste0("UMAP dim 1")) +
            ylab(paste0("UMAP dim 2")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
                  axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
                  legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        ps_tfb[[tf]] <<- p_TFB  
        
        embed_TFB2 = embed_TFB2[!embed_TFB2$Clusters_TileMatrix_named %in% c("Resistant_side", "Sensitive_side"),]
        embed_TFB2$Clusters_TileMatrix_named = factor(embed_TFB2$Clusters_TileMatrix_named, unique(embed_gex2$Clusters_TileMatrix_named)) 
        p_tfb2 = NULL; p_tfb2 = ggplot(embed_TFB2[embed_TFB2$symbol == tf,], aes(x=Clusters_TileMatrix_named, y=log2(value + 1), fill=Clusters_TileMatrix_named), fill="black") +
            geom_violin() +  geom_boxplot(outlier.shape=NA, width=0.1, color="grey70") +
            scale_fill_manual(values=clusters_named_vColors) +
            ggtitle(paste0(tf)) +
            ylab(paste0("chromVAR (z-score)")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            theme_minimal() + 
            stat_compare_means(comparisons=my_comparisons, method="wilcox.test", paired=F, size=3) +
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
                  axis.text=element_text(size=8, angle=90, hjust=1, vjust=0.5), strip.text.y=element_text(size=8, angle=0), legend.position="right",
                  legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        ps_tfb2[[tf]] <<- p_tfb2
    })
    saveRDS(file="../OUTPUT/chromvar_ordinations_YAP_FOSL_TEAD.rds", ps_tfb)
    saveRDS(file="../OUTPUT/chromvar_violin_YAP_FOSL_TEAD.rds", ps_tfb2)
    p_tfb_all = ggarrange(plotlist=ps_tfb)
    p_tfb_all2 = ggarrange(plotlist=ps_tfb2)
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_UMAP_TileMatrix_chromvar", ".pdf"), plot=p_tfb_all, device="pdf", height=5.48, width=10.12, units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_violin_TileMatrix_chromvar", ".pdf"), plot=p_tfb_all2, device="pdf", height=5.48, width=10.12, units="in", limitsize=FALSE)
    
    
    ##############################################################
    # Investigate with other annotation DBs for motifs           #
    #                                                            #
    ##############################################################
    embed_MM2 = NULL
    local(for(k in 1:length(markersMM)){
        curr_marker_gene = markersMM[k]
        tmp_df = embed_MM[,c("UMAP_dim_1","UMAP_dim_2", "Clusters_TileMatrix_named")]
        tmp_df = cbind(tmp_df, embed_MM[[curr_marker_gene]])
        colnames(tmp_df)[ncol(tmp_df)] = "value"
        tmp_df$symbol = curr_marker_gene
        if(is.null(embed_MM2)){
            embed_MM2 <<- tmp_df
        }else{
            embed_MM2 <<- rbind(embed_MM2, tmp_df)
        }
    })
    embed_MM2$Clusters_TileMatrix_named = factor(embed_MM2$Clusters_TileMatrix_named, 
                                                 levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))
    
    p_MM = NULL; p_MM = ggplot(embed_MM2, aes(x=UMAP_dim_1, y=UMAP_dim_2, color=value), fill="black") +
        facet_wrap(. ~ symbol) +
        geom_point(size=0.3, alpha=1) +
        scale_color_gradientn(colors=colorRampPalette(c("gray95", "red"))(100), limit=c(0,4), oob=scales::squish) +
        ggtitle(paste0("cisbp motif (chromVAR, z score)")) +
        xlab(paste0("UMAP dim 1")) +
        ylab(paste0("UMAP dim 2")) +
        guides(fill=guide_legend("")) +
        labs(color=NULL) +
        theme_minimal() + 
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
    p_MM
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_UMAP_TilleMatrix_chromvar_cisbp_motif", ".pdf"), plot=p_MM, device="pdf", height=7.85, width=12.4, units="in", limitsize=FALSE)
    
    p_MM2 = NULL; p_MM2 = ggplot(embed_MM2, aes(x=value, y=Clusters_TileMatrix_named), fill="black") +
        facet_grid(. ~ symbol, space="fixed", scales="free") +
        geom_density_ridges(alpha=0.5, fill="#00AFBB") +
        ggtitle(paste0("ENCODE_and_cistromeDB_TF_peaks (chromVAR, z-score)")) +
        ylab(paste0("Cluster")) +
        xlab(paste0("z-score")) +
        guides(fill=guide_legend("")) +
        labs(color=NULL) +
        theme_minimal() + 
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
    p_MM2
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_density_TilleMatrix_chromvar_encode_cistromeDB", ".pdf"), plot=p_MM2, device="pdf", height=3.7, width=15, units="in", limitsize=FALSE)
    
    embed_TF2 = NULL
    local(for(k in 1:length(markersTF)){
        curr_marker_gene = markersTF[k]
        tmp_df = embed_TF[,c("UMAP_dim_1","UMAP_dim_2", "Clusters_TileMatrix_named")]
        tmp_df = cbind(tmp_df, embed_TF[[curr_marker_gene]])
        colnames(tmp_df)[ncol(tmp_df)] = "value"
        tmp_df$symbol = curr_marker_gene
        if(is.null(embed_TF2)){
            embed_TF2 <<- tmp_df
        }else{
            embed_TF2 <<- rbind(embed_TF2, tmp_df)
        }
    })
    embed_TF2$Clusters_TileMatrix_named = factor(embed_TF2$Clusters_TileMatrix_named, 
                                                 levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))
    
    p_TF = NULL; p_TF = ggplot(embed_TF2, aes(x=UMAP_dim_1, y=UMAP_dim_2, color=value), fill="black") +
        facet_wrap(. ~ symbol) +
        geom_point(size=0.3, alpha=1) +
        scale_color_gradientn(colors=colorRampPalette(c("gray95", "red"))(100), limit=c(0,4), oob=scales::squish) +
        ggtitle(paste0("ENCODE_and_cistromeDB_TF_peaks (chromVAR, z score)")) +
        xlab(paste0("UMAP dim 1")) +
        ylab(paste0("UMAP dim 2")) +
        guides(fill=guide_legend("")) +
        labs(color=NULL) +
        theme_minimal() + 
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
    p_TF
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_UMAP_TilleMatrix_chromvar_encode_cistromeDB", ".pdf"), plot=p_TF, device="pdf", height=7.85, width=8.68, units="in", limitsize=FALSE)
    
    p_TF2 = NULL; p_TF2 = ggplot(embed_TF2, aes(x=value, y=Clusters_TileMatrix_named), fill="black") +
        facet_grid(. ~ symbol, space="fixed", scales="free") +
        geom_density_ridges(alpha=0.5, fill="#00AFBB") +
        ggtitle(paste0("ENCODE_and_cistromeDB_TF_peaks (chromVAR, z-score)")) +
        ylab(paste0("Cluster")) +
        xlab(paste0("z-score")) +
        guides(fill=guide_legend("")) +
        labs(color=NULL) +
        theme_minimal() + 
        theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
    p_TF2
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_density_TilleMatrix_chromvar_cisbp_motif", ".pdf"), plot=p_TF2, device="pdf", height=3.7, width=24, units="in", limitsize=FALSE)
    
    ############################################################################
    # Gene signatures                                                          #
    #                                                                          #
    ############################################################################
    # Load gene signatures of interest.
    se <- dsassembly::getExperiment("DS000000950")
    hippo = getFeatureSetCollection("GMTY188:analysis/hippo.gmt.bz2@REVISION-2")
    names(hippo) = hippo@elementMetadata@listData[["name"]]
    search.results <- searchFiles("hallmark")
    hallmark <- genomitory::getFeatureSetCollection(search.results[search.results$id=="GMTY42:human/H.gmt.bz2@REVISION-1",]$id)
    names(hallmark) <- hallmark@elementMetadata@listData[["name"]]
  
    E2F <- hallmark$HALLMARK_E2F_TARGETS
    E2F.symbol <- rowData(se)$symbol[match(E2F, rowData(se)$ID)]
    
    KRAS_up <- hallmark$HALLMARK_KRAS_SIGNALING_UP
    KRAS_up.symbol <- rowData(se)$symbol[match(KRAS_up, rowData(se)$ID)]
    
    KRAS_down <- hallmark$HALLMARK_KRAS_SIGNALING_DN
    KRAS_down.symbol <- rowData(se)$symbol[match(KRAS_down, rowData(se)$ID)]
    
    EMT <- hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
    EMT.symbol <- rowData(se)$symbol[match(EMT, rowData(se)$ID)]
    
    hippo[["EMT"]] = EMT.symbol
    hippo[["cluster4"]] = cluster4
    hippo[["E2F"]] = E2F.symbol
    hippo[["KRAS_up"]] = KRAS_up.symbol
    hippo[["KRAS_down"]] = KRAS_down.symbol
    
    my_comparisons = list(
        c("Sensitive_DMSO", "Sensitive_GNE7883"),
        c("Sensitive_DMSO", "Resistant_main"),
        c("Resistant_main", "Sensitive_GNE7883")
    )
    heights = list("hippo3"=11,
                   "hippo9"=11,
                   "hippo52"=7.5,
                   "hippo145"=6,
                   "mapk"=6.5,
                   "apoptosis"=12,
                   "proliferation"=6,
                   "hippo_mapk"=7.5,
                   "EMT"=4,
                   "cluster4"=14,
                   "E2F"=4,
                   "KRAS_up"=4,
                   "KRAS_down"=2) 
    ps = list()
    ps2 = list()
    local(for(k in 1:length(hippo)){
        curr_marker_genes = hippo[[k]]
        curr_name = names(hippo)[k]
        
        # GeneExpressionMatrix
        idx = which(rowData(seGEX)$name %in% curr_marker_genes)
        if(length(idx) > 2){
            seGEX2 = seGEX[idx,]
            embed_gex = getEmbedding(archr_proj,  embedding="UMAP_TileMatrix")
            selected_genes_cells = data.frame(as.matrix(t(assays(seGEX2)$GeneExpressionMatrix)))
            colnames(selected_genes_cells) = rowData(seGEX2)$name
            identical(row.names(embed_gex), row.names(selected_genes_cells))
            embed_gex = cbind(embed_gex, selected_genes_cells)
            embed_gex$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
            colnames(embed_gex)[1] = "UMAP_dim_1"
            colnames(embed_gex)[2] = "UMAP_dim_2"
            embed_gex$score = rowMeans(embed_gex[,c(3:(ncol(embed_gex)-1))])
            embed_gex$Clusters_TileMatrix_named = factor(embed_gex$Clusters_TileMatrix_named, levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))
            
            p_gex = NULL; p_gex = ggplot(embed_gex, aes(x=UMAP_dim_1, y=UMAP_dim_2, color=score), fill="black") +
                geom_point(size=0.3, alpha=1) +
                scale_color_gradientn(colors=colorRampPalette(c("gray", "darkred"))(100), labels=function(x) paste0(format(round(x, 2), nsmall = 2))) +
                ggtitle(paste0(curr_name)) +
                xlab(paste0("UMAP dim 1")) +
                ylab(paste0("UMAP dim 2")) +
                guides(fill=guide_legend("")) +
                labs(color=NULL) +
                theme_minimal() + 
                theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
                      axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
                      legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
            ps[[curr_name]] <<- p_gex
            
            embed_gex2 = embed_gex[!embed_gex$Clusters_TileMatrix_named %in% c("Resistant_side", "Sensitive_side"),]
            embed_gex2$Clusters_TileMatrix_named = factor(embed_gex2$Clusters_TileMatrix_named, unique(embed_gex2$Clusters_TileMatrix_named))
            p_gex2 = NULL; p_gex2 = ggplot(embed_gex2, aes(x=Clusters_TileMatrix_named, y=score, fill=Clusters_TileMatrix_named), fill="black") +
                geom_violin() + geom_boxplot(outlier.shape=NA, width=0.1, color="grey70") +
                scale_fill_manual(values=clusters_named_vColors) +
                ggtitle(paste0(curr_name)) +
                ylim(c(0,heights[[curr_name]])) +
                ylab(paste0("score")) +
                guides(fill=guide_legend("")) +
                labs(color=NULL) +
                stat_compare_means(comparisons=my_comparisons, method="wilcox.test", paired=F, size=3) +
                stat_compare_means(label.y = 45) +
                theme_minimal() + 
                theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5),
                      axis.title.y=element_text(size=8), 
                      axis.title.x=element_blank(),
                      axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5), 
                      axis.text.y=element_text(size=10, angle=0, hjust=1, vjust=0.5), 
                      strip.text.y=element_text(size=8, angle=0), 
                      legend.position="right",
                      legend.text=element_text(size=8, angle=0), 
                      plot.title=element_text(size=10), 
                      legend.key.size=unit(0.35, 'cm')
                )
            ps2[[curr_name]] <<- p_gex2
        }
    })
    p_signatures = ggarrange(plotlist=ps, ncol=4, nrow=3)
    ggsave(paste0(curr_outdir, "/gene_signatures_TileMatrix", ".pdf"), plot=p_signatures, device="pdf", height=6.8, width=9.8, units="in", limitsize=FALSE)
    p_signatures2 = ggarrange(plotlist=ps2, ncol=4, nrow=3)
    ggsave(paste0(curr_outdir, "/gene_signatures_TileMatrix_for_fig_3C_3H", ".pdf"), plot=p_signatures2, device="pdf", height=10.8, width=17.4, units="in", limitsize=FALSE)
    saveRDS(file="../OUTPUT/hippo_signatures_gex_scores_ordinations.rds", ps)
    saveRDS(file="../OUTPUT/hippo_signatures_gex_scores_violin.rds", ps2)

    message("Completed archr_marker_analysis_hippo.R!")
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
archr_marker_genes(indir=indir, num_threads=num_threads, my_FDR, my_log2FC)
