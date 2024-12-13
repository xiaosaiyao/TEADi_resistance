#!/usr/bin/env Rscript

# Workflow to execute epiregulon
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

do_epiregulon <- function(indir=NULL, num_threads=num_threads, my_FDR=0.1, my_log2FC=1){ #, gene_signature="tead"){ #}, my_FDR=0.1, my_log2FC=0.5) {

    #setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    #indir = "/gstore/project/tead/scMultiome/NGS4467_pilot_JT/OUTPUT/archr/"
    #num_threads = 4
    #my_FDR = 0.1
    #my_log2FC = 1
    #one_panel_x = 2.11
    #one_panel_y = 2.35
    
    my_FDR = as.numeric(my_FDR)
    my_log2FC = as.numeric(my_log2FC)
    num_threads = as.numeric(num_threads)
    message("Using num_threads=", num_threads)#, ". log2FC=", my_log2FC, ". FDR=", my_FDR)
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
    library(epiregulon)
    library(epiregulon.archr)
    library(epiregulon.extra)
    library(scater)
    source("./utils.R")
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(msigdbr)
    library(RColorBrewer)
    library(igraph)
    library(ggrepel)
    
    archr_proj = ArchR::loadArchRProject(indir)
    
    curr_outdir = "../OUTPUT/regulon/"
    dir.create(curr_outdir)
    
    # Loading metadata
    mapping = getCellColData(archr_proj, select = c("hash_assignment2", "CNAME", "TEST_ARTICLE", "Treatment", "library", "HTO", "Treatment_HTO", "Treatment_library", "Clusters_TileMatrix_named"))
    mapping$Treatment_library_HTO = paste0(mapping$Treatment_library, "_", mapping$HTO)
    archr_proj$Treatment_library_HTO = mapping$Treatment_library_HTO
    mapping = unique(mapping)
    mapping = data.frame(mapping)
    
    # Define list of colors to be used for Treatment_HTO
    curr_vColors = vColors[1:length(unique(mapping$Treatment_HTO))]
    names(curr_vColors) = unique(mapping$Treatment_HTO)
    
    ###########################################################
    # Start epiregulon                                        #
    #                                                         #
    ###########################################################
    gene_expression_matrix = getMatrixFromProject(ArchRProj=archr_proj ,useMatrix="GeneExpressionMatrix", useSeqnames=NULL, verbose=TRUE, binarize=FALSE, threads=getArchRThreads())
    gene_expression_matrix = ArchRMatrix2SCE(gene_expression_matrix, rename = "normalizedCounts")
    rownames(gene_expression_matrix) = rowData(gene_expression_matrix)$name
    reducedDim(gene_expression_matrix, "UMAP_TileMatrix") = getEmbedding(ArchRProj=archr_proj, embedding="UMAP_TileMatrix", returnDF=TRUE)[colnames(gene_expression_matrix),]
    plotReducedDim(gene_expression_matrix, dimred = "UMAP_TileMatrix", text_by = "Clusters_TileMatrix_named", colour_by = "Clusters_TileMatrix_named")
    reducedDim(gene_expression_matrix, "IterativeLSI_TileMatrix") = getReducedDims(ArchRProj=archr_proj, reducedDims="IterativeLSI_TileMatrix")
    
    # load other matrices
    peak_matrix = getMatrixFromProject(ArchRProj=archr_proj, useMatrix="PeakMatrix", useSeqnames=NULL, verbose=TRUE, binarize=FALSE, threads=1)
    peak_matrix = as(peak_matrix, "SingleCellExperiment")
    peak_matrix = peak_matrix[, colnames(gene_expression_matrix)]
    names(assays(peak_matrix)) = "counts"
    
    TF_bindingMatrix <- getMatrixFromProject(ArchRProj=archr_proj, useMatrix="TF_bindingMatrix", useSeqnames=NULL, verbose=TRUE, binarize=FALSE, threads=1)
    
    pdf("../OUTPUT/regulon/cluster_plots.pdf")
    plotReducedDim(gene_expression_matrix, dimred = "UMAP_TileMatrix", text_by = "Clusters_TileMatrix_named", colour_by = "Clusters_TileMatrix_named", rasterise = TRUE)
    plotReducedDim(gene_expression_matrix, dimred = "UMAP_TileMatrix", text_by = "Clusters_TileMatrix_named", colour_by = "Gex_MitoRatio", rasterise = TRUE)
    plotReducedDim(gene_expression_matrix, dimred = "UMAP_TileMatrix", text_by = "Clusters_TileMatrix_named", colour_by = "Gex_RiboRatio", rasterise = TRUE)
    plotReducedDim(gene_expression_matrix, dimred = "UMAP_TileMatrix", text_by = "Clusters_TileMatrix_named", colour_by = "DoubletScore", rasterise = TRUE)
    plotReducedDim(gene_expression_matrix, dimred = "UMAP_TileMatrix", text_by = "Clusters_TileMatrix_named", colour_by = "TSSEnrichment", rasterise = TRUE)
    dev.off()
    
    # load H226 peaks
    TEAD1 <- read.delim("/gstore/data/genomics/congee_rest_runs/649086a939473e4888ac65ca/SAM24431959/croo_output/LIB5470293_SAM24431959_macs_peaks.narrowPeak", header = FALSE)
    colnames(TEAD1) <- c("chr","start","end")
    TEAD1 <- makeGRangesFromDataFrame(TEAD1)
    
    # Retrieve bulk TF ChIP-seq binding sites
    grl = getTFMotifInfo(genome = "hg38")
    grl$TEAD1 <- TEAD1
    
    set.seed(1010)
    p2g = calculateP2G(peakMatrix = peak_matrix,
                        peak_assay = "counts",
                        expMatrix = gene_expression_matrix,
                        exp_assay = "normalizedCounts",
                        reducedDim = reducedDim(gene_expression_matrix, "IterativeLSI_TileMatrix"))
    head(p2g)
    
    # Construct Regulon
    overlap <- addTFMotifInfo(archR_project_path=indir, grl=grl, p2g=p2g)
    regulon_df <- getRegulon(p2g, overlap)
    
    # prune regulon
    pruned.regulon <- pruneRegulon(expMatrix = gene_expression_matrix,
                                   exp_assay = "normalizedCounts",
                                   peakMatrix = peak_matrix,
                                   peak_assay = "counts",
                                   test = "chi.sq",
                                   regulon_df,
                                   clusters = gene_expression_matrix$Clusters_TileMatrix_named,
                                   prune_value = "pval",
                                   regulon_cutoff = 0.05
    )
    
    # Add weights
    regulon.w.wilcox <- addWeights(regulon=pruned.regulon,
                                   expMatrix=gene_expression_matrix,
                                   peakMatrix = peak_matrix,
                                   exp_assay="normalizedCounts",
                                   peak_assay="counts",
                                   method="wilcox",
                                   exp_cutoff=0.5, min_targets=3,
                                   clusters=gene_expression_matrix$Clusters_TileMatrix_named
    )
    regulon.w.wilcox
    
    # add motif
    regulon.w.wilcox.motif <- addMotifScore(regulon=regulon.w.corr, archr_path=indir)
    regulon.w.wilcox.motif$weight.motif <- regulon.w.wilcox.motif$weight*regulon.w.wilcox.motif$motif
    
    #############################
    # calculate activity
    # use TF gene expression
    score.combine.wilcox <- calculateActivity(expMatrix = gene_expression_matrix,
                                              regulon = regulon.w.wilcox,
                                              mode = "weight",
                                              method ="weightedMean",
                                              exp_assay="normalizedCounts")
    
    epiregulon_output = list(
        "pruned.regulon" = pruned.regulon,
        "regulon.w.wilcox" = regulon.w.wilcox,
        "score.combine.wilcox"=score.combine.wilcox
    )
    saveRDS(file="../OUTPUT/regulon/epiregulon.rds", epiregulon_output)
    
    # checkpoint
    #epiregulon_output = readRDS(file="../OUTPUT/regulon/epiregulon.rds")
    #pruned.regulon = epiregulon_output[["pruned.regulon"]]
    #regulon.w.wilcox = epiregulon_output[["regulon.w.wilcox"]]
    #score.combine.wilcox = epiregulon_output[["score.combine.wilcox"]]
    
    #########################################
    # hippo signatures                      #
    #                                       #
    #########################################
    library(genomitory)
    hippo <- getFeatureSetCollection("GMTY188:analysis/hippo.gmt.bz2@REVISION-2")
    names(hippo) <- hippo@elementMetadata@listData[["name"]]
    
    # EMT
    se <- dsassembly::getExperiment("DS000000950")
    search.results <- searchFiles("hallmark")
    hallmark <- genomitory::getFeatureSetCollection(search.results[3,"id"])
    names(hallmark) <- hallmark@elementMetadata@listData[["name"]]
    EMT <- hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
    EMT.symbol <- rowData(se)$symbol[match(EMT, rowData(se)$ID)]
    hippo[["EMT"]] = EMT.symbol
    hippo[["cluster4"]] = cluster4
    
    regulon.hippo.signature <- epiregulon:::genesets2regulon(hippo)
    score.combine.hippo <- calculateActivity(expMatrix = gene_expression_matrix,
                                             regulon = regulon.hippo.signature,
                                             mode = "weight",
                                             method ="weightedMean",
                                             exp_assay="normalizedCounts")
    
    tfs <- c("YAP1","WWTR1","TEAD1","TEAD4","FOSL1","SATB2")
    
    pdf("../OUTPUT/regulon/violin.plots.pdf")
    plotActivityViolin(score.combine.wilcox,
                       tf=tfs,
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="wilcox",
                       nrow=2,
                       boxplot=TRUE)
    
    plotActivityViolin(score.combine.hippo,
                       tf=c("hippo52","hippo145","mapk","proliferation"),
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="hippo signature",
                       boxplot=TRUE,
                       nrow=2,
                       ncol=2)
    
    dev.off()
    
    # differential activity using weights
    markers <- findDifferentialActivity(score.combine.wilcox,
                                        gene_expression_matrix$Clusters_TileMatrix_named,
                                        pval.type="some",
                                        direction="up",
                                        test.type= "t")
    markers.sig <- getSigGenes(markers, topgenes=5 )
    
    
    pdf("../OUTPUT/regulon/bubble.plots.pdf")
    plotBubble(activity_matrix = score.combine.wilcox,
               tf = markers.sig$tf,
               clusters = gene_expression_matrix$Clusters_TileMatrix_named)
    dev.off()
    
    #umap
    pdf("../OUTPUT/regulon/umap.plots.pdf", width=12)
    plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=score.combine.wilcox,
                    tf= tfs, dimtype="UMAP_TileMatrix",
                    label = "Clusters_TileMatrix_named",
                    point_size=0.1,
                    rasterise = TRUE,
                    title = "wilcox",
                    nrow=2)
    
    dev.off()
    
    # Plot gene expresssion
    gex = assay(gene_expression_matrix, "normalizedCounts")
    pdf("../OUTPUT/regulon/geneexpression.pdf", width=12)
    plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=gex,
                    tf=tfs,
                    dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    legend.label="Gene Exp",
                    point_size=0.1,
                    color=c("grey","blue"),
                    limit = c(0,2),
                    rasterise=TRUE)
    plotActivityViolin(gex,
                       tf = tfs,
                       clusters =  gene_expression_matrix$Clusters_TileMatrix_named,
                       legend.label = "gene expression")
    dev.off()
    
    # Plot chromVar
    chromvar <- assay(tf_binding_matrix, "z")
    chromvar <- chromvar[tfs,]
    
    pdf("../OUTPUT/regulon/chromvar.pdf", width=12)
    plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=chromvar,
                    tf= tfs,
                    dimtype="UMAP_TileMatrix",
                    label = "Clusters_TileMatrix_named",
                    legend.label = "chromvar",
                    point_size=0.1,
                    color = c("grey","red"),
                    limit = c(0,3),
                    rasterise = TRUE)
    plotActivityViolin(chromvar,
                       tf = tfs,
                       clusters = gene_expression_matrix$Clusters_TileMatrix_named,
                       legend.label = "chromvar",
                       nrow=2,
                       boxplot = TRUE)
    dev.off()
    
    # Plot hippo score
    pdf("../OUTPUT/regulon/hippo.pdf", width=12)
    plotActivityDim(sce = gene_expression_matrix,
                    activity_matrix = score.combine.hippo,
                    tf = rownames( score.combine.hippo),
                    dimtype="UMAP_TileMatrix",
                    label = "Clusters_TileMatrix_named",
                    legend.label = "score",
                    point_size=0.1,
                    color = c("grey","brown"),
                    limit = c(1,3),
                    rasterise = TRUE)
    plotActivityViolin(score.combine.hippo,
                       tf = rownames(score.combine.hippo),
                       clusters = gene_expression_matrix$Clusters_TileMatrix_named,
                       legend.label = "score",
                       nrow=2,
                       boxplot = TRUE)
    dev.off()
    
    #######################################
    #  GeneSet enrichment
    H <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb", cat = "H", gene.id.type = "SYMBOL" )
    C2 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb", cat = "C2", gene.id.type = "SYMBOL" )
    C6 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb", cat = "C6", gene.id.type = "SYMBOL" )
    C8 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb", cat = "C8", gene.id.type = "SYMBOL" )
    
    # only H and C6
    gs <- c(H, C6)
    gs.list <- do.call(rbind,lapply(names(gs), function(x) {data.frame(gs=x, genes=gs[[x]])}))
    
    enrichresults <- regulonEnrich(tfs,
                                   regulon=regulon.w.wilcox,
                                   weight ="weight",
                                   weight_cutoff=0,
                                   genesets=gs.list)
    
    # Write table
    gsea_df = NULL
    local(for(name in names(enrichresults)){
        print("name:")
        print(name)
        curr_df = enrichresults[[name]]
        print(head(curr_df))
        curr_df$tf = name
        curr_df = na.omit(curr_df)
    
        
        if(is.null(gsea_df)){
            gsea_df <<- curr_df
        }else{
            gsea_df <<- rbind(gsea_df, curr_df)
        }
    })
    write.table(gsea_df, "../OUTPUT/regulon/GSEA_table_epiregulon.tsv", row.names=F, sep="\t", quote=F)
    
    #plot results
    p1 = enrichPlot(results=enrichresults)
    p2 = plotGseaNetwork(tf=c("TEAD1", "FOSL1", "SATB2"), enrichresults, p.adj_cutoff = 0.1, ntop_pathways = 10)
    p1_p2 = ggarrange(p1, p2, ncol=1, nrow=2)
    ggsave("../OUTPUT/regulon/GSEA_TFs_regulon.pdf", plot=p1_p2, device="pdf", height=12, width=25, units="in", limitsize=FALSE)
    
    # Do with more categories
    gs <- c(H, C2, C6, C8)
    gs.list <- do.call(rbind,lapply(names(gs), function(x) {data.frame(gs=x, genes=gs[[x]])}))
    
    enrichresults <- regulonEnrich(tfs,
                                   regulon=regulon.w.wilcox,
                                   weight ="weight",
                                   weight_cutoff=0,
                                   genesets=gs.list)
    
    #plot results
    pdf("../OUTPUT/regulon/GSEA_TFs_regulon_H_C2_C6_C8.pdf", width = 28, height = 10)
    enrichPlot(results=enrichresults)
    plotGseaNetwork(tf=names(enrichresults), enrichresults, p.adj_cutoff = 0.1, ntop_pathways = 10)
    dev.off()
    
    # create graphs
    Resistance_main_network <- buildGraph(regulon.w.wilcox, weights = "weight", cluster="Resistant_main")
    Sensitive_GNE7883_network <- buildGraph(regulon.w.wilcox, weights = "weight", cluster="Sensitive_GNE7883")
    
    # construct a difference graph
    diff_graph_title <- "differential TFs (Resistance_main - Sensitive_GNE7883) ranked by degree centrality"
    diff_graph <- buildDiffGraph(Resistance_main_network, Sensitive_GNE7883_network, abs_diff = FALSE)
    diff_graph <- addCentrality(diff_graph)
    diff_graph <- normalizeCentrality(diff_graph)
    rank_table <- rankTfs(diff_graph)
    
    pdf("../OUTPUT/regulon/ranktable.pdf")
    ggplot(rank_table, aes(x = rank, y = centrality)) +
        geom_point() +
        ggrepel::geom_text_repel(data = rbind(head(rank_table,10),
                                              tail(rank_table,10)),
                                 max.overlaps = Inf,
                                 aes(label = tf),
                                 nudge_x = 0, nudge_y = 0, box.padding = 0.5) +
        theme_classic() + ggtitle(diff_graph_title)
    dev.off()
    
    pdf("../OUTPUT/regulon/ranktable_v2.pdf", width=6, height=6)
    ggplot(rank_table, aes(x = rank, y = centrality)) +
        geom_point(size=1) +
        ggrepel::geom_text_repel(data = rbind(head(rank_table,63),
                                              tail(rank_table,63)),
                                 max.overlaps=30,
                                 size=2,
                                 aes(label=tf),
                                 nudge_x=0.5, nudge_y=0.5, box.padding=0.5,
                                 label.r=0.005, label.padding=0.09,) +
        theme_classic() + 
        ggtitle(diff_graph_title)
    dev.off()
    
    #################################################
    # Then try resistant main vs sensitive dmso.    #
    #################################################
    Resistance_main_network <- buildGraph(regulon.w.wilcox, weights = "weight", cluster="Resistant_main")
    Sensitive_DMSO_network <- buildGraph(regulon.w.wilcox, weights = "weight", cluster="Sensitive_DMSO")
    
    # construct a difference graph
    diff_graph_title <- "differential TFs (Resistance_main - Sensitive_DMSO) ranked by degree centrality"
    diff_graph <- buildDiffGraph(Resistance_main_network, Sensitive_DMSO_network, abs_diff = FALSE)
    diff_graph <- addCentrality(diff_graph)
    diff_graph <- normalizeCentrality(diff_graph)
    rank_table <- rankTfs(diff_graph)
    
    pdf("../OUTPUT/regulon/ranktable_resistant_main_vs_sensitive_DMSO.pdf")
    ggplot(rank_table, aes(x = rank, y = centrality)) +
        geom_point(size=1) +
        ggrepel::geom_text_repel(data = rbind(head(rank_table,43),
                                              tail(rank_table,43)),
                                 max.overlaps=30,
                                 size=2,
                                 aes(label=tf),
                                 nudge_x=0.5, nudge_y=0.5, box.padding=0.5,
                                 label.r=0.005, label.padding=0.09,) +
        theme_classic() + 
        ggtitle(diff_graph_title)
    dev.off()
    
    # compute similarity
    TFs_interest = c("TEAD1", "FOSL1", "SATB2")
    for(TF_interest in TFs_interest){
        diff_graph_filter <- igraph::subgraph.edges(diff_graph,
                                                    E(diff_graph)[E(diff_graph)$weight>0],
                                                    del=TRUE)
        
        
        # compute a similarity matrix of all TFs
        similarity_score <- calculateJaccardSimilarity(diff_graph_filter)
        
        # Focus on TF of interest
        similarity_score_TF_interest <- similarity_score[, TF_interest]
        similarity_df <- data.frame(similarity = head(sort(similarity_score_TF_interest,
                                                           decreasing = TRUE),20),
                                    TF = names(head(sort(similarity_score_TF_interest,
                                                         decreasing = TRUE),20)))
        
        similarity_df$TF <- factor(similarity_df$TF,
                                   levels = rev(unique(similarity_df$TF)))
        
        # plot top TFs most similar to TF of interest
        topTFplot <- ggplot(similarity_df, aes(x=TF, y=similarity)) +
            geom_bar(stat="identity") +
            coord_flip() +
            ggtitle(paste(TF_interest, "similarity")) +
            theme_classic()
        
        print(topTFplot)
        
        # account for background
        permute_matrix <- permuteGraph(diff_graph_filter, TF_interest, 100, p=1)
        permute_matrix <- permute_matrix[names(similarity_score_TF_interest),]
        diff_matrix <- similarity_score_TF_interest-rowMeans(permute_matrix)
        
        diff_matrix_df <- data.frame(similarity = head(sort(diff_matrix,
                                                            decreasing = TRUE),20),
                                     TF = names(head(sort(diff_matrix,
                                                          decreasing = TRUE),20)))
        
        diff_matrix_df$TF <- factor(diff_matrix_df$TF, levels = rev(unique(diff_matrix_df$TF)))
        
        # plot top TFs most similar to EBF1
        topTFplot <- ggplot(diff_matrix_df, aes(x=TF, y=similarity)) +
            geom_bar(stat="identity") +
            coord_flip() +
            ggtitle(paste(TF_interest, "similarity")) +
            theme_classic()
        
        pdf(paste0("../OUTPUT/regulon/similarity.background.subtracted_", TF_interest, ".pdf"), width=2, height=5)
        print(topTFplot)
        dev.off()
    }
    
    ####################################################
    # Generate final plots:                            #
    ####################################################
    # TF
    markersTF = c("TEAD1", "YAP1", "FOSL1", "SATB2")
    idx = which(row.names(score.combine.wilcox) %in% markersTF)
    score.combine.wilcox_select = score.combine.wilcox[idx,]
    embed_TFAct = getEmbedding(archr_proj,  embedding="UMAP_TileMatrix")
    dim(score.combine.wilcox_select)
    dim(embed_TF)
    score.combine.wilcox_select_t = data.frame(as.matrix(t(score.combine.wilcox_select)))
    identical(row.names(embed_TFAct), row.names(score.combine.wilcox_select_t))
    embed_TFAct = cbind(embed_TFAct, score.combine.wilcox_select_t)
    embed_TFAct$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
    colnames(embed_TFAct)[1] = "UMAP_dim_1"
    colnames(embed_TFAct)[2] = "UMAP_dim_2"
    
    embed_TFAct2 = NULL
    local(for(k in 1:length(markersTF)){
        curr_marker_gene = markersTF[k]
        tmp_df = embed_TF[,c("UMAP_dim_1","UMAP_dim_2", "Clusters_TileMatrix_named")]
        tmp_df = cbind(tmp_df, embed_TFAct[[curr_marker_gene]])
        colnames(tmp_df)[ncol(tmp_df)] = "value"
        tmp_df$symbol = curr_marker_gene
        if(is.null(embed_TFAct2)){
            embed_TFAct2 <<- tmp_df
        }else{
            embed_TFAct2 <<- rbind(embed_TFAct2, tmp_df)
        }
    })
    embed_TFAct2$Clusters_TileMatrix_named = factor(embed_TFAct2$Clusters_TileMatrix_named, levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))
    
    my_labels_df3 = data.frame(Clusters_TileMatrix_named=unique(embed_TFAct2$Clusters_TileMatrix_named),
                               label=paste0("", unique(embed_TFAct2$Clusters_TileMatrix_named)))
    my_labels3 = embed_TFAct2 %>%
        dplyr::group_by(Clusters_TileMatrix_named) %>%
        dplyr::summarize(UMAP_dim_1 = mean(UMAP_dim_1), UMAP_dim_2 = mean(UMAP_dim_2)) %>%
        dplyr::left_join(my_labels_df3) %>% as.data.frame()
    
    clusters_named_vColors = c("Resistant_main" = "#000000", 
                               "Resistant_side" ="#696969", 
                               "Sensitive_DMSO" = "#8b4513", 
                               "Sensitive_GNE7883" = "#006400",
                               "Sensitive_side" = "#808000")
    
    embed_TFAct2$Clusters_TileMatrix_named = factor(embed_TFAct2$Clusters_TileMatrix_named, levels=c("Sensitive_DMSO", "Sensitive_GNE7883", "Sensitive_side", "Resistant_main", "Resistant_side"))
    
    p_tf_list = list()
    p_tf_list2 = list()
    local(for(gene in markersTF){
        p_TF = NULL; p_TF = ggplot(embed_TFAct2[embed_TFAct2$symbol == gene,], aes(x=UMAP_dim_1, y=UMAP_dim_2, color=value), fill="black") +
            geom_point(size=0.2, alpha=1) +
            scale_color_gradientn(colors=c("#352A86", "#343DAE", "orange", "#F6DA23", "#F8FA0D" )) +
            ggtitle(paste0(gene)) +
            xlab(paste0("UMAP dim 1")) +
            ylab(paste0("UMAP dim 2")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            geom_label_repel(data=my_labels3, alpha=0.8, label.size=0.05, fill="gray90", color="black", size=3, fontface='plain', 
                             box.padding=0.01, aes(label=my_labels3[["label"]])) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        #p_TF
        p_tf_list[[gene]] <<- p_TF
        
        p_TF2 = NULL; p_TF2 = ggplot(embed_TFAct2[embed_TFAct2$symbol == gene,], aes(x=Clusters_TileMatrix_named, y=value, fill=Clusters_TileMatrix_named), fill="black") +
            geom_violin() + geom_boxplot(outlier.shape=NA, width=0.1, color="grey70") + #geom_point(alpha=0.5, size=0.05, position=position_jitter(0.1)) +
            scale_fill_manual(values=clusters_named_vColors) +
            ggtitle(paste0(gene)) +
            xlab(paste0("")) +
            ylab(paste0("")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), axis.text.x=element_text(size=8, angle=90, vjust=0.5, hjust=1), 
                  axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right", 
                  legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        p_tf_list2[[gene]] <<- p_TF2
    })
    p_tfs = ggarrange(plotlist=p_tf_list, ncol=4)
    p_tfs2 = ggarrange(plotlist=p_tf_list2, ncol=1, nrow=4)
    saveRDS(file="../OUTPUT/epiregulon_activity_ordinations.rds", p_tf_list)
    saveRDS(file="../OUTPUT/epiregulon_activity_violin.rds", p_tf_list2)
    ggsave(paste0(curr_outdir, "/FOSL_TEAD_YAP_SATB2_UMAP_TileMatrix_epiregulon_activity", ".pdf"), plot=p_tfs, device="pdf", height=4, width=15, units="in", limitsize=FALSE)

    ###########################################
    # Generate GSEA enrichment plots.
    #
    msigdb_results_df = NULL
    local(for(symbol in names(enrichresults)){
        curr_df = enrichresults[[symbol]]
        curr_df$geneID = NULL
        curr_df$symbol = symbol
        idxs = which(!is.na(curr_df$p.adjust))
        curr_df = curr_df[idxs,]
        if(nrow(curr_df) != 0){
            if(is.null(msigdb_results_df)){
                msigdb_results_df <<- curr_df
            }else{
                msigdb_results_df <<- rbind(msigdb_results_df, curr_df)
            }
        }
    })
    
    # order by NES
    curr_msigdb_results_df_all = NULL
    local(for(curr_category in c("H", "C2", "C5", "C6", "C8")){
        
        pathways_in_category = EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
                                          cat = curr_category, gene.id.type = "SYMBOL" )
        pathways_in_category = names(pathways_in_category)
        
        curr_msigdb_results_df = msigdb_results_df[msigdb_results_df$ID %in% pathways_in_category,]
        order_df = curr_msigdb_results_df %>%
            dplyr::group_by(ID) %>%
            dplyr::summarise_at(vars(p.adjust), list(SUM = ~ sum(.),
                                                SD = ~sd(.))) %>%
            as.data.frame()
        order_df = order_df[order(-order_df$SUM),]
        
        curr_msigdb_results_df$ID = factor(curr_msigdb_results_df$ID, levels=unique(order_df$ID))
        IDs_vector = unique(curr_msigdb_results_df$ID)
        selection = as.character(IDs_vector[grep("*", IDs_vector, ignore.case=TRUE)])
       
        # #Store in another df for further use downstream.
        curr_msigdb_results_df = curr_msigdb_results_df[curr_msigdb_results_df$ID %in% selection,]
        if(nrow(curr_msigdb_results_df != 0)){
            curr_msigdb_results_df$category = curr_category
            if(is.null(curr_msigdb_results_df_all)){
                curr_msigdb_results_df_all <<- curr_msigdb_results_df
            }else{
                curr_msigdb_results_df_all <<- rbind(curr_msigdb_results_df_all, curr_msigdb_results_df)
            }
            
            curr_p = ggplot(curr_msigdb_results_df, aes(x=(-1*log10(p.adjust)), y=ID, fill=GeneRatio)) +
                geom_bar(stat="identity") +
                facet_grid(. ~ symbol, scales="free_y", space="free_y") +
                geom_vline(xintercept=0, linetype="dashed",  color="red", linewidth=0.4) +
                ylab("") +
                ggtitle(paste0("GSEA analysis; selected functions. Category:", curr_category,"\nNES:Normalized Enrichment Score")) +
                theme_minimal() + theme(
                    panel.border=element_rect(linewidth=0.5),
                    panel.background=element_rect(fill='transparent', color=NA),
                    rect = element_rect(fill = "transparent"),
                    axis.text.x=element_text(angle=90),
                    axis.text.y=element_text(angle=0, size=7),
                    strip.text.y=element_text(angle=0, hjust=0),
                    strip.text.x=element_text(angle=0, hjust=0.5),
                    plot.title=element_text(size=8)
                )
        }
    })
    
    IDs_vector = curr_msigdb_results_df_all$ID
    curr_p = ggplot(curr_msigdb_results_df_all, aes(x=log10(p.adjust), y=ID, fill=GeneRatio)) + # do -1*
        geom_bar(stat="identity") +
        facet_grid(category ~ symbol, scales="free_y", space="free_y") +
        geom_vline(xintercept=0, linetype="dashed",  color="red", linewidth=0.4) +
        ylab("") +
        ggtitle(paste0("GSEA analysis; selected functions. Category:", "H,C2,C6,C8" ,"\nNES:Normalized Enrichment Score")) +
        theme_minimal() + theme(
            panel.border=element_rect(linewidth=0.5),
            panel.background=element_rect(fill='transparent', color=NA),
            rect = element_rect(fill = "transparent"),
            axis.text.x=element_text(angle=90),
            axis.text.y=element_text(angle=0, size=7),
            strip.text.y=element_text(angle=0, hjust=0),
            strip.text.x=element_text(angle=0, hjust=0.5),
            plot.title=element_text(size=8)
        )
    #curr_p
    curr_height = generate_plot_dimensions(curr_msigdb_results_df_all, by_variable="category")[[1]]
    curr_width = generate_plot_dimensions(curr_msigdb_results_df_all, by_variable="category")[[2]]
    ggsave(paste0("../OUTPUT/regulon/GSEA_epiregulon_barplot.pdf"), curr_p, width=curr_width+4, height=curr_height, units="in", device="pdf", limitsize=FALSE)
    
    message("Completed epiregulon.R!")
}

usage=function(errM) {
    cat("\nUsage : Rscript do_epiregulon.R [option] <Value>\n")
    cat("       -i        : indir containing the archR data structure.")
    cat("       -n        : num_threads.")
    cat("       -s        : FDR cutoff. Default=0.1")
    cat("       -f        : log2FC. Default=1")
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
    }else if(ARG[i] == "-c") {
        comparisons_file=ARG[i+1]
    }
}
do_epiregulon(indir=indir, num_threads=num_threads, comparisons_file=comparisons_file, my_FDR=my_FDR, my_log2FC=my_log2FC)
