#!/usr/bin/env Rscript

# Workflow to generate figures from epiregulon results.
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

do_epiregulon_plots <- function(indir=NULL, tfs=tfs, use_groups=use_groups, background_groups=background_groups){ #, gene_signature="tead"){ #}, my_FDR=0.1, my_log2FC=0.5) {

    #setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    #indir = "/gstore/project/tead/scMultiome/NGS4467_pilot_JT/OUTPUT/archr/"
    #tfs = "TEAD1, FOSL1"
    #use_groups        =  "Resistant_main, Resistant_main,    Sensitive_GNE7883"
    #background_groups =  "Sensitive_DMSO, Sensitive_GNE7883, Sensitive_DMSO"
   
    # Define variables for dimensions when generating plots.
    one_panel_x = 4.11
    one_panel_y = 4.35
    my_df = data.frame(ratio=c(1.05, 1.55), nrow=c(2,3))
    my_lm = lm(ratio ~ nrow, data=my_df)
    
    if(is.null(indir)){ # TODO check if dir
        stop("an indir has to be included")
    }
    if(is.null(tfs)){
        stop("A vector of marker TFs have to be specified.")
    }else{
        tfs = unlist(strsplit(tfs, split=","))
        tfs = gsub(" ", "", tfs)
    }
    if(is.null(use_groups)){
        stop("Please specify at least one use group. You can also specify multiple use groups seperated by a ','. -u Group1,Group2 and background groups: -b Group3,Group4.")
    }
    if(is.null(background_groups)){
        stop("Please specify at least one background_groups group. You can also specify multiple background_groups groups seperated by a ','. -u Group1,Group2 and use groups: -b Group3,Group4.")
    }
    use_groups = unlist(strsplit(use_groups, split=",")) 
    use_groups = gsub(" ", "", use_groups)
    background_groups = unlist(strsplit(background_groups, split=",")) 
    background_groups = gsub(" ", "", background_groups)
    
    if(length(background_groups) != length(use_groups)){
        stop("background groups have to be the same length as use groups")  
    }
        
    library(ArchR)
    library(parallel)
    library(BSgenome.Hsapiens.Genentech.GRCh38)
    library(gridExtra)
    library(genomitory)
    library(ggpubr)
    library(grid)
    library(epiregulon)
    library(epiregulon.extra)
    library(epiregulon.archr)
    library(scater)
    source("./utils.R")
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(msigdbr)
    library(RColorBrewer)
    library(igraph)
    library(ggrepel)
    
    ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
    
    generate_plot_dimensions <- function(curr_msigdb_results_df, WIDTH=7.9, HEIGHT=3.4, by_variable=NULL){
        if(is.null(by_variable)){ stop("specify by_variable=<char> arg.") }
        legend_width = WIDTH * 1/8
        ID_labels_max = max(nchar(as.character(curr_msigdb_results_df$ID)))
        ID_labels_width = (ID_labels_max * (WIDTH/3) / 42)
        
        my_variable_labels_max = max(nchar(as.character(curr_msigdb_results_df[[by_variable]])))
        my_variable_width = (my_variable_labels_max * ((WIDTH / 7) / 8))
        base_width = legend_width + ID_labels_width + my_variable_width
        one_panel_x = (WIDTH / 7)
        number_of_x_panels = length(unique(curr_msigdb_results_df$comparison))
        panels_x_width = one_panel_x * number_of_x_panels
        curr_width = base_width + panels_x_width
        
        base_height = HEIGHT * 3/7
        number_of_my_variables = length(unique(curr_msigdb_results_df[[by_variable]]))
        number_of_ID_labels = 0
        for(curr_variable in unique(curr_msigdb_results_df[[by_variable]])){
            number_of_ID_labels = number_of_ID_labels + length(unique(curr_msigdb_results_df[curr_msigdb_results_df[[by_variable]] == curr_variable,]$ID))
        }
        panels_y_height = number_of_ID_labels * ((HEIGHT * 4/7) / 18)
        curr_height = panels_y_height + base_height
        
        return(list(curr_height, curr_width))
    }
    
    archr_proj = ArchR::loadArchRProject(indir)
    
    select_column = "Clusters_TileMatrix_named"
    select_variables_from_column = c(use_groups, background_groups)
    if(!is.null(select_column) & !is.null(select_variables_from_column)){
        #select_variables_from_column = unlist(strsplit(select_variables_from_column, split=","))
        if(select_column %in% colnames(getCellColData(archr_proj))){
            for(select_variable_from_column in select_variables_from_column){
                if(!select_variable_from_column %in% getCellColData(archr_proj)[[select_column]]){
                    stop(select_variable_from_column, " is not found in ", select_column)
                }
            }
        }
    }
    
    message("tfs: ", paste0(tfs, sep=", "))
    message("use_groups: ", paste0(use_groups, sep=", "))
    message("background_groups: ", paste0(background_groups, sep=", "))
    
    curr_outdir = paste0(indir, "/regulon/")
    dir.create(curr_outdir)
    
    # Loading metadata
    mapping = getCellColData(archr_proj, select = c("hash_assignment2", "CNAME", "TEST_ARTICLE", "Treatment", "library", "HTO", "Treatment_HTO", "Treatment_library", "Clusters_TileMatrix_named"))
    mapping$Treatment_library_HTO = paste0(mapping$Treatment_library, "_", mapping$HTO)
    archr_proj$Treatment_library_HTO = mapping$Treatment_library_HTO
    mapping = unique(mapping)
    mapping = data.frame(mapping)
    
    reference_ordinations = readRDS(file=paste0("../OUTPUT/", "reference_ordinations.rds"))
    clusters_named_vColors = reference_ordinations$clusters_named_vColors
    labels_order = reference_ordinations$labels_order
    
    epiregulon_output = readRDS(file=paste0("../OUTPUT/", "regulon/epiregulon.rds"))
    pruned.regulon = epiregulon_output[["pruned.regulon"]]
    regulon.w.corr = epiregulon_output[["regulon.w.corr"]]
    regulon.w.wilcox = epiregulon_output[["regulon.w.wilcox"]]
    score.combine.corr.merge = epiregulon_output[["score.combine.corr.merge"]]
    score.combine.wilcox = epiregulon_output[["score.combine.wilcox"]]
    
    # Load gex matrix
    gene_expression_matrix = getMatrixFromProject(ArchRProj=archr_proj ,useMatrix="GeneExpressionMatrix", useSeqnames=NULL, verbose=TRUE, binarize=FALSE, threads=1)
    gene_expression_matrix = ArchRMatrix2SCE(gene_expression_matrix, rename="normalizedCounts")
    
    #names(assays(gene_expression_matrix)) = "normalizedCounts"
    gene_expression_matrix = as(gene_expression_matrix, "SingleCellExperiment")
    rownames(gene_expression_matrix) = rowData(gene_expression_matrix)$name
    reducedDim(gene_expression_matrix, "UMAP_TileMatrix") = getEmbedding(ArchRProj=archr_proj, embedding="UMAP_TileMatrix", returnDF=TRUE)[colnames(gene_expression_matrix),]
    reducedDim(gene_expression_matrix, "IterativeLSI_TileMatrix") = getReducedDims(ArchRProj=archr_proj, reducedDims="IterativeLSI_TileMatrix")
    
    # load tf binding matrix
    tf_binding_matrix = getMatrixFromProject(ArchRProj=archr_proj, useMatrix="TF_bindingMatrix", useSeqnames=NULL, verbose=TRUE, binarize=FALSE, threads=1)
    #########################################
    # hippo signatures                      #
    #                                       #
    #########################################
    library(genomitory)
    hippo = getFeatureSetCollection("GMTY188:analysis/hippo.gmt.bz2@REVISION-2")
    names(hippo) = hippo@elementMetadata@listData[["name"]]
    
    se = dsassembly::getExperiment("DS000000950")
    search.results = searchFiles("hallmark")
    hallmark <- genomitory::getFeatureSetCollection(search.results[search.results$id=="GMTY42:human/H.gmt.bz2@REVISION-1",]$id)
    names(hallmark) = hallmark@elementMetadata@listData[["name"]]
    EMT = hallmark$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
    EMT.symbol = rowData(se)$symbol[match(EMT, rowData(se)$ID)]
    hippo[["EMT"]] = EMT.symbol
    hippo[["cluster4"]] = cluster4
    
    regulon.hippo.signature = epiregulon:::genesets2regulon(hippo)
    score.combine.hippo = calculateActivity(expMatrix=gene_expression_matrix,
                                             regulon=regulon.hippo.signature,
                                             mode="weight",
                                             method="weightedMean",
                                             exp_assay="normalizedCounts")
    
    pdf(paste0(indir, "/regulon/violin.plots.pdf"), width=8, height=11)
    p1 = plotActivityViolin(score.combine.corr,
                       tf=tfs,
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="corr",
                       nrow=3,
                       boxplot=TRUE)
    p2 = plotActivityViolin(score.combine.corr.merge,
                       tf=tfs,
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="corr merge TF and RE",
                       nrow=3,
                       boxplot=TRUE)
    p3 = plotActivityViolin(score.combine.wilcox,
                       tf=tfs,
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="wilcox",
                       nrow=3,
                       boxplot=TRUE)
    p4 = plotActivityViolin(score.combine.wilcox.motif,
                       tf=tfs,
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="wilcox motif",
                       boxplot=TRUE,
                       nrow=3)
    p5 = plotActivityViolin(score.combine.hippo,
                       tf=c("hippo52","hippo145","mapk","proliferation"),
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       title="hippo signature",
                       boxplot=TRUE,
                       nrow=3)
    print(p1); print(p2); print(p3); print(p4); print(p5)
    dev.off()

    # differential activity using weights
    markers = findDifferentialActivity(score.combine.wilcox,
                                        gene_expression_matrix$Clusters_TileMatrix_named,
                                        pval.type="some",
                                        direction="up",
                                        test.type="t")
    markers.sig = getSigGenes(markers, topgenes=NULL)
    
    curr_tfs=c("ZBTB2",
                    "E2F4",
                    "TRRAP",
                    "ZBTB38",
                    "NFE2L1",
                    "ARID5B",
                    "KLF3",
                    "SMAD1",
                    "EGR1",
                    "ZIC2",
                    "ZNF503",
                    "C17orf49",
                    "TEAD1",
                    "SKP2",
                    "MBD3",
                    "PRMT1",
                    "PBX3",
                    "LMNA",
                    "E2F7",
                    "FOSL1",
                    "VDR",
                    "ATF3",
                    "RUNX2")
    idx = which(gene_expression_matrix$Clusters_TileMatrix_named %in% c("Sensitive_DMSO", "Sensitive_GNE7883", "Resistant_main"))
    gene_expression_matrix2 = gene_expression_matrix[,idx]
    score.combine.wilcox2 = score.combine.wilcox[row.names(score.combine.wilcox) %in% curr_tfs,]
    score.combine.wilcox2 = score.combine.wilcox2[,idx] 
    pdf(paste0(indir, "/regulon/bubble.plots.pdf"), width=6.3, height=7.36)
    p = NULL; p = plotBubble(activity_matrix=score.combine.wilcox2,
               tf=curr_tfs,
               clusters=gene_expression_matrix2$Clusters_TileMatrix_named
               )
    print(p)
    dev.off()
    ggsave(paste0(indir, "/regulon/bubble.plots.png"), plot=p, device="png") 
    
    #umap
    pdf(paste0(indir, "/regulon/umap.plots.pdf"), width=12)
    p1 = NULL; p1 = plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=score.combine.corr,
                    tf=tfs,
                    dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    point_size=0.1,
                    rasterise=TRUE,
                    title="corr",
                    nrow=2)
    p2 = NULL; p2 = plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=score.combine.corr.merge,
                    tf=tfs, dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    point_size=0.1,
                    rasterise=TRUE,
                    title="corr TF RE merge",
                    nrow=2)
    p3 = NULL; p3 = plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=score.combine.wilcox,
                    tf= tfs, dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    point_size=0.1,
                    rasterise=TRUE,
                    title="wilcox",
                    nrow=2)
    print(p1); print(p2); print(p3);
    dev.off()
    
    # Plot gene expresssion
    gex = assay(gene_expression_matrix, "normalizedCounts")
    pdf(file=paste0(indir, "/regulon/geneexpression.pdf"), width=12)
    p1 = NULL; p1 = plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=gex,
                    tf=tfs,
                    dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    legend.label="Gene Exp",
                    point_size=0.1,
                    color=c("grey","blue"),
                    limit=c(0,2),
                    rasterise=TRUE)
    p2 = NULL; p2 = plotActivityViolin(gex,
                       tf=tfs,
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       legend.label="gene expression")
    print(p1); print(p2);
    dev.off()
    
    # Plot chromVar
    chromvar = assay(tf_binding_matrix, "z")
    chromvar = chromvar[tfs,]
    
    pdf(paste0(indir, "/regulon/chromvar.pdf"), width=12)
    p1 = NULL; p1 = plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=chromvar,
                    tf=tfs,
                    dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    legend.label="chromvar",
                    point_size=0.1,
                    color=c("grey","red"),
                    limit=c(0,3),
                    rasterise=TRUE)
    p2 = NULL; p2 = plotActivityViolin(chromvar,
                       tf=tfs,
                       clusters = gene_expression_matrix$Clusters_TileMatrix_named,
                       legend.label="chromvar",
                       nrow=2,
                       boxplot=TRUE)
    print(p1); print(p2);
    dev.off()
    
    # Plot hippo score
    pdf(paste0(indir, "/regulon/hippo.pdf"), width=12)
    p1 = NULL; p1 = plotActivityDim(sce=gene_expression_matrix,
                    activity_matrix=score.combine.hippo,
                    tf=rownames(score.combine.hippo),
                    dimtype="UMAP_TileMatrix",
                    label="Clusters_TileMatrix_named",
                    legend.label="score",
                    point_size=0.1,
                    color=c("grey","brown"),
                    limit=c(1,3),
                    rasterise=TRUE)
    p2 = NULL; p2 = plotActivityViolin(score.combine.hippo,
                       tf=rownames(score.combine.hippo),
                       clusters=gene_expression_matrix$Clusters_TileMatrix_named,
                       legend.label="score",
                       nrow=2,
                       boxplot=TRUE)
    print(p1); print(p2);
    dev.off()

    #######################################
    #  GeneSet enrichment                 #
    #                                     #
    #######################################
    #retrieve genesets
    H = EnrichmentBrowser::getGenesets(org="hsa", db="msigdb", cat="H", gene.id.type="SYMBOL")
    C2 = EnrichmentBrowser::getGenesets(org="hsa", db="msigdb", cat="C2", gene.id.type="SYMBOL")
    C6 = EnrichmentBrowser::getGenesets(org="hsa", db="msigdb", cat="C6", gene.id.type="SYMBOL")
    C8 = EnrichmentBrowser::getGenesets(org="hsa", db="msigdb", cat="C8", gene.id.type="SYMBOL")
    
    # combine genesets and convert genesets to be compatible with enricher
    # only H and C6
    gs = c(H, C6)
    gs.list = do.call(rbind,lapply(names(gs), function(x) {data.frame(gs=x, genes=gs[[x]])}))
    
    enrichresults = regulonEnrich(tfs,
                                   regulon=regulon.w.wilcox,
                                   weight="weight",
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
    write.table(gsea_df, paste0(indir, "/regulon/GSEA_table_epiregulon.tsv"), row.names=F, sep="\t", quote=F)
    
    #plot results
    p1 = enrichPlot(results=enrichresults)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon.pdf"), plot=p1, device="pdf", height=22, width=25, units="in", limitsize=FALSE)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon.png"), plot=p1, device="png", height=22, width=25, units="in", limitsize=FALSE)
    
    # Write gene table for paper
    df = ggplot_build(p1)$plot$data
    df = df[,c("ID", "geneID")]
    gene_list = list()
    for(ID in unique(df$ID)){
        curr_genes = df[df$ID == ID,]$geneID
        curr_genes2 = strsplit(curr_genes, split="/")[[1]]
        gene_list[[ID]] = curr_genes2
    }
    df2 = stack(gene_list)
    colnames(df2) = c("symbol", "pathway")
    df2 = df2[,c("pathway", "symbol")]
    df2 = unique(df2)
    write.table(df2, "../OUTPUT/archr/regulon/gene_table_for_fig3F.tsv", sep="\t", quote=F, row.names=F)
    
    p2 = plotGseaNetwork(tf=tfs, enrichresults, p.adj_cutoff=0.1, ntop_pathways=10)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon_network.pdf"), plot=p2, device="pdf", height=12, width=12, units="in", limitsize=FALSE)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon_network.png"), plot=p2, device="png", height=12, width=12, units="in", limitsize=FALSE)
    
    p1_p2 = ggarrange(p1, p2, ncol=1, nrow=2)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon_combined.pdf"), plot=p1_p2, device="pdf", height=22, width=25, units="in", limitsize=FALSE)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon_combined.png"), plot=p1_p2, device="png", height=22, width=25, units="in", limitsize=FALSE)
    
    # Do with more categories
    gs = c(H, C2, C6, C8)
    gs.list = do.call(rbind,lapply(names(gs), function(x) {data.frame(gs=x, genes=gs[[x]])}))
    
    enrichresults = regulonEnrich(tfs,
                                   regulon=regulon.w.wilcox,
                                   weight ="weight",
                                   weight_cutoff=0,
                                   genesets=gs.list)
    
    #plot results PDFs
    pdf(paste0(indir, "/regulon/GSEA_TFs_regulon_H_C2_C6_C8_enrich.pdf"), width=28, height=10)
    print(enrichPlot(results=enrichresults))
    dev.off()
    
    pdf(paste0(indir, "/regulon/GSEA_TFs_regulon_H_C2_C6_C8_enrich_network.pdf"), width=28, height=10)
    print(plotGseaNetwork(tf=names(enrichresults), enrichresults, p.adj_cutoff=0.1, ntop_pathways=10))
    dev.off()
    
    # PNGs
    p1 = enrichPlot(results=enrichresults)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon_H_C2_C6_C8_enrich.png"), plot=p1, device="png", height=22, width=25, units="in", limitsize=FALSE)
    
    p2 = plotGseaNetwork(tf=tfs, enrichresults, p.adj_cutoff=0.1, ntop_pathways=10)
    ggsave(paste0(indir, "/regulon/GSEA_TFs_regulon_H_C2_C6_C8_enrich_network.png"), plot=p2, device="png", height=22, width=25, units="in", limitsize=FALSE)
    
    # create graphs
    tfs = c("unknown")
    p_network_diffs = list()
    p_top_tfs = list()
    rank_tables = list()
    df_ranks = NULL
    for(j in 1:length(use_groups)){
        use_group = use_groups[j]
        background_group = background_groups[j]
        
        use_network = buildGraph(regulon.w.wilcox, weights= "weight", cluster=use_group)
        background_network = buildGraph(regulon.w.wilcox, weights= "weight", cluster=background_group)
        
        diff_graph_title = paste0("Differential TFs (", use_group, " vs ", background_group, "\nranked by degree centrality")
        diff_graph = buildDiffGraph(use_network, background_network, abs_diff = FALSE)
        diff_graph = addCentrality(diff_graph)
        diff_graph = normalizeCentrality(diff_graph)
        rank_table = rankTfs(diff_graph)
        rank_tables[[paste0(use_group, " vs ", background_group)]] = rank_table
        
        rank_table$comparison = paste0(use_group, "_vs_", background_group)
        
        if(is.null(df_ranks)){
            df_ranks = rank_table
        }else{
            df_ranks = rbind(df_ranks, rank_table)
        }
        
        p = ggplot(rank_table, aes(x=rank, y=centrality)) +
            geom_point() +
            ggtitle(diff_graph_title) +
            ggrepel::geom_text_repel(data=rbind(head(rank_table,10),
                                                tail(rank_table,10)),
                                     max.overlaps=Inf,
                                     aes(label=tf),
                                     nudge_x=0, nudge_y=0, box.padding=0.5, size=3) +
            theme_classic() + theme(title=element_text(size=7))
        p_network_diffs[[paste0(use_group, " vs ", background_group)]] = p
        
        # Then compute similarity.
        diff_graph_filter = igraph::subgraph.edges(diff_graph, E(diff_graph)[E(diff_graph)$weight>0], del=TRUE)
        
        # compute a similarity matrix of all TFs
        similarity_score = calculateJaccardSimilarity(diff_graph_filter)
        
        # Focus on TF of interest
        for(TF_interest in tfs){
            if(TF_interest %in% colnames(similarity_score)){
                message("Processing ", TF_interest)
                similarity_score_TF_interest = similarity_score[, TF_interest]
                similarity_df = data.frame(similarity = head(sort(similarity_score_TF_interest,
                                                                  decreasing=TRUE), 20),
                                           TF = names(head(sort(similarity_score_TF_interest,
                                                                decreasing=TRUE), 20)))
                
                similarity_df$TF = factor(similarity_df$TF, levels=rev(unique(similarity_df$TF)))
                
                # plot top TFs most similar to TF of interest
                topTFplot = ggplot(similarity_df, aes(x=TF, y=similarity)) +
                    geom_bar(stat="identity") +
                    coord_flip() +
                    ggtitle(paste(TF_interest, "")) +
                    theme_classic() + theme(title=element_text(size=6))
                p_top_tfs[[paste0(use_group, " vs ", background_group)]][[TF_interest]][["raw"]] = topTFplot
                
                # account for background
                permute_matrix = permuteGraph(diff_graph_filter, TF_interest, 100, p=1)
                permute_matrix = permute_matrix[names(similarity_score_TF_interest),]
                diff_matrix = similarity_score_TF_interest-rowMeans(permute_matrix)
                diff_matrix_df = data.frame(similarity=head(sort(diff_matrix, decreasing=TRUE), 20),
                                            TF=names(head(sort(diff_matrix, decreasing=TRUE), 20)))
                
                diff_matrix_df$TF = factor(diff_matrix_df$TF, levels=rev(unique(diff_matrix_df$TF)))
                
                # plot top TFs most similar to EBF1
                topTFplot2 = ggplot(diff_matrix_df, aes(x=TF, y=similarity)) +
                    geom_bar(stat="identity") +
                    coord_flip() +
                    ggtitle(paste(TF_interest, "")) +
                    theme_classic() + theme(title=element_text(size=6), axis.text.y=element_text(size=6), axis.text.x=element_text(size=6, angl=90))
                p_top_tfs[[paste0(use_group, " vs ", background_group)]][[TF_interest]][["corrected"]] = topTFplot2
            }else{
                message("Skipping ", TF_interest)
            }
        }
    }
    df_ranks2 = reshape2::dcast(df_ranks, tf ~ comparison, value.var="centrality")
    write.table(df_ranks2, "../OUTPUT/activity_ranks.tsv", sep="\t", quote=F, row.names=F)
    
    saveRDS(rank_tables, file=paste0(curr_outdir, "/rank_tables.rds"))
    
    my_ncol = 3
    my_nrow = ceiling(length(use_groups)/my_ncol)
    p1 = ggarrange(plotlist=p_network_diffs, ncol=my_ncol, nrow=my_nrow)
    
    p2 = list()
    i = 1
    for(curr_name in names(p_network_diffs)){
        j = 1
        p_inner = list()
        for(tf in unique(names(p_top_tfs[[1]]))){
            p_inner[[j]] =  p_top_tfs[[curr_name]][[tf]][["corrected"]] + theme(axis.text.x=element_text(size=5), axis.text.y=element_text(size=5))
            j = j + 1
        }
        p2[[i]] = annotate_figure(ggarrange(plotlist=p_inner, ncol=length(unique(names(p_top_tfs[[1]])))), top=text_grob(paste0(curr_name, " "), color = "black", face="plain", size=7))
        i = i+1
    }
    p2 = ggarrange(plotlist=p2, ncol=1)#, nrow=length(unique(names(p_top_tfs[[1]]))))
    my_ncol = length(unique(names(p_top_tfs[[1]])))
    my_nrow = length(names(p_network_diffs))
    one_panel_y = 2; one_panel_x = 1
    ggsave(paste0(curr_outdir, "/diff_ranks_background_corrected_more", ".pdf"), plot=p2, device="pdf", 
           height=my_nrow*one_panel_y, 
           width=my_ncol*one_panel_x, 
           units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/diff_ranks_background_corrected_more", ".png"), plot=p2, device="png", 
           height=my_nrow*one_panel_y, 
           width=my_ncol*one_panel_x, 
           units="in", limitsize=FALSE)
    
    
    ####################################################
    # Pairwise diff activity                           #
    ####################################################
    # differential activity using weights
    df_activity = NULL
    for(j in 1:length(use_groups)){
        use_group = use_groups[j]
        background_group = background_groups[j]
        
        idx = which(gene_expression_matrix$Clusters_TileMatrix_named %in% c(use_group, background_group))
        gene_expression_matrix_sub = gene_expression_matrix[,idx]
        score.combine.wilcox.sub = score.combine.wilcox[,idx]
        
        # up
        markers_up = findDifferentialActivity(score.combine.wilcox.sub,
                                       gene_expression_matrix_sub$Clusters_TileMatrix_named,
                                       pval.type="some",
                                       direction="up",
                                       test.type="t")
        
        markers_up = markers_up[c(use_group, background_group)]
        markers.sig.up = getSigGenes(markers_up, topgenes=NULL, direction="up")
        if(nrow(data.frame(markers.sig.up)) != 0){
            markers.sig.up$comparison = paste0(use_group, "_vs_", background_group)
            markers.sig.up$direction = "up"
            markers.sig.up = markers.sig.up[markers.sig.up$class == use_group,]
        }
        
        markers_down = findDifferentialActivity(score.combine.wilcox.sub,
                                              gene_expression_matrix_sub$Clusters_TileMatrix_named,
                                              pval.type="some",
                                              direction="down",
                                              test.type="t")
        markers_down = markers_down[c(use_group, background_group)]
        markers.sig.down = getSigGenes(markers_down, topgenes=NULL, direction="down")
        
        if(nrow(data.frame(markers.sig.down)) != 0){
            markers.sig.down$comparison = paste0(use_group, "_vs_", background_group)
            markers.sig.down$direction = "down"
            markers.sig.down = markers.sig.down[markers.sig.down$class == use_group,]
        }
        
        markers.sig = NULL
        if(nrow(data.frame(markers.sig.down)) != 0 & nrow(data.frame(markers.sig.down)) != 0)  {
            markers.sig = rbind(markers.sig.up, markers.sig.down)
        }else if(nrow(data.frame(markers.sig.up)) != 0){
            markers.sig = markers.sig.up
        }else if(nrow(data.frame(markers.sig.down)) != 0){
            markers.sig = rbind(markers.sig.down)
        }
        
        if(is.null(df_activity)){
            df_activity = markers.sig
        }else{
            df_activity = rbind(df_activity, markers.sig) 
        }
    }
    head(df_activity)
    df_activity2 = reshape2::dcast(df_activity, tf ~ comparison, value.var="summary.logFC")
    write.table(df_activity2, "../OUTPUT/activity_diff.tsv", sep="\t", quote=F, row.names=F)

    
    ####################################################
    # Generate final plots:                            #
    ####################################################
    # TF
    #score.combine.wilcox 
    markersTF = tfs
    idx = which(row.names(score.combine.wilcox) %in% markersTF)
    score.combine.wilcox_select = score.combine.wilcox[idx,]
    embed_TFAct = getEmbedding(archr_proj,  embedding="UMAP_TileMatrix")
    dim(score.combine.wilcox_select)
    dim(embed_TFAct)
    score.combine.wilcox_select_t = data.frame(as.matrix(t(score.combine.wilcox_select)))
    identical(row.names(embed_TFAct), row.names(score.combine.wilcox_select_t))
    embed_TFAct = cbind(embed_TFAct, score.combine.wilcox_select_t)
    embed_TFAct$Clusters_TileMatrix_named = archr_proj$Clusters_TileMatrix_named
    colnames(embed_TFAct)[1] = "UMAP_dim_1"
    colnames(embed_TFAct)[2] = "UMAP_dim_2"
    
    embed_TFAct2 = NULL
    local(for(k in 1:length(markersTF)){
        curr_marker_gene = markersTF[k]
        if(!is.null(embed_TFAct[[curr_marker_gene]])){
            tmp_df = embed_TFAct[,c("UMAP_dim_1","UMAP_dim_2", "Clusters_TileMatrix_named")]
            tmp_df = cbind(tmp_df, embed_TFAct[[curr_marker_gene]])
            colnames(tmp_df)[ncol(tmp_df)] = "value"
            tmp_df$symbol = curr_marker_gene
            if(is.null(embed_TFAct2)){
                embed_TFAct2 <<- tmp_df
            }else{
                embed_TFAct2 <<- rbind(embed_TFAct2, tmp_df)
            }
        }
    })
    
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
    clusters_named_vColors
    
    my_comparisons = list(
        c("Sensitive_DMSO", "Sensitive_GNE7883"),
        c("Sensitive_DMSO", "Resistant_main"),
        c("Resistant_main", "Sensitive_GNE7883")
    )
    
    p_tf_list = list()
    p_tf_list2 = list()
    embed_TFAct2 = embed_TFAct2[!embed_TFAct2$Clusters_TileMatrix_named %in% c("Resistant_side", "Sensitive_side"),]
    embed_TFAct2$Clusters_TileMatrix_named = factor(embed_TFAct2$Clusters_TileMatrix_named, unique(embed_TFAct2$Clusters_TileMatrix_named)) 
    local(
    for(gene in unique(embed_TFAct2$symbol)){
        p_TF = NULL; p_TF = ggplot(embed_TFAct2[embed_TFAct2$symbol == gene,], aes(x=UMAP_dim_1, y=UMAP_dim_2, color=value), fill="black") +
            geom_point(size=0.2, alpha=1) +
            scale_color_gradientn(colors=c("#352A86", "#343DAE", "orange", "#F6DA23", "#F8FA0D" )) +
            ggtitle(paste0(gene)) +
            xlab(paste0("UMAP dim 1")) +
            ylab(paste0("UMAP dim 2")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), 
              axis.text=element_text(size=8, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right",
              legend.text=element_text(size=8, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        p_tf_list[[gene]] <<- p_TF
        
        p_TF2 = NULL; p_TF2 = ggplot(embed_TFAct2[embed_TFAct2$symbol == gene,], aes(x=Clusters_TileMatrix_named, y=value, fill=Clusters_TileMatrix_named), fill="black") +
            geom_violin() + geom_boxplot(outlier.shape=NA, width=0.1, color="grey70") + #geom_point(alpha=0.5, size=0.05, position=position_jitter(0.1)) +
            scale_fill_manual(values=clusters_named_vColors) +
            ggtitle(paste0(gene)) +
            xlab(paste0("")) +
            ylab(paste0("")) +
            guides(fill=guide_legend("")) +
            labs(color=NULL) +
            stat_compare_means(comparisons=my_comparisons, method="wilcox.test", paired=F, size=3) +
            theme_minimal() + 
            theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.title=element_text(size=8), axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1), 
                  axis.text=element_text(size=9, angle=0), strip.text.y=element_text(size=8, angle=0), legend.position="right", 
                  legend.text=element_text(size=9, angle=0), plot.title=element_text(size=10), legend.key.size=unit(0.35, 'cm'))
        p_tf_list2[[gene]] <<- p_TF2
    }
    )
    saveRDS(file=paste0(indir, "/regulon/epiregulon_activity_ordinations.rds"), p_tf_list)
    saveRDS(file=paste0(indir, "/regulon/epiregulon_activity_violin.rds"), p_tf_list2)
    
    my_ncol=4; my_nrow=ceiling(length(names(p_tf_list))/my_ncol)
    my_top_height_ratio = predict(my_lm, newdata=data.frame(nrow=my_nrow))[[1]]
    p_tfs = ggarrange(plotlist=p_tf_list, ncol=my_ncol, nrow=my_nrow)             #guides(color = guide_legend(override.aes = list(size = 10))) 
    p_ref = ggarrange(reference_ordinations$init_clusters_named_no_labels + guides(color=guide_legend(ncol=2, override.aes=list(size=2))) + theme(legend.position="right"), 
                      ncol=my_ncol, widths=c(1.225, 0.5, 0.5, 0.5))
    p_tfs = ggarrange(p_tfs, p_ref, ncol=1, nrow=2, heights=c(my_top_height_ratio, 0.5))
    #p_tfs2 = ggarrange(plotlist=p_tf_list2, ncol=my_ncol, nrow=my_nrow)
    p_tfs2 = do.call("grid_arrange_shared_legend", c(p_tf_list2, ncol=my_ncol, nrow=my_nrow, position="right"))
    one_panel_y = 2; one_panel_x = 3;
    ggsave(paste0(curr_outdir, "/epiregulon_activity_umap_TileMatrix_named", ".pdf"), plot=p_tfs, device="pdf", 
           height=(my_nrow+1)*one_panel_y, 
           width=my_ncol*one_panel_x, 
           units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/epiregulon_activity_umap_TileMatrix_named", ".png"), plot=p_tfs, device="png", 
           height=(my_nrow+1)*one_panel_y, 
           width=my_ncol*one_panel_x, 
           units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/epiregulon_activity_violin_TileMatrix_named", ".pdf"), plot=p_tfs2, device="pdf", 
           height=(my_nrow+0)*(one_panel_y+1.5), 
           width=my_ncol*(one_panel_x+0), 
           units="in", limitsize=FALSE)
    ggsave(paste0(curr_outdir, "/epiregulon_activity_violin_TileMatrix_named", ".png"), plot=p_tfs2, device="png", 
           height=(my_nrow+0)*(one_panel_y+1.5), 
           width=my_ncol*(one_panel_x+0), 
           units="in", limitsize=FALSE)
    
    ###########################################
    # Generate GSEA enrichment plots.         #
    #                                         #
    ###########################################
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
        pathways_in_category = EnrichmentBrowser::getGenesets(org="hsa", db="msigdb",
                                          cat=curr_category, gene.id.type="SYMBOL")
        pathways_in_category = names(pathways_in_category)
        
        curr_msigdb_results_df = msigdb_results_df[msigdb_results_df$ID %in% pathways_in_category,]
        order_df = curr_msigdb_results_df %>%
            dplyr::group_by(ID) %>%
            dplyr::summarise_at(vars(p.adjust), list(SUM = ~ sum(.),
                                                SD = ~sd(.))) %>% as.data.frame()
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
            axis.text.x=element_text(angle=90, hjust=1),
            axis.text.y=element_text(angle=0, size=7),
            strip.text.y=element_text(angle=0, hjust=0),
            strip.text.x=element_text(angle=0, hjust=0.5),
            plot.title=element_text(size=8)
        )
    #curr_p
    curr_height = generate_plot_dimensions(curr_msigdb_results_df_all, by_variable="category")[[1]]
    curr_width = generate_plot_dimensions(curr_msigdb_results_df_all, by_variable="category")[[2]]
    ggsave(paste0(indir, "/regulon/GSEA_all_epiregulon_barplot.pdf"), curr_p, width=curr_width+4, height=curr_height, units="in", device="pdf", limitsize=FALSE)
    ggsave(paste0(indir, "/regulon/GSEA_all_epiregulon_barplot.png"), curr_p, width=curr_width+4, height=curr_height, units="in", device="png", limitsize=FALSE)
    
    curr_p_H = ggplot(curr_msigdb_results_df_all[curr_msigdb_results_df_all$category=="H",], aes(x=log10(p.adjust), y=ID, fill=GeneRatio)) + # do -1*
        geom_bar(stat="identity") +
        facet_grid(category ~ symbol, scales="free_y", space="free_y") +
        geom_vline(xintercept=0, linetype="dashed",  color="red", linewidth=0.4) +
        ylab("") +
        ggtitle(paste0("GSEA analysis; selected functions. Category:", "H,C2,C6,C8" ,"\nNES:Normalized Enrichment Score")) +
        theme_minimal() + theme(
            panel.border=element_rect(linewidth=0.5),
            panel.background=element_rect(fill='transparent', color=NA),
            rect = element_rect(fill = "transparent"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
            axis.text.y=element_text(angle=0, size=7),
            strip.text.y=element_text(angle=0, hjust=0),
            strip.text.x=element_text(angle=0, hjust=0.5),
            plot.title=element_text(size=8)
        )
    #curr_p_H
    curr_height = generate_plot_dimensions(curr_msigdb_results_df_all[curr_msigdb_results_df_all$category=="H",], by_variable="category")[[1]]
    curr_width = generate_plot_dimensions(curr_msigdb_results_df_all[curr_msigdb_results_df_all$category=="H",], by_variable="category")[[2]]
    ggsave(paste0(indir, "/regulon/GSEA_H_epiregulon_barplot.pdf"), curr_p_H, width=curr_width+4, height=curr_height, units="in", device="pdf", limitsize=FALSE)
    ggsave(paste0(indir, "/regulon/GSEA_H_epiregulon_barplot.png"), curr_p_H, width=curr_width+4, height=curr_height, units="in", device="png", limitsize=FALSE)
    
    message("Completed epiregulon.R!")
}

is.defined <- function(sym) {
    sym <- deparse(substitute(sym))
    env <- parent.frame()
    exists(sym, env)
}

usage=function(errM) {
    cat("\nUsage : Rscript do_epiregulon.R [option] <Value>\n")
    cat("       -i        : indir containing the archR data structure.")
    cat("       -t        : transcription factors. ex: -t TEAD1,FOSL1")
    cat("       -u        : Use groups. Values have to be separated by a ','. -u and -b args have to contain the same number of string elements separated by a ','")
    cat("       -b        : Background groups. Values have to be separated by a ','. -u and -b args have to contain the same number of string elements separated by a ','")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
    usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        indir=ARG[i+1]
    }else if(ARG[i] == "-t") {
        tfs=ARG[i+1]
    }else if(ARG[i] == "-u") {
        use_groups=ARG[i+1]
    }else if(ARG[i] == "-b") {
        background_groups=ARG[i+1]
    }
}

if(!is.defined(tfs))               { tfs = NULL }
if(!is.defined(indir))             { indir = NULL }
if(!is.defined(use_groups))        { use_groups = NULL }
if(!is.defined(background_groups)) { background_groups = NULL }

do_epiregulon_plots(indir=indir, tfs=tfs, use_groups=use_groups, background_groups=background_groups)
