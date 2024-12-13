library(gp.sa)
library(gp.sa.diff)
library(DataSetDB)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(plyr)
library(edgeR)
library(grid)
library(gridExtra)
library(R.devices)
library(data.table)
library(pheatmap)
library(SummarizedExperiment)
library(edgeR)
library(vegan)
library(kableExtra)
library(ggpubr)
library(eulerr)
library(genomitory)
source("./utils.R")

setwd("/path/to/my/manuscript/RNAseq/NGS4372/analysis/")
outdir <- "OUTPUT/"

se <- getDatasetAsSE("DS000012520")
se$CELL_LINE <- sapply(strsplit(se$SAMPLE_LABEL, " "), "[", 2)
se$TREATMENT_NAME <- stringr::str_trim(gsub("SP_EXP1 ", "", se$TREATMENT_NAME))
se$TREATMENT <- paste0(se$CELL_LINE, "_", se$TREATMENT_NAME)
assays(se)$RPKM <- gp.sa.diff::normalizedRPKM(assays(se)$count, rowData(se)$size)
assays(se)$CPM <- gp.sa.diff::normalizedCPM(assays(se)$count)
expr <- as.matrix(assay(se,"RPKM"))
rownames(expr) <- rowData(se)$symbol

pathway <- getFeatureSetCollection("GMTY188:analysis/hippo.gmt.bz2@REVISION-2")
names(pathway) <- pathway@elementMetadata@listData[["name"]]
pathway[["hippo_KRAS"]] <- c("IRS1","EVA1A","ANLN","UBASH3B","FOSL1","TGM2")
names(pathway)[5] <- "MAPK"
names(pathway)[8] <- "hippo_MAPK"
selected_pathway <- pathway[c(1,2,3,4,5,6,7,8,9)]

gene_index <- match(unique(data.frame(unlist(pathway))[,1]), rowData(se)[,c("symbol")])
nas_index <- which(is.na(gene_index))
gene_index <- gene_index[-c(nas_index)]
backup <- rowData(se)[gene_index, c("symbol")]

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compute signature scores
# Build plots and data frame for all signatures
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
plotlist <- list()
plotlist2 <- list()
final_dfs <- list()
score_table <- matrix(NA, nrow=ncol(se), ncol=length(selected_pathway))
score_table_stdev <- matrix(NA, nrow=ncol(se), ncol=length(selected_pathway))
colnames(score_table) <- names(selected_pathway)
colnames(score_table_stdev) <- names(selected_pathway)
for (i in 1:length(selected_pathway)){
    final_df <- NULL
    for(j in 1:length(1)){
        curr_CELL_LINE <- unique(se$CELL_LINE)[j]
        message("Processing cname: ", curr_CELL_LINE)
        score_name <- names(selected_pathway)[[i]]
        gene_index <- match(selected_pathway[[score_name]], rownames(expr))
        cell_index <- which(se$CELL_LINE %in% unique(se$CELL_LINE))
        selected <- t(expr[gene_index, cell_index])
        Vehicle_index <- colnames(se)[which(se$TREATMENT == "H226_DMSO")]
        Vehicle <- colMeans(selected[Vehicle_index,])
        Vehicle_matrix <- matrix(Vehicle, nrow=nrow(selected), ncol=ncol(selected), byrow=T)
        logFC <- selected - Vehicle_matrix
        # mean
        score <- rowMeans(logFC, na.rm=T)
        score_table[cell_index, score_name] <- score
        # stdev
        score_stdev <- rowSds(logFC, na.rm=F)
        score_table_stdev[cell_index, score_name] <- score_stdev
        table <- data.frame(cbind(colData(se[, cell_index]), score))
        table$TREATMENT_NAME <- factor(table$TREATMENT_NAME, levels=c("DMSO", "G7883"))

        table$signature <- names(selected_pathway)[[i]]
        if(is.null(final_df)){
            final_df <- table
        }else{
            final_df <- rbind(final_df, table)
        }
    }
    final_dfs[[names(selected_pathway)[[i]]]] <- final_df
    final_df <- final_df[final_df$CELL_LINE != "DMSO(83)",]
    dotplot <- ggplot(final_df, aes_string(x="TREATMENT", y="score", color=NULL) ) +
        geom_boxplot(width=0.5) +
        geom_dotplot(binaxis="y", stackdir="center") +
        xlab("") +
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect=element_rect(fill="transparent"),
              axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              plot.title=element_text(size=10, hjust=0.5), legend.position="bottom") +
        ggtitle(paste(score_name)) + ylab ("average logFC") +
        facet_wrap(. ~ CELL_LINE, scales="free_x")
    plotlist[[i]] <- dotplot

    final_df2 <- final_df[final_df$TREATMENT != "G7883(CL3)_DMSO",]
    dotplot2 <- ggplot(final_df2, aes_string(x="TREATMENT", y="score", color=NULL) ) +
        geom_boxplot(width=0.5) +
        geom_dotplot(binaxis="y", stackdir="center") +
        xlab("") +
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect=element_rect(fill="transparent"),
              axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              plot.title=element_text(size=10, hjust=0.5), legend.position="bottom") +
        ggtitle(paste(score_name)) + ylab ("average logFC")
    plotlist2[[i]] <- dotplot2
}

# Build boxplot figures.
p <- ggarrange(plotlist=plotlist, ncol=3, nrow=3)
p2 <- annotate_figure(
    p,
    fig.lab="\n\n*Average logFC values are of each treatment are all relative to the H226_DMSO libraries.",
    fig.lab.pos=c("bottom.left"),
    fig.lab.size=7
)
ggsave("./OUTPUT/figure_signature_all_signatures_for_paper.pdf", plot=p2, device="pdf", units="cm", width=17, height=22, dpi=300)

# Then generate plots with p-values use final_df to do so
# Add p-values to Fig 2C (red vs blue, and grey vs green) a
# H226 + DMSO vs H226-7883R+DMSO
# H226 + GNE-7883 vs H226-7883R + GNE-7883
ylims <- list()
ylims[["hippo145"]] <- c(-1, 0.5)
ylims[["apoptosis"]] <- c(-0.5, 1.2)
ylims[["MAPK"]] <- c(-0.4,2)
ps_stats <- list()
for(signature in c("hippo145", "apoptosis", "MAPK")){
    curr_df <- final_dfs[[signature]]
    curr_df <- curr_df[,c("TREATMENT", "score", "signature")]
    curr_df <- curr_df[curr_df$TREATMENT %in% c("H226_DMSO", "G7883(CL3)_DMSO", "H226_G7883", "G7883(CL3)_G7883"),]
    curr_df$TREATMENT <- factor(curr_df$TREATMENT, levels=c("H226_DMSO", "H226_G7883", "G7883(CL3)_DMSO", "G7883(CL3)_G7883"))
    my_comparisons <- list( c("H226_DMSO", "G7883(CL3)_DMSO"), c("H226_G7883", "G7883(CL3)_G7883") )

    p <- ggplot(curr_df, aes_string(x="TREATMENT", y="score", color=NULL) ) +
        geom_boxplot(width=0.3) +
        geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5) +
        #scale_colour_manual(values=c("black", "blue", "darkgreen", "red")) +
        ylim(ylims[[signature]]) +
        xlab("") +
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect=element_rect(fill="transparent"),
              axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              plot.title=element_text(size=10, hjust=0.5), legend.position="bottom") +
        ggtitle(paste(signature)) + ylab ("average logFC") +
        stat_compare_means(comparisons=my_comparisons, method="t.test", paired=T) +
        stat_compare_means(label.y=45)
    ps_stats[[signature]] <- p
}

p <- ggarrange(ps_stats[[1]],
              ps_stats[[2]],
              ps_stats[[3]], ncol=1, nrow=3)
p
ggsave("./OUTPUT/figure_signature_selected_signatures_for_paper_pvalues.pdf", plot=p, device="pdf", units="in", width=2.75, height=11, dpi=300)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# selected signatures
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Build plots and data frame for selected signatures
selected_pathway <- pathway[c(3,8,6,7)]
gene_index <- match(unique(data.frame(unlist(pathway))[,1]), rowData(se)[,c("symbol")])
nas_index <- which(is.na(gene_index))
gene_index <- gene_index[-c(nas_index)]
backup <- rowData(se)[gene_index, c("symbol")]

plotlist <- list()
plotlist2 <- list()
score_table <- matrix(NA, nrow=ncol(se), ncol=length(selected_pathway))
score_table_stdev <- matrix(NA, nrow=ncol(se), ncol=length(selected_pathway))
colnames(score_table) <- names(selected_pathway)
colnames(score_table_stdev) <- names(selected_pathway)
for (i in 1:length(selected_pathway)){
    final_df <- NULL
    final_df2 <- NULL
    for(j in 1:length(unique(se$CELL_LINE))){
        curr_CELL_LINE <- unique(se$CELL_LINE)[j]
        message("Processing cname: ", curr_CELL_LINE)
        score_name <- names(selected_pathway)[[i]]
        gene_index <- match(selected_pathway[[score_name]], rownames(expr))
        cell_index <- which(se$CELL_LINE %in% unique(se$CELL_LINE))
        selected <- t(expr[gene_index, cell_index])
        Vehicle_index <- which(se$TREATMENT == "H226_DMSO")
        Vehicle <- colMeans(selected[Vehicle_index,])
        Vehicle_matrix <- matrix(Vehicle, nrow=nrow(selected), ncol=ncol(selected), byrow=T)
        logFC <- selected - Vehicle_matrix
        # mean
        score <- rowMeans(logFC, na.rm=T)
        score_table[cell_index, score_name] <- score
        # stdev
        score_stdev <- rowSds(logFC, na.rm=F)
        score_table_stdev[cell_index, score_name] <- score_stdev
        table <- data.frame(cbind(colData(se[, cell_index]), score))
        table$TREATMENT_NAME <- factor(table$TREATMENT_NAME, levels=c("DMSO", "G7883"))

        if(is.null(final_df)){
            final_df <- table
        }else{
            final_df <- rbind(final_df, table)
        }
    }
    final_df <- final_df[final_df$CELL_LINE != "DMSO(83)",]
    final_df <- final_df[!duplicated(final_df), ]
    dotplot <- ggplot(final_df, aes_string(x="TREATMENT", y="score", color=NULL) ) +
        geom_boxplot(width=0.5) +
        geom_dotplot(binaxis="y", stackdir="center") +
        xlab("") +
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect=element_rect(fill="transparent"),
              axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              plot.title=element_text(size=10, hjust=0.5), legend.position="bottom") +
        ggtitle(paste(score_name)) + ylab ("average logFC") +
        facet_wrap(. ~ CELL_LINE, scales="free_x")
    plotlist[[i]] <- dotplot

    final_df2 <- final_df[final_df$TREATMENT != "G7883(CL3)_DMSO",]
    final_df2 <- final_df2[!duplicated(final_df2), ]
    dotplot2 <- ggplot(final_df2, aes_string(x="TREATMENT", y="score", color=NULL) ) +
        geom_boxplot(width=0.5) +
        geom_dotplot(binaxis="y", stackdir="center") +
        xlab("") +
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect=element_rect(fill="transparent"),
              axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              plot.title=element_text(size=10, hjust=0.5), legend.position="bottom") +
        ggtitle(paste(score_name)) + ylab ("average logFC")  #facet_wrap(~ SAM_COMS) #+
    plotlist2[[i]] <- dotplot2
}

# box plots
p <- ggarrange(plotlist=plotlist, ncol=2, nrow=2)
p <- annotate_figure(
    p,
    fig.lab="\n\n*Average logFC values are of each treatment are all relative to the H226_DMSO libraries.",
    fig.lab.pos=c("bottom.left"),
    fig.lab.size=5.5
)
ggsave("./OUTPUT/figure_signature_for_paper_v1.pdf", plot=p, device="pdf", units="cm", width=11, height=15, dpi=300)
p2 <- ggarrange(plotlist=plotlist2, ncol=2, nrow=2)
p2 <- annotate_figure(
    p2,
    fig.lab="\n\n*Average logFC values are of each treatment are all relative to the H226_DMSO libraries.",
    fig.lab.pos=c("bottom.left"),
    fig.lab.size=5.5
)
ggsave("./OUTPUT/figure_signature_for_paper_v2.pdf", plot=p2, device="pdf", units="cm", width=9, height=15, dpi=300)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Heatmaps
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
curr_list <- list()
x=1
mapping <- data.frame(colData(se)[,c("TREATMENT"),drop=F])
mapping <- mapping[mapping$TREATMENT %in% c("H226_DMSO", "H226_G7883", "G7883(CL3)_G7883", "G7883(CL3)_DMSO"),,drop=F]
colnames(mapping)[1] <- "Clones"
mapping2 <- mapping
for(j in 1:ncol(mapping2)){
    curr_col_name <- names(mapping2)[j]
    curr_var_names <- unique(mapping2[,j])
    curr_colors <- vColors[x:(x+(length(curr_var_names)-1))]
    names(curr_colors) <- curr_var_names
    x <- length(curr_var_names) + 1 + x
    curr_list[[curr_col_name]] <- curr_colors
}
print(curr_list)

rpkm_df <- NULL
pathway[["additional_1"]] <- c("VIM", "EGFR", "CD44")
signatures <- names(pathway)
for(k in 1:length(signatures)){
    curr_signature <- signatures[k]
    expr2 <- expr[row.names(expr) %in% pathway[[curr_signature]], colnames(expr) %in% row.names(mapping2)]
    dim(expr2)
    gene_std <- apply(na.omit(expr2), 1, sd)
    expr2 <- expr2[!rownames(expr2) %in% names(which(gene_std == 0)),]

    tmp_df <- data.frame(expr2)
    tmp_df$signature <- curr_signature
    if(is.null(rpkm_df)){
        rpkm_df <- tmp_df
    }else{
        rpkm_df <- rbind(rpkm_df, tmp_df)
    }

    pheatmap::pheatmap(
        expr2,
        main=paste0("heatmap (RPKM) of the ", curr_signature, " gene signature."),
        show_rownames=T,
        show_colnames=F,
        border_color="gray",
        cellwidth=9,
        cellheight=10,
        cluster_cols=F,cluster_rows=T,
        annotation=mapping2,
        annotation_colors=curr_list,
        color=colorRampPalette(c("#2C397F", "white", "darkorange4"))(100),
        clustering_method="average",
        fontsize=7, fontsize_row=9, fontsize_col=9,
        file=paste0("./OUTPUT/figure_heatmap_", curr_signature, "_for_paper.pdf"),
        scale="row"
    )
}
dev.off()

rpkm_df <- NULL
pathway[["AP1s"]] <- c("FOS", "FOSL1", "FOSL2", "FOSB", "JUN", "JUNB", "JUND")
signatures <- names(pathway)
for(k in 1:length(signatures)){
    curr_signature <- signatures[k]
    expr2 <- expr[row.names(expr) %in% pathway[[curr_signature]], colnames(expr) %in% row.names(mapping)]
    dim(expr2)
    gene_std <- apply(na.omit(expr2), 1, sd)
    expr2 <- expr2[!rownames(expr2) %in% names(which(gene_std == 0)),]

    tmp_df <- data.frame(expr2)
    tmp_df$signature <- curr_signature
    tmp_df$symbol <- row.names(tmp_df)
    row.names(tmp_df) <- seq(1, nrow(tmp_df), 1)
    if(is.null(rpkm_df)){
        rpkm_df <- tmp_df
    }else{
        rpkm_df <- rbind(rpkm_df, tmp_df)
    }
}
mapping <- data.frame(colData(se)[,c("TREATMENT", "CELL_LINE")])
df <- melt(rpkm_df)
df <- merge(df, mapping, by.x="variable", by.y="row.names")
df <- df %>%
    dplyr::group_by(signature, symbol, TREATMENT) %>%
    dplyr::summarise_at(vars(value), list(Mean = ~mean(.),
                                          Stdev = ~sd(.))) %>% as.data.frame() %>% as.data.frame()

p <- ggplot(df, aes(x=TREATMENT, y=Mean)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
    facet_grid(. ~ signature) +
    theme_minimal() +
    ylab("average RPKM") +
    labs(caption="Each point = average (of RPKMs) of 3 replicates of one gene belonging to its gene signature set.\nEach panel = one gene signature set.") +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect=element_rect(fill="transparent"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(angle=0, size=7),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8),
        plot.caption=element_text(hjust=0, size=7)
    )
p
ggsave(paste0(outdir,"/RPKM_boxplot_gene_signatures_for_paper.pdf"), p, width=14.17, height=3.83, units="in", device="pdf")

df$TREATMENT <- factor(df$TREATMENT, levels=c("H226_DMSO", "H226_G7883", "G7883(CL3)_DMSO", "G7883(CL3)_G7883"))
p <- ggplot(df[df$signature == "AP1s",], aes(x=symbol, y=Mean, fill=TREATMENT)) +
    geom_bar(stat="identity", position=position_dodge(0.9), color="black") +
    geom_errorbar(aes(x=symbol, ymin=Mean-Stdev, ymax=Mean+Stdev), size=0.5,
                   width=.25, position=position_dodge(.9)) +
    scale_fill_manual(values=curr_list[[1]]) +
    facet_grid(. ~ signature) +
    theme_minimal() +
    ylab("average RPKM") +
    labs(caption="Each point = average (of RPKMs) of 3 replicates of one gene belonging to its gene signature set.\nEach panel = one gene signature set.") +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect=element_rect(fill="transparent"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(angle=0, size=9),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8),
        plot.caption=element_text(hjust=0, size=7)
    )
p
ggsave(paste0(outdir,"/RPKM_boxplot_gene_signatures_for_paper_AP1s.pdf"), p, width=5, height=3.9, units="in", device="pdf")

df <- melt(rpkm_df)
df <- merge(df, mapping, by.x="variable", by.y="row.names")
df <- df[df$signature == "additional_1",]

p <- ggplot(df, aes(x=TREATMENT, y=Mean)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
    facet_grid(. ~ symbol) +
    theme_minimal() +
    ylab("average RPKM") +
    labs(caption="Each point = average (of RPKMs) of 3 replicates of one gene belonging to its gene signature set.\nEach panel = one gene signature set.") +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect=element_rect(fill="transparent"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(angle=0, size=7),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8),
        plot.caption=element_text(hjust=0, size=7)
    )
p

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# GSEA 
#                                    
# H226_G7883 vs H226_DMSO            
# G7883(CL3)_G7883 H226_G7883        
# H226-P-DMSO vs H226-7883R CL3 7883 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom_output <- readRDS(paste0(outdir, "voom.output_JT.rds"))
voom_output_selected <- list()
voom_output_selected$`NCI-H226`[["Difference between `H226_G7883` vs `H226_DMSO`"]] <- voom_output$`NCI-H226`[["Difference between `H226_G7883` vs `H226_DMSO`"]]
voom_output_selected$`NCI-H226`[["Difference between `G7883(CL3)_G7883` vs `H226_G7883"]] <- voom_output$`NCI-H226`[["Difference between `G7883(CL3)_G7883` vs `H226_G7883`"]]
voom_output_selected$`NCI-H226`[["Difference between `G7883(CL3)_DMSO` vs `H226_DMSO`"]] <- voom_output$`NCI-H226`[["Difference between `G7883(CL3)_DMSO` vs `H226_DMSO`"]]
voom_output_selected$`NCI-H226`[["Difference between `G7883(CL3)_DMSO` vs `H226_DMSO`"]] <- voom_output$`NCI-H226`[["Difference between `G7883(CL3)_DMSO` vs `H226_DMSO`"]]
voom_output_selected$`NCI-H226`[["Difference between `G7883(CL3)_G7883` vs `H226_DMSO`"]] <- voom_output$`NCI-H226`[["Difference between `G7883(CL3)_G7883` vs `H226_DMSO`"]]

# KEGG_ABC_TRANSPORTERS
#C2: Curated
#CP: Canonical Pathways
#CP:KEGG_LEGACY: KEGG Legacy Pathways
#M11911
msigdb_results <- list()
msigplot <- list()
msigdb_results_df <- NULL
cells <- unique(se$CNAME)
x <- 0
qvalue_cutoff <- 0.05
local(for(i in 1:length(cells)){
    cell <- cells[i]

    category=c("H","C2", "C5", "C6","C8")

    for (j in 1:length(category)) {
        if(category[j] == "C2_KEGG"){
            m_t2g <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
            m_t2g <- m_t2g[m_t2g$gs_id == "M11911",]
            m_t2g <- m_t2g %>% dplyr::select(gs_name, gene_symbol)
        }else{
            m_t2g <- msigdbr(species="Homo sapiens", category=category[j]) %>%
                dplyr::select(gs_name, gene_symbol)
        }

        for (k in 1:length(voom_output_selected[[cell]])) {
            data1 <- voom_output_selected[[cell]][[k]]$LogFC
            names(data1) <- rowData(se)$symbol
            data1 <- sort(data1, decreasing=TRUE)
            data1_clean <- data1[which(!is.na(names(data1)))]
            message(paste("analyzing", category[j],  names(voom_output_selected[[cell]][k]), cell))

            msigdb <- GSEA(data1_clean, TERM2GENE=m_t2g)
            # NES: (normalized enrichment score)
            if(nrow(msigdb@result) > 0){
                curr_results <- msigdb@result[order(msigdb@result$NES), ]
                curr_results$core_enrichment <- NULL
                curr_results$direction <- ifelse(curr_results$NES > 0, "UP", "DOWN")
                curr_results$cell <- cell
                curr_results$category <- category[j]
                curr_results$comparison <- names(voom_output_selected[[cell]][k])
                curr_results <- curr_results[curr_results$qvalue < qvalue_cutoff,] #0.05
                msigdb_results[[ cell ]][[names(voom_output_selected[[cell]][k])]][[ category[j] ]][[ k ]] <- curr_results
            }else{
                message("...no term enriched under specific pvalueCutoff for ", paste0(cell, " - ",  names(voom_output_selected[[cell]][k])))
            }

            if(x == 0){
                msigdb_results_df <<- curr_results
            }else{
                msigdb_results_df <<- rbind(msigdb_results_df, curr_results)
            }
            x <- x + 1
        }
    }
})

# Write table
msigdb_results_df2 <- msigdb_results_df
msigdb_results_df2$comparison <- gsub("\n", " ", msigdb_results_df2$comparison)
write.table(msigdb_results_df2, "./OUTPUT/GSEA_table.tsv", sep="\t", row.names=F, quote=F)

# Then process the resulting df to generate plots.
head(msigdb_results_df)
msigdb_results_df$comparison <- gsub("Difference between ", "", msigdb_results_df$comparison)
msigdb_results_df$comparison <- gsub("`", "", msigdb_results_df$comparison)
msigdb_results_df$comparison <- gsub(" vs ", "\nvs\n", msigdb_results_df$comparison)
msigdb_results_df$comparison <- factor(msigdb_results_df$comparison, levels=unique(msigdb_results_df$comparison))

# order by NES
curr_msigdb_results_df_all <- NULL
for(curr_category in c("H", "C2", "C5", "C6", "C8")){
    curr_msigdb_results_df <- msigdb_results_df[msigdb_results_df$category == curr_category,]
    order_df <- curr_msigdb_results_df %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise_at(vars(NES), list(SUM = ~ sum(.),
                                            SD = ~sd(.))) %>%
        as.data.frame()
    order_df <- order_df[order(-order_df$SUM),]

    curr_msigdb_results_df$ID <- factor(curr_msigdb_results_df$ID, levels=unique(order_df$ID))
    IDs_vector <- unique(curr_msigdb_results_df$ID)
    if(is.null(curr_msigdb_results_df_all)){
        curr_msigdb_results_df_all <- curr_msigdb_results_df
    }else{
        curr_msigdb_results_df_all <- rbind(curr_msigdb_results_df_all, curr_msigdb_results_df)
    }

    curr_msigdb_results_df$comparison <- factor(curr_msigdb_results_df$comparison, levels=c("G7883(CL3)_G7883\nvs\nH226_G7883", "H226_G7883\nvs\nH226_DMSO", "G7883(CL3)_DMSO\nvs\nH226_DMSO", "G7883(CL3)_G7883\nvs\nH226_DMSO"))
    curr_p <- ggplot(curr_msigdb_results_df, aes(x=NES, y=ID, fill=qvalue)) +
        geom_bar(stat="identity") +
        facet_grid(. ~ comparison, scales="free_y", space="free_y") +
        geom_vline(xintercept=0, linetype="dashed",  color="red", size=0.4) +
        ylab("") +
        ggtitle(paste0("GSEA analysis; selected functions. Category:", curr_category,"\nNES:Normalized Enrichment Score")) +
        theme_minimal() + theme(
            panel.border=element_rect(linewidth=0.5),
            panel.background=element_rect(fill='transparent', color=NA),
            rect=element_rect(fill="transparent"),
            axis.text.x=element_text(angle=90),
            axis.text.y=element_text(angle=0, size=7),
            strip.text.y=element_text(angle=0, hjust=0),
            strip.text.x=element_text(angle=0, hjust=0.5),
            plot.title=element_text(size=8)
        )
    curr_height <- generate_plot_dimensions(curr_msigdb_results_df, by_variable="cell")[[1]]
    curr_width <- generate_plot_dimensions(curr_msigdb_results_df, by_variable="cell")[[2]]

    ggsave(paste0(outdir,"/GSEA_", curr_category, "_for_paper2.pdf"), curr_p, width=curr_width, height=curr_height,
           units="in", device="pdf", limitsize=FALSE)

    df <- curr_msigdb_results_df[curr_msigdb_results_df$comparison == "H226_G7883\nvs\nH226_DMSO",]
    signatures <- as.character(unique(df$ID))
    signatures <- signatures[!signatures %in% c("HALLMARK_G2M_CHECKPOINT",
                                               "HALLMARK_ESTROGEN_RESPONSE_LATE",
                                               "HALLMARK_BILE_ACID_METABOLISM",
                                               "HALLMARK_COMPLEMENT",
                                               "HALLMARK_HEME_METABOLISM",
                                               "HALLMARK_XENOBIOTIC_METABOLISM",
                                               "HALLMARK_COAGULATION"
                                               )]
    m_t2g <- msigdbr(species="Homo sapiens", category="H") %>%
        dplyr::select(gs_name, gene_symbol)
    m_t2g <- m_t2g[m_t2g$gs_name %in% signatures,]
    m_t2g <- unique(m_t2g)
    write.table(m_t2g, paste0(outdir, "gene_tables_for_paper_fig1E_and_2G.tsv"), sep="\t", quote=F, row.names=F)
}

selected_gseas <- c("OXIDATIVE_PHOSPHORYLATION",
                   "GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",
                   "ORGANELLE_INNER_MEMBRANE",
                   "RESPONSE_TO_CYTOKINE",
                   "DEFENSE_RESPONSE",
                   "CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
                   "ADAPTATIVE_IMMUNE_RESPONSE",
                   "DEFENSE_RESPONSE_TO_OTHER_ORGANISM",
                   "INNATE_IMMUNE_RESPONSE",
                   "OXIDANT",
                   "WNT")

curr_msigdb_results_df_all_selected <- curr_msigdb_results_df_all[grepl(paste(selected_gseas, collapse="|"), curr_msigdb_results_df_all$Description),]
curr_msigdb_results_df_all_selected$comparison <- gsub("Difference between `", "", curr_msigdb_results_df_all_selected$comparison)
curr_msigdb_results_df_all_selected$comparison <- gsub("`", "", curr_msigdb_results_df_all_selected$comparison)
curr_msigdb_results_df_all_selected$comparison <- factor(curr_msigdb_results_df_all_selected$comparison, levels=c("G7883(CL3)_G7883\nvs\nH226_G7883", "H226_G7883\nvs\nH226_DMSO", "G7883(CL3)_DMSO\nvs\nH226_DMSO", "G7883(CL3)_G7883\nvs\nH226_DMSO"))
curr_p <- ggplot(curr_msigdb_results_df_all_selected, aes(x=NES, y=ID, fill=qvalue)) +
    geom_bar(stat="identity") +
    facet_grid(. ~ comparison, scales="free_y", space="free_y") +
    geom_vline(xintercept=0, linetype="dashed",  color="red", size=0.4) +
    ylab("") +
    ggtitle(paste0("GSEA analysis; selected functions. Category:", "[H] and [C5]","\nNES:Normalized Enrichment Score")) +
    theme_minimal() + theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect=element_rect(fill="transparent"),
        axis.text.x=element_text(angle=90),
        axis.text.y=element_text(angle=0, size=6.5),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8)
    )
curr_p

ggsave(paste0(outdir,"/GSEA_", "selection_for_paper", ".pdf"), curr_p, width=8, height=3.55, units="in", device="pdf")

