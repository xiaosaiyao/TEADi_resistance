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
library(aplot)
source("analysis/utils.R")

setwd("/path/to/my/manuscript/Quantseq/DS000015488/")
outdir="OUTPUT/"

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

    curr_height = panels_y_height + base_height

    return(list(curr_height, curr_width))
}


se <- getDatasetAsSE("DS000015488")
se <- se[,se$cell_name == "MSTO-211H"]
se <- se[,grep("DMSO|7883", se$SAM_COMS)]
se$TEST_ARTICLE = "bla"
se$TEST_ARTICLE = ifelse(grepl("DMSO", se$SAM_COMS), "DMSO", se$TEST_ARTICLE)
se$TEST_ARTICLE = ifelse(grepl("7883", se$SAM_COMS), "GNE-7883", se$TEST_ARTICLE)
se$CELLTYPE = ifelse(grepl("par", se$SAM_COMS), "MSTO_211H",       se$SAM_COMS)
se$CELLTYPE = ifelse(grepl("res", se$SAM_COMS), "MSTO_211H-7883R", se$CELLTYPE)
se$TREATMENT = paste0(se$CELLTYPE, " + ", se$TEST_ARTICLE)
colData(se)
outdir="OUTPUT/"

# I kept the $RPKM variable name, but these are really CPMs.
assays(se)$RPKM = gp.sa.diff::normalizedCPM(assays(se)$count)
expr=as.matrix(assay(se,"RPKM"))
rownames(expr)=rowData(se)$symbol
head(expr)

pathway <- getFeatureSetCollection("GMTY188:analysis/hippo.gmt.bz2@REVISION-2")
names(pathway) <- pathway@elementMetadata@listData[["name"]]
pathway[["hippo_KRAS"]] = c("IRS1","EVA1A","ANLN","UBASH3B","FOSL1","TGM2")
names(pathway)[5] = "MAPK"
names(pathway)[8] = "hippo_MAPK"

selected_pathway = pathway[c(1,2,3,4,5,6,7,8,9)]

gene_index = match(unique(data.frame(unlist(pathway))[,1]), rowData(se)[,c("symbol")])
nas_index = which(is.na(gene_index))
gene_index = gene_index[-c(nas_index)]
backup = rowData(se)[gene_index, c("symbol")]

# Build plots and data frame for all signatures
se$CELL_LINE = se$CELLTYPE
plotlist = list()
plotlist2 = list()
final_dfs = list()
score_table = matrix(NA, nrow=ncol(se), ncol=length(selected_pathway)) # row = 36 samples, col = number of signature (8)
score_table_stdev = matrix(NA, nrow=ncol(se), ncol=length(selected_pathway)) # row = 36 samples, col = number of signature (8)
colnames(score_table) = names(selected_pathway)
colnames(score_table_stdev) = names(selected_pathway)
x = 1
for (i in 1:length(selected_pathway)){
    final_df = NULL
    for(j in 1:length(unique(se$CELL_LINE))){ # so since we always compute against H226_DMSO, no need to loop through cell line types.
        curr_CELL_LINE = unique(se$CELL_LINE)[j]
        message("Processing cname: ", curr_CELL_LINE)
        score_name = names(selected_pathway)[[i]]
        gene_index = match(selected_pathway[[score_name]], rownames(expr))
        cell_index = which(se$CELL_LINE == curr_CELL_LINE)
        cell_index_ctrl = which(se$TREATMENT == "MSTO_211H + DMSO")
        cell_index = unique(c(cell_index, cell_index_ctrl))
        selected = t(expr[gene_index, cell_index])
        Vehicle_index = colnames(se)[which(se$TEST_ARTICLE == "DMSO" & se$CELL_LINE == "MSTO_211H")]
        Vehicle = colMeans(selected[Vehicle_index,])
        Vehicle_matrix = matrix(Vehicle, nrow=nrow(selected), ncol=ncol(selected), byrow=T)
        logFC = selected - Vehicle_matrix
        # mean
        score = rowMeans(logFC, na.rm=T)
        score_table[cell_index, score_name] = score
        # stdev
        score_stdev = rowSds(logFC, na.rm=F)
        score_table_stdev[cell_index, score_name] = score_stdev
        table = data.frame(cbind(colData(se[, cell_index]), score))
        table$TEST_ARTICLE = factor(table$TEST_ARTICLE, levels=c("DMSO", "GNE-7883"))

        table$signature = names(selected_pathway)[[i]]
        if(is.null(final_df)){
            final_df = table
        }else{
            final_df = rbind(final_df, table)
        }
    }
    final_df = unique(final_df)
    final_dfs[[names(selected_pathway)[[i]]]] = final_df
    final_df2 = final_df
    dotplot = ggplot(final_df2, aes(x=TREATMENT, y=score, color=NULL) ) +
        geom_boxplot(width = 0.5) +
        geom_dotplot(binaxis = "y", stackdir = "center") +
        facet_wrap(. ~ CELL_LINE) +
        xlab("") +
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect = element_rect(fill = "transparent"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.x = element_text(size=7),
              plot.title = element_text(size=10, hjust = 0.5), legend.position = "bottom") +
        ggtitle(paste(score_name)) + ylab ("average logFC") +  #facet_wrap(~ SAM_COMS) #+
        facet_wrap(. ~ CELL_LINE, scales = "free_x") #, space='free')
    plotlist[[x]] = dotplot

    x = x + 1
}

# box plots
p = ggarrange(plotlist=c(plotlist), ncol=3, nrow=3)
p2 = annotate_figure(
    p,
    fig.lab = "\n\n*Average logFC values are of each treatment are all relative to the MSTO-211H DMSO libraries.",
    fig.lab.pos = c("bottom.left"),
    fig.lab.size=7
)
ggsave("./OUTPUT/figure_signature_all_signatures_for_paper.pdf", plot=p2, device="pdf", units="in", 
       width=8, height=12, dpi=300)

CELL_LINE = "MSTO-211H"
ylims = list()
ylims[["hippo145"]] = c(-0.45, 0.75)
ylims[["apoptosis"]] = c(-0.5, 2)
ylims[["MAPK"]] = c(-0.4,0.8)
ps_stats = list()
for(signature in c("hippo145", "apoptosis", "MAPK")){
    curr_df = final_dfs[[signature]]
   
        curr_df2 = curr_df
        curr_df2 = curr_df2[,c("TREATMENT", "TEST_ARTICLE", "score", "signature")]
        curr_df2 = curr_df2[curr_df2$TEST_ARTICLE %in% c("DMSO", "GNE-7883"),]
        curr_df2$TEST_ARTICLE = factor(curr_df2$TEST_ARTICLE, levels=c("DMSO", "GNE-7883"))
        my_comparisons = list()
        
        my_comparisons[[1]] = c("MSTO_211H + GNE-7883",          "MSTO_211H + DMSO")
        my_comparisons[[2]] = c("MSTO_211H-7883R + GNE-7883",    "MSTO_211H-7883R + DMSO")
        my_comparisons[[3]] = c("MSTO_211H + DMSO",              "MSTO_211H-7883R + DMSO")
        my_comparisons[[4]] = c("MSTO_211H + DMSO",              "MSTO_211H-7883R + GNE-7883")
        my_comparisons[[5]] = c("MSTO_211H + GNE-7883",          "MSTO_211H-7883R + GNE-7883")
    
        p = ggplot(curr_df2, aes(x=TREATMENT, y=score, color=NULL) ) +
            geom_boxplot(width = 0.3) +
            geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.5) +
            ylim(ylims[[signature]]) +
            xlab("") +
            theme_minimal() +
            theme(panel.border=element_rect(linewidth=0.5),
                  panel.background=element_rect(fill='transparent', color=NA),
                  rect = element_rect(fill = "transparent"),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  strip.text.x = element_text(size=7),
                  plot.title = element_text(size=10, hjust = 0.5), legend.position = "bottom") +
            ggtitle(paste(signature)) + ylab ("average logFC") + #facet_wrap(~ SAM_COMS) #+
            stat_compare_means(comparisons=my_comparisons, method="t.test", paired=T, size=2.5) +
            stat_compare_means(label.y = 45)
        ps_stats[[paste0(signature, "_", CELL_LINE)]] = p
}

p = ggarrange(ps_stats[[1]],
              ps_stats[[2]],
              ps_stats[[3]], ncol=3, nrow=1, align="h")
ggsave("./OUTPUT/figure_signature_selected_signatures_for_paper_pvalues.pdf", plot=p, device="pdf", 
       units="in", width=6.6, height=4.4, dpi=300)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# selected signatures
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Build plots and data frame for selected signatures
selected_pathway = pathway[c(3,8,6,7)]
gene_index = match(unique(data.frame(unlist(pathway))[,1]), rowData(se)[,c("symbol")])
nas_index = which(is.na(gene_index))
gene_index = gene_index[-c(nas_index)]
backup = rowData(se)[gene_index, c("symbol")]

gene_index = match(unique(data.frame(unlist(pathway))[,1]), rowData(se)[,c("symbol")])
nas_index = which(is.na(gene_index))
gene_index = gene_index[-c(nas_index)]
backup = rowData(se)[gene_index, c("symbol")]

plotlist = list()
plotlist2 = list()
score_table = matrix(NA, nrow=ncol(se), ncol=length(selected_pathway))
score_table_stdev = matrix(NA, nrow=ncol(se), ncol=length(selected_pathway))
colnames(score_table) = names(selected_pathway)
colnames(score_table_stdev) = names(selected_pathway)
for (i in 1:length(selected_pathway)){
    final_df = NULL
    for(j in 1:length(unique(se$CELL_LINE))){
        curr_CELL_LINE = unique(se$CELL_LINE)[j]
        message("Processing cname: ", curr_CELL_LINE)
        score_name = names(selected_pathway)[[i]]
        gene_index = match(selected_pathway[[score_name]], rownames(expr))
        gene_index = gene_index[!is.na(gene_index)]
        # here do this, because we want to compute everything against MSTO_211H + DMSO 
        cell_index = which(se$CELL_LINE == curr_CELL_LINE)
        cell_index_ctrl = which(se$TREATMENT == "MSTO_211H + DMSO")
        cell_index = unique(c(cell_index, cell_index_ctrl))
        selected = t(expr[gene_index, cell_index])
        Vehicle_index = colnames(se)[which(se$TEST_ARTICLE == "DMSO" & se$CELL_LINE == "MSTO_211H")]
        Vehicle = colMeans(selected[Vehicle_index,])
        Vehicle_matrix = matrix(Vehicle, nrow=nrow(selected), ncol=ncol(selected), byrow=T)
        logFC = selected - Vehicle_matrix
        # mean
        score = rowMeans(logFC, na.rm=T)
        score_table[cell_index, score_name] = score
        # stdev
        score_stdev = rowSds(logFC, na.rm=F)
        score_table_stdev[cell_index, score_name] = score_stdev
        table = data.frame(cbind(colData(se[, cell_index]), score)) # A merge would have technically been safer, but ok.
        table$TEST_ARTICLE = factor(table$TEST_ARTICLE, levels=c("DMSO", "GNE-7883"))

        if(is.null(final_df)){
            final_df = table
        }else{
            final_df = rbind(final_df, table)
        }
    }
    final_df = unique(final_df)
    final_df2 = final_df
    dotplot = ggplot(final_df2, aes(x=TREATMENT, y=score, color=NULL) ) +
        geom_boxplot(width = 0.5) +
        geom_dotplot(binaxis = "y", stackdir = "center") +
        xlab("") +
        facet_wrap(. ~ CELL_LINE) + 
        theme_minimal() +
        theme(panel.border=element_rect(linewidth=0.5),
              panel.background=element_rect(fill='transparent', color=NA),
              rect = element_rect(fill = "transparent"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.x = element_text(size=7),
              plot.title = element_text(size=10, hjust = 0.5), legend.position = "bottom") +
        ggtitle(paste(score_name)) + ylab ("average logFC") +  #facet_wrap(~ SAM_COMS) #+
        facet_wrap(. ~ CELL_LINE, scales = "free_x") #, space='free')
    plotlist[[i]] = dotplot
}

# box plots
p = ggarrange(plotlist=plotlist, ncol=2, nrow=2)
p = annotate_figure(
    p,
    fig.lab = "\n\n*Average logFC values are of each treatment are all relative to the H226_DMSO libraries.",
    fig.lab.pos = c("bottom.left"),
    fig.lab.size=5.5
)
ggsave("./OUTPUT/figure_signature_for_paper_v1.pdf", plot=p, device="pdf", units="in", width=6, height=9, dpi=300)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Heatmap
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
curr_list = list()
x=1
mapping = data.frame(colData(se)[,c("TREATMENT"), drop=F])
mapping2 = mapping
for(j in 1:ncol(mapping2)){
    curr_col_name = names(mapping2)[j]
    curr_var_names = unique(mapping2[,j])
    curr_colors = vColors[x:(x+(length(curr_var_names)-1))]
    names(curr_colors) = curr_var_names
    x = length(curr_var_names) + 1 + x
    curr_list[[curr_col_name]] = curr_colors
}
print(curr_list)
curr_list[["TREATMENT"]][["MSTO_211H + GNE-7883"]] = "gray"
curr_list[["TREATMENT"]][["MSTO_211H + DMSO"]] = "red"
curr_list[["TREATMENT"]][["MSTO_211H-7883R + GNE-7883"]] = "green"
curr_list[["TREATMENT"]][["MSTO_211H-7883R + DMSO"]] = "blue"

mapping2 = mapping2[order(mapping$TREATMENT),,drop=F]

rpkm_df = NULL
pathway[["additional_1"]] = c("VIM", "EGFR", "CD44")
signatures = names(pathway)
for(k in 1:length(signatures)){
    curr_signature = signatures[k]
    expr2 = expr[row.names(expr) %in% pathway[[curr_signature]], colnames(expr) %in% row.names(mapping2)]
    dim(expr2)
    gene_std = apply(na.omit(expr2), 1, sd)
    expr2 = expr2[!rownames(expr2) %in% names(which(gene_std == 0)),]

    tmp_df = data.frame(expr2)
    tmp_df$signature = curr_signature
    if(is.null(rpkm_df)){
        rpkm_df = tmp_df
    }else{
        rpkm_df = rbind(rpkm_df, tmp_df)
    }

    pheatmap::pheatmap(
        expr2[,row.names(mapping2)],
        main=paste0("heatmap (RPKM) of the ", curr_signature, " gene signature."),
        show_rownames=T,
        show_colnames=F,
        border_color="gray",
        cellwidth=9,
        cellheight=10,
        cluster_cols=F,cluster_rows=T,
        annotation=mapping2,
        annotation_colors=curr_list,
        #color=colorRampPalette(c("#2C397F", "#46B2E5", "#9FCE63", "#F0E921", "#EE3128", "#552E31"))(100),
        color=colorRampPalette(c("#2C397F", "white", "darkorange4"))(100),
        clustering_method="average",
        fontsize=7, fontsize_row=9, fontsize_col=9,
        file=paste0("./OUTPUT/figure_heatmap_", curr_signature, "_for_paper.pdf"),
        scale="row"
    )
}
dev.off()

rpkm_df = NULL
pathway[["AP1s"]] = c("FOS", "FOSL1", "FOSL2", "FOSB", "JUN", "JUNB", "JUND") 
signatures = names(pathway)
for(k in 1:length(signatures)){
    curr_signature = signatures[k]
    expr2 = expr[row.names(expr) %in% pathway[[curr_signature]], colnames(expr) %in% row.names(mapping)]
    dim(expr2)
    gene_std = apply(na.omit(expr2), 1, sd)
    expr2 = expr2[!rownames(expr2) %in% names(which(gene_std == 0)),]

    tmp_df = data.frame(expr2)
    tmp_df$signature = curr_signature
    tmp_df$symbol = row.names(tmp_df)
    row.names(tmp_df) = seq(1, nrow(tmp_df), 1)
    if(is.null(rpkm_df)){
        rpkm_df = tmp_df
    }else{
        rpkm_df = rbind(rpkm_df, tmp_df)
    }
}
mapping = data.frame(colData(se)[,c("TEST_ARTICLE", "CELL_LINE")])
df = reshape2::melt(rpkm_df)
df = merge(df, mapping, by.x="variable", by.y="row.names")
df = df %>%
    dplyr::group_by(signature, symbol, TEST_ARTICLE, CELL_LINE) %>%
    dplyr::summarise_at(vars(value), list(Mean = ~mean(.),
                                          Stdev = ~sd(.))) %>% as.data.frame() %>% as.data.frame()
df$TREATMENT = paste0(df$CELL_LINE, " + ", df$TEST_ARTICLE)

p = ggplot(df, aes(x=TREATMENT, y=Mean)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
    facet_grid(. ~ signature) +
    theme_minimal() +
    ylab("average RPKM") +
    labs(caption="Each point = average (of RPKMs) of 3 replicates of one gene belonging to its gene signature set.\nEach panel = one gene signature set.") +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(angle=0, size=9),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=9),
        plot.title=element_text(size=8),
        plot.caption=element_text(hjust=0, size=7)
    )
p
ggsave(paste0(outdir,"/RPKM_boxplot_gene_signatures_for_paper.pdf"), p, width=11.0, height=5.5, units="in", device="pdf")

df = reshape2::melt(rpkm_df)
df = merge(df, mapping, by.x="variable", by.y="row.names")
df = df[df$signature == "AP1s",]
df = df %>%
    dplyr::group_by(signature, symbol, TEST_ARTICLE, CELL_LINE) %>%
    dplyr::summarise_at(vars(value), list(Mean = ~mean(.),
                                          Stdev = ~sd(.))) %>% as.data.frame() %>% as.data.frame()
df$TREATMENT = paste0(df$CELL_LINE, " + ", df$TEST_ARTICLE)

color_list = list(
 "MSTO_211H + GNE-7883" = "gray",
 "MSTO_211H + DMSO" = "red",
 "MSTO_211H-7883R + GNE-7883" = "green",
 "MSTO_211H-7883R + DMSO" = "blue")

p = ggplot(df, aes(x=symbol, y=Mean, fill=TREATMENT)) +
    geom_bar(stat="identity", position=position_dodge(0.9), color="black") +
    geom_errorbar(aes(x=symbol, ymin=Mean-Stdev, ymax=Mean+Stdev), size=0.5,
                  width=.25,position=position_dodge(.9)) +
    scale_fill_manual(values=color_list) + 
    facet_grid(. ~ signature) +
    theme_minimal() +
    ylab("average RPKM") +
    labs(caption="Each point = average (of RPKMs) of 3 replicates of one gene belonging to its gene signature set.\nEach panel = one gene signature set.") +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(angle=0, size=9),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8),
        plot.caption=element_text(hjust=0, size=7)
    )
p
ggsave(paste0(outdir,"/RPKM_boxplot_gene_signatures_for_paper_AP1s.pdf"), p, width=7, height=4.0, units="in", device="pdf")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# AP1s 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
df = reshape2::melt(rpkm_df)
df = merge(df, mapping, by.x="variable", by.y="row.names")
df = df[df$signature == "additional_1",]
df = df %>%
    dplyr::group_by(signature, symbol, TEST_ARTICLE, CELL_LINE) %>%
    dplyr::summarise_at(vars(value), list(Mean = ~mean(.),
                                          Stdev = ~sd(.))) %>% as.data.frame() %>% as.data.frame()
df$TREATMENT = paste0(df$CELL_LINE, " + ", df$TEST_ARTICLE)

p = ggplot(df, aes(x=TREATMENT, y=Mean)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(shape=16, position=position_jitter(0.2), size=0.5) +
    facet_grid(. ~ symbol) +
    theme_minimal() +
    ylab("average RPKM") +
    labs(caption="Each point = average (of RPKMs) of 3 replicates of one gene belonging to its gene signature set.\nEach panel = one gene signature set.") +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(angle=0, size=7),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8),
        plot.caption=element_text(hjust=0, size=7)
    )
p
ggsave(paste0(outdir,"/RPKM_boxplot_gene_signatures_for_paper_cd44_egfr_vim.pdf"), p, width=5.0, height=5.5, units="in", device="pdf")



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# GSEA with filtering                
#                                    
# H226_G7883 vs H226_DMSO            
# G7883(CL3)_G7883 H226_G7883        
# H226-P-DMSO vs H226-7883R CL3 7883 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom_output = readRDS(paste0(outdir, "voom.output_JT.rds"))
voom_output_selected = voom_output

msigdb_results = list()
msigplot = list()
msigdb_results_df = NULL
cells = unique(se$CELLTYPE)
x = 0
qvalue_cutoff = 0.05
local(for(i in 1:length(cells)){
    cell = cells[i]
    cell2 = NULL
    if(cell == "MSTO_211H"){
        cell2 = "parental"
    }else{
        cell2 = "resistant"
    }
    #Hallmark
    category=c("H","C2", "C5", "C6","C8")

    for (j in 1:length(category)) {
        if(j ==1){
            m_t2g = msigdbr(species = "Homo sapiens", category = "H") %>%
                dplyr::select(gs_name, gene_symbol)
        }else{
            m_t2g = msigdbr(species = "Homo sapiens", category = category[[j]]) %>%
                dplyr::select(gs_name, gene_symbol)
        }

        for(k in 1:length(voom_output_selected[[cell2]])) {
            data1 = voom_output_selected[[cell2]][[k]]$LogFC
            names(data1) = rowData(se)$symbol
            data1 = sort(data1, decreasing = TRUE)
            data1_clean = data1[which(!is.na(names(data1)))]
            message(paste("analyzing", category[j],  names(voom_output_selected[[cell2]][k]), cell2))

            msigdb = GSEA(data1_clean, TERM2GENE=m_t2g)
            # NES: (normalized enrichment score)
            if(nrow(msigdb@result) > 0){
                curr_results = msigdb@result[order(msigdb@result$NES), ]
                curr_results$core_enrichment = NULL
                curr_results$direction = ifelse(curr_results$NES > 0, "UP", "DOWN")
                curr_results$cell = cell
                curr_results$category = category[j]
                curr_results$comparison = names(voom_output_selected[[cell2]][k])
                curr_results = curr_results[curr_results$qvalue < qvalue_cutoff,] #0.05
                msigdb_results[[ cell ]][[names(voom_output_selected[[cell2]][k])]][[ category[j] ]][[ k ]] = curr_results
            }else{
                message("...no term enriched under specific pvalueCutoff for ", paste0(cell, " - ",  names(voom_output_selected[[cell2]][k])))
            }

            if(x == 0){
                msigdb_results_df <<- curr_results
            }else{
                msigdb_results_df <<- rbind(msigdb_results_df, curr_results)
            }
            x = x + 1
        }
    }
})

# Write table
msigdb_results_df$comparison = gsub("\\`", "", msigdb_results_df$comparison)
msigdb_results_df$comparison = gsub("\\n", " ", msigdb_results_df$comparison)
msigdb_results_df2 = msigdb_results_df
msigdb_results_df2$comparison = gsub("\n", " ", msigdb_results_df2$comparison)
write.table(msigdb_results_df2, "./OUTPUT/GSEA_table.tsv", sep="\t", row.names=F, quote=F)

# Then process the resulting df to generate plots.
head(msigdb_results_df)
msigdb_results_df$comparison = gsub("Difference between ", "", msigdb_results_df$comparison)
msigdb_results_df$comparison = gsub("`", "", msigdb_results_df$comparison)
msigdb_results_df$comparison = gsub(" vs ", "\nvs\n", msigdb_results_df$comparison)
msigdb_results_df$comparison = factor(msigdb_results_df$comparison, levels=unique(msigdb_results_df$comparison))

# order by NES
curr_msigdb_results_df_all = NULL
for(curr_category in c("H", "C2", "C5", "C6", "C8")){
    curr_msigdb_results_df = msigdb_results_df[msigdb_results_df$category == curr_category,]
    order_df = curr_msigdb_results_df %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise_at(vars(NES), list(SUM = ~ sum(.),
                                            SD = ~sd(.))) %>%
        as.data.frame()
    order_df = order_df[order(-order_df$SUM),]

    curr_msigdb_results_df$ID = factor(curr_msigdb_results_df$ID, levels=unique(order_df$ID))
    IDs_vector = unique(curr_msigdb_results_df$ID)
    # Store in df for further use downstream.
    if(is.null(curr_msigdb_results_df_all)){
        curr_msigdb_results_df_all = curr_msigdb_results_df
    }else{
        curr_msigdb_results_df_all = rbind(curr_msigdb_results_df_all, curr_msigdb_results_df)
    }
    curr_msigdb_results_df$comparison2 = paste0(curr_msigdb_results_df$cell, "\n", curr_msigdb_results_df$comparison)
    curr_p = ggplot(curr_msigdb_results_df, aes(x=NES, y=ID, fill=qvalue)) +
        geom_bar(stat="identity") +
        facet_grid(. ~ comparison2, scales="free_y", space="free_y") +
        geom_vline(xintercept=0, linetype="dashed",  color="red", size=0.4) +
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
    curr_height = generate_plot_dimensions(curr_msigdb_results_df, by_variable="cell")[[1]]
    curr_width = generate_plot_dimensions(curr_msigdb_results_df, by_variable="cell")[[2]]

    ggsave(paste0(outdir,"/GSEA_", curr_category, "_for_paper2.pdf"), curr_p, width=curr_width, height=curr_height,
           units="in", device="pdf", limitsize=FALSE)

    df = curr_msigdb_results_df[curr_msigdb_results_df$comparison == "G7883\nvs\nDMSO",]
    signatures = as.character(unique(df$ID))
    signatures = signatures[!signatures %in% c("HALLMARK_G2M_CHECKPOINT",
                                               "HALLMARK_ESTROGEN_RESPONSE_LATE",
                                               "HALLMARK_BILE_ACID_METABOLISM",
                                               "HALLMARK_COMPLEMENT",
                                               "HALLMARK_HEME_METABOLISM",
                                               "HALLMARK_XENOBIOTIC_METABOLISM",
                                               "HALLMARK_COAGULATION"
                                               )]
    m_t2g = msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, gene_symbol)
    m_t2g = m_t2g[m_t2g$gs_name %in% signatures,]
    m_t2g = unique(m_t2g)
    write.table(m_t2g, paste0(outdir, "gene_tables_for_paper_fig1E_and_2G.tsv"), sep="\t", quote=F, row.names=F)


}

selected_gseas = c("OXIDATIVE_PHOSPHORYLATION",
                   "GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",
                   "ORGANELLE_INNER_MEMBRANE",
                   "RESPONSE_TO_CYTOKINE",
                   "DEFENSE_RESPONSE",
                   "CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
                   "ADAPTATIVE_IMMUNE_RESPONSE", # ADAPTIVE_IMMUNE_RESPONSE
                   "DEFENSE_RESPONSE_TO_OTHER_ORGANISM",
                   "INNATE_IMMUNE_RESPONSE",
                   "OXIDANT",
                   "WNT")

curr_msigdb_results_df_all_selected = curr_msigdb_results_df_all[grepl(paste(selected_gseas, collapse="|"), curr_msigdb_results_df_all$Description),]
curr_msigdb_results_df_all_selected$comparison = gsub("Difference between `", "", curr_msigdb_results_df_all_selected$comparison)
curr_msigdb_results_df_all_selected$comparison = gsub("`", "", curr_msigdb_results_df_all_selected$comparison)
curr_msigdb_results_df_all_selected$comparison2 = paste0(curr_msigdb_results_df_all_selected$cell, "\n", curr_msigdb_results_df_all_selected$comparison)
curr_p = ggplot(curr_msigdb_results_df_all_selected, aes(x=NES, y=ID, fill=qvalue)) +
    geom_bar(stat="identity") +
    facet_grid(. ~ comparison2, scales="free_y", space="free_y") +
    geom_vline(xintercept=0, linetype="dashed",  color="red", size=0.4) +
    ylab("") +
    ggtitle(paste0("GSEA analysis; selected functions. Category:", "[H] and [C5]","\nNES:Normalized Enrichment Score")) +
    theme_minimal() + theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle=90),
        axis.text.y=element_text(angle=0, size=6.5),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8)
    )
curr_p

ggsave(paste0(outdir,"/GSEA_", "selection_for_paper", ".pdf"), curr_p, width=7, height=4.34, units="in", device="pdf")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# up and down genes. Venn diagrams.
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom_df = readRDS(paste0(outdir,"voom.output_JT.rds"))

cell_lines = unique(colData(se)$CELLTYPE)
df_out = NULL
deg_results = list()
deg_gene_names_list_up = list() # for venn diagram.
deg_gene_names_list_down = list()
x = 1
fdr_cutoff = 0.05
logFC_cutoff = 1.0
for(i in 1:length(cell_lines)){
    curr_cell_line =  cell_lines[i]
    df = voom_df[[ curr_cell_line ]]

    comparisons = names(df)
    comparisons = c("Difference between `7883` vs `DMSO`")
    for(j in 1:length(comparisons)){
        curr_comparison = comparisons[j]
        curr_df = data.frame(df[[ curr_comparison ]])
        curr_df$symbol = rowData(se)$symbol
        # add annotations.
        curr_df = curr_df[!is.na(curr_df$LogFC),]
        curr_df = curr_df[!is.na(curr_df$FDR),]
        curr_df = curr_df[curr_df$FDR < fdr_cutoff, ]
        curr_df = curr_df[abs(curr_df$LogFC) >= logFC_cutoff, ]
        dim(curr_df)

        # split between up and down
        curr_df$direction = ifelse(curr_df$LogFC > 0, "up", "down")
        up = nrow(curr_df[curr_df$direction == "up", ])
        down = nrow(curr_df[curr_df$direction == "down", ])
        #store in list for later interrogation
        deg_results[[curr_cell_line]][[curr_comparison]] = curr_df
        deg_gene_names_list_up[[curr_cell_line]][[curr_comparison]] = curr_df[curr_df$direction == "up",]$symbol
        deg_gene_names_list_down[[curr_cell_line]][[curr_comparison]] = curr_df[curr_df$direction == "down",]$symbol

        tmp_df = data.frame(cell_line=curr_cell_line, comparison=curr_comparison, up=up, down=down)
        if(x == 1){
            df_out = tmp_df
        }else{
            df_out = rbind(df_out, tmp_df)
        }
        x = x + 1
    }
}

## Venn diagram
venn_diagrams_list = list()
for(cell_line in names(deg_gene_names_list_up)){
    all_genes = unique(as.vector(unlist(deg_gene_names_list_up[[cell_line]])))

    # UP
    venn_df = data.frame(row.names=all_genes)
    venn_df$symbol = row.names(venn_df)
    for(comp in names(deg_gene_names_list_up[[cell_line]])){
        tmp_df = data.frame(deg_gene_names_list_up[[cell_line]][[comp]])
        colnames(tmp_df)[1] = "V1"
        tmp_df$value = TRUE
        venn_df = merge(venn_df, tmp_df, by.x="symbol", by.y="V1", all=TRUE)
        colnames(venn_df)[ncol(venn_df)] = comp
    }
    venn_df[is.na(venn_df)] = FALSE
    venn_df$symbol = NULL
    colnames(venn_df) = gsub("Difference between " ,"",  colnames(venn_df))
    colnames(venn_df) = gsub("`" ,"",  colnames(venn_df))
    colnames(venn_df) = gsub(" vs " ,"\nvs\n",  colnames(venn_df))
    fit = euler(venn_df)
    print(plot(fit, quantities=TRUE, main=paste0(cell_line, "; Up"))    )
    venn_diagrams_list[["up"]] = plot(fit, quantities=TRUE, main=paste0(cell_line, "; Up; logFC >= ", logFC_cutoff, "; FDR < ", fdr_cutoff))

    # DOWN
    venn_df = data.frame(row.names=all_genes)
    venn_df$symbol = row.names(venn_df)
    for(comp in names(deg_gene_names_list_down[[cell_line]])){
        tmp_df = data.frame(deg_gene_names_list_down[[cell_line]][[comp]])
        colnames(tmp_df)[1] = "V1"
        tmp_df$value = TRUE
        venn_df = merge(venn_df, tmp_df, by.x="symbol", by.y="V1", all=TRUE)
        colnames(venn_df)[ncol(venn_df)] = comp
    }
    venn_df[is.na(venn_df)] = FALSE
    venn_df$symbol = NULL
    colnames(venn_df) = gsub("Difference between " ,"",  colnames(venn_df))
    colnames(venn_df) = gsub("`" ,"",  colnames(venn_df))
    colnames(venn_df) = gsub(" vs " ,"\nvs\n",  colnames(venn_df))
    fit = euler(venn_df)
    print(plot(fit, quantities=TRUE, main=paste0(cell_line, "; Down")) )
    venn_diagrams_list[["down"]] = plot(fit, quantities=TRUE, main=paste0(cell_line, "; Up; logFC >= ", logFC_cutoff, "; FDR < ", fdr_cutoff))
}
pdf(file=paste0("./OUTPUT/venn_diagram_for_paper.pdf"), width=8.58, height=9)
print(venn_diagrams_list[["up"]])
print(venn_diagrams_list[["down"]])
dev.off()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# boxplot of RPKMs for gene signatures                                         
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
expr2 = expr[row.names(expr) %in% pathway[["MAPK"]], ]
df = reshape2::melt(expr2)
colnames(df) = c("symbol", "Sample", "value")
df = merge(df, colData(se)[,c("SAMPLE_ID", "TREATMENT")], by.x="Sample", by.y="SAMPLE_ID")
order = as.data.frame(df) %>%
    dplyr::select(symbol, TREATMENT, value) %>%
    dplyr::group_by(symbol, TREATMENT) %>%
    dplyr::summarise_at(vars(value), list(Mean = ~mean(.),
                                          Stdev = ~sd(.))) %>% as.data.frame()
order
order2 = dcast(order, TREATMENT ~ symbol, value.var="Mean")
order2$mean = rowMeans(order2[2:ncol(order2)])
order3 = dcast(order, TREATMENT ~ symbol, value.var="Stdev")
order3$mean_stdev = rowMeans(order3[2:ncol(order3)])

expr2 = expr[row.names(expr) %in% pathway[["MAPK"]], ]
mapping = data.frame(colData(se)[,c("TREATMENT", "CELLTYPE")])
df = reshape2::melt(expr2)
colnames(df) = c("symbol", "SampleID", "value")
df = merge(df, mapping, by.x="SampleID", by.y="row.names")
df$TREATMENT = factor(df$TREATMENT, levels=c("MSTO_211H + DMSO", "MSTO_211H + GNE-7883", "MSTO_211H-7883R + DMSO", "MSTO_211H-7883R + GNE-7883"))
p = ggplot(df, aes(x=TREATMENT, y=value)) +
    geom_boxplot() +
    facet_grid(. ~ symbol) +
    theme_minimal() +
    theme(
        panel.border=element_rect(linewidth=0.5),
        panel.background=element_rect(fill='transparent', color=NA),
        rect = element_rect(fill = "transparent"),
        axis.text.x=element_text(angle=90),
        axis.text.y=element_text(angle=0, size=7),
        strip.text.y=element_text(angle=0, hjust=0),
        strip.text.x=element_text(angle=0, hjust=0.5, size=8),
        plot.title=element_text(size=8)
    )
p
pdf(file=paste0("./OUTPUT/RPMK_boxplots_MAPK.pdf"), width=8.58, height=5)
print(p)
dev.off()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check specific gene signatures                     
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
hippo = getFeatureSetCollection("GMTY188:analysis/hippo.gmt.bz2@REVISION-2")
names(hippo) = hippo@elementMetadata@listData[["name"]]
search.results <- searchFiles("hallmark")
hallmark <- genomitory::getFeatureSetCollection(search.results[2,"id"])
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

hippo2 = hippo[c("hippo145", "apoptosis", "mapk")]
df = stack(hippo2)
df = unique(df)
write.table(df, "./OUTPUT/gene_tables_for_MSTO_quantseq.tsv", sep="\t", quote=F, row.names=F)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Generate df to have list of up/down genes per      
# condition                                          
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom_df = readRDS(paste0(outdir,"voom.output.rds"))

cell_lines = unique(colData(se)$CELLTYPE)
df_out = NULL
deg_results = list()
deg_gene_names_list_up = list() # for venn diagram.
deg_gene_names_list_down = list()
fdr_cutoff = 0.05
logFC_cutoff = 1.0
final_df = NULL
for(i in 1:length(cell_lines)){
    curr_cell_line =  cell_lines[i]
    df = voom_df[[ curr_cell_line ]]

    comparisons = names(df)

    comparisons = c(
        "Difference between `7883` vs `DMSO`"
    )


    for(j in 1:length(comparisons)){
        curr_comparison = comparisons[j]
        curr_df = data.frame(df[[ curr_comparison ]])
        curr_df$symbol = rowData(se)$symbol
        # add annotations.
        curr_df = curr_df[!is.na(curr_df$LogFC),]
        curr_df = curr_df[!is.na(curr_df$FDR),]
        curr_df = curr_df[curr_df$FDR < fdr_cutoff, ]
        curr_df = curr_df[abs(curr_df$LogFC) >= logFC_cutoff, ]
        curr_df = curr_df[order(-curr_df$LogFC),]
        curr_df = curr_df[, c("symbol", "LogFC")]
        curr_df$comparison = curr_comparison
        curr_df$comparison = gsub("Difference between ", "", curr_df$comparison)
        curr_df$comparison = gsub("`", "", curr_df$comparison)

        curr_df$symbol = ifelse(curr_df$symbol == "SPATA13", paste0(curr_df$symbol, "_", row.names(curr_df)), curr_df$symbol)
        curr_df$CELLTYPE = curr_cell_line

        dim(curr_df)
        head(curr_df)

        if(is.null(final_df)){
            final_df = curr_df
        }else{
            final_df = rbind(final_df, curr_df)
        }
    }
}
final_df2 = reshape2::dcast(final_df, "symbol ~ comparison + CELLTYPE", value.var="LogFC", fun.aggregate=NULL, fill=0)
write.table(final_df2, "./OUTPUT/logFC_summary_foreach_gene_per_comparison.tsv", sep="\t", quote=F)
