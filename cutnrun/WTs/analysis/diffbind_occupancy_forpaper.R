library(DiffBind)
library(rtracklayer)
library(ChIPpeakAnno)
library(dplyr)
library(pheatmap)
library(vegan)
library(ChIPQC)
library(ggplot2)

ggsave <- function(..., bg='white') ggplot2::ggsave(..., bg=bg)

# This code is solely intented provide a reference on the processing of the CUTNRUN data in reference to the results of Figure 4

setwd("/path/to/my/manuscript/cutnrun/WTs/analysis")
source("./utils.R")
sampleSheet <- read.table("../OUTPUT/samplesheet_WTs.csv", sep="\t", header=T, check.names=F)
colnames(sampleSheet)[5] <- "Replicate"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# dba (load samples)                                                           
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
results <- dba(sampleSheet=sampleSheet)
results

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# dba (blacklist).                                                             
#                                                                              
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
results <- dba.blacklist(results, blacklist=DBA_BLACKLIST_GRCH38, greylist=FALSE, cores=1)
head(results$peaks[[1]])
saveRDS(results, file="../OUTPUT/dba.rds")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# dba (count).                                                                 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
results$config$singleEnd <- TRUE
results <- dba.count(results, summits=250,  bParallel=3)
saveRDS(results, "../OUTPUT/dba_count.rds")

p_heatmap1 <- dba.plotHeatmap(results, attributes=DBA_TREATMENT, ColAttributes=c(DBA_FACTOR, DBA_CONDITION), colSideCols=list(c("orange","orchid1","mediumpurple1", "blue"), c("magenta", "cyan", "red")))
p_heatmap2 <- dba.plotHeatmap(results, attributes=DBA_FACTOR, ColAttributes=c(DBA_FACTOR, DBA_CONDITION), colSideCols=list(c("orange","orchid1","mediumpurple1", "blue"), c("magenta", "cyan", "red")))
p_pca1 <- dba.plotPCA(results, attributes=DBA_FACTOR, label=DBA_TREATMENT)
p_pca2 <- dba.plotPCA(results, attributes=DBA_TREATMENT, label=DBA_TREATMENT)

# Normalize background regions using csaw gives best differential
results_norm <- dba.normalize(results,
                         library=DBA_LIBSIZE_BACKGROUND,
                         method=DBA_ALL_METHODS,
                         normalize=DBA_NORM_NATIVE,
                         background=TRUE)
saveRDS(results_norm, "../OUTPUT/dba_count_norm.rds")

# write clean version of plots:
binding_df <- data.frame(results_norm$binding)
row.names(binding_df) <- paste0(binding_df$CHR, ":", binding_df$START, "-", binding_df$END)
binding_df$CHR <- NULL
binding_df$START <- NULL
binding_df$END <- NULL

binding_df_cor <- cor((binding_df))
binding_df_cor2 <- data.frame(data.matrix(binding_df_cor))
mapping3 <- sampleSheet[,c(1,2,3,4,5)]
row.names(mapping3) <- sampleSheet$SampleID
mapping3$sampleName <- paste0(mapping3$Treatment, "_", mapping3$Factor, "_Rep", mapping3$Rep)
samplenames_map <- mapping3[,c("sampleName")]
names(samplenames_map) <- row.names(mapping3)
row.names(binding_df_cor2 ) <- samplenames_map[row.names(binding_df_cor2)]
colnames(binding_df_cor2 ) <- samplenames_map[colnames(binding_df_cor2)]

curr_list <- list()
mapping4 <- mapping3[,c(2,3,4,5)]
row.names(mapping4) <- mapping3$sampleName; mapping4$sampleName=NULL;
mapping4$Rep <- as.character(mapping4$Rep)
x=1
for(j in 1:ncol(mapping4)){
    curr_col_name <- names(mapping4)[j]
    curr_var_names <- unique(mapping4[,j])
    curr_colors <- vColors[x:(x+(length(curr_var_names)-1))]
    names(curr_colors) <- curr_var_names
    x <- length(curr_var_names) + 1 + x
    curr_list[[curr_col_name]] <- curr_colors
}
print(curr_list)
dir.create("../OUTPUT/diffbind/")
pheatmap::pheatmap(
    binding_df_cor2,
    main=paste0("Heatmap of correlation between samples by normalized binding sites scores"),
    filename=paste0("../OUTPUT/diffbind/dba_norm_heatmap.png"),
    show_rownames=TRUE, show_colnames=TRUE,
    border_color="gray",
    fontsize=9, fontsize_row=9,fontsize_col=9,
    cellwidth=12,
    cellheight=12,
    annotation=mapping4,
    annotation_colors=curr_list,
    color=colorRampPalette(c("white", "green4"))(100),
    clustering_method="average",
    clustering_distance_rows="correlation",
    clustering_distance_cols="correlation"
)
dev.off()

binding_df_cor <- cor((binding_df))
binding_df_cor2 <- data.frame(data.matrix(binding_df_cor))
mapping <- sampleSheet[,c(1,2,3,4,5)]
row.names(mapping) <- mapping$SampleID; mapping$SampleID <- NULL;
pcoa <- prcomp(binding_df_cor2)
summary <- data.frame(summary(pcoa)[[6]])
perc1 <- round(summary$PC1[2] * 100, digits=2)
perc2 <- round(summary$PC2[2] * 100, digits=2)
pcoa <- pcoa$rotation[,c(1,2)]
pcoa <- merge(pcoa, mapping, by.x="row.names", by.y="row.names")
if(any(colnames(pcoa) == "Rep") == FALSE){ pcoa$Rep <- pcoa$Replicate}
pcoa$Replicate <- as.character(pcoa$Replicate)

my_colors <- vColors[1:length(unique(pcoa$Treatment))]
p_pcoa <- ggplot(pcoa, aes(x=PC1, y=PC2, shape=Condition, color=Replicate, fill=Replicate)) +
    geom_point(size=4) +
    scale_shape_manual(values=c(21,22,23,24,25)) +
    scale_color_manual(values=vColors) +
    scale_fill_manual(values=vColors) +
    xlab(paste0("PC1 (", perc1, "%)")) + ylab(paste0("PC2 (", perc2, "%)")) +
    theme_minimal() +
    theme(panel.border=element_rect(linewidth=0.5),
          panel.background=element_rect(fill='transparent', color=NA),
          rect=element_rect(fill="transparent"),
          axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5),
          plot.title=element_text(size=10, hjust=0.5), legend.position="right") +
    ggtitle(paste("PCoA of correlations of peaks; split by replicates."))
p_pcoa
ggsave("../OUTPUT/diffbind/pcoa_norm.pdf", plot=p_pcoa, device="pdf", units="in", width=4.38, height=4, dpi=300)
ggsave("../OUTPUT/diffbind/pcoa_norm.png", plot=p_pcoa, device="png", units="in", width=4.38, height=4, dpi=300)

# check lib size factors
results_norm$norm$DESeq2$norm.facs
results_norm$norm$edgeR$norm.facs

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# annotate peaks                                                               
#                                                                              
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v86)
ensembl.hs86.transcript <- transcripts(EnsDb.Hsapiens.v86)
macs_peak_gr <- toGRanges(sampleSheet[1,]$Peaks, format="MACS2")
head(macs_peak_gr, n=2)

macs_peak_ensembl <- annotatePeakInBatch(macs_peak_gr, AnnotationData=ensembl.hs86.transcript)
head(macs_peak_ensembl)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Peaks feature distrib 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks_anno_df <- NULL
for(i in 1:nrow(sampleSheet)){
    peak_file <- sampleSheet[i,]$Peaks
    tmp <- tryCatch(read.table(peak_file), error=function(e) NULL)
    if(is.null(tmp)){ next }
    macs_peak_gr <- toGRanges(peak_file, format="MACS2")
    
    res <- annotatePeak(macs_peak_gr, TxDb=txdb, verbose=F)
    p <- plotAnnoBar(res)
    tmp_df <- ggplot_build(p)$plot$data
    tmp_df$factor <- sampleSheet[i,]$Factor
    tmp_df$condition <- sampleSheet[i,]$Condition
    tmp_df$treatment <- sampleSheet[i,]$Treatment
    tmp_df$rep <- sampleSheet[i,]$Rep
    
    if(is.null(peaks_anno_df)){
        peaks_anno_df <- tmp_df
    }else{
        peaks_anno_df <- rbind(peaks_anno_df, tmp_df)
    }
}

peaks_anno_df$condition <- factor(peaks_anno_df$condition, levels=unique(peaks_anno_df$condition))
p <- ggplot(peaks_anno_df, aes(y=Frequency, x=factor, fill=Feature)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values=vColors) +
    facet_grid(rep ~ condition) +
    theme_minimal() +
    theme(
        panel.border=element_rect(fill=NA, linetype="solid", colour="black", linewidth=0.25),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90, hjust=0)
    )
p
ggsave("../OUTPUT/peaks_anno.pdf", plot=p, device="pdf", units="in", width=5.5, height=5.6, dpi=300)
ggsave("../OUTPUT/peaks_anno.png", plot=p, device="png", units="in", width=5.5, height=5.6, dpi=300)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Differential analysis                                                        #
#                                                                              #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
dir.create("../OUTPUT/diffbind/")

# Perform contrast
contrast_list <- list(
    c("7883RCL3_G7883",  "P_G7883"),
    c("7883RCL3_G7883",  "P_DMSO"),
    c("P_G7883",         "P_DMSO")
)

for(comparison in contrast_list){
    message("comparing ", comparison[1], " vs ", comparison[2])
    comparison_name <-  paste(c(comparison[1], comparison[2]), collapse="_vs_")

    for(factor in unique(sampleSheet$Factor)){
        results_block <- dba(results_norm, mask=results_norm$masks[[factor]])
        results_block[results_block$samples$Replicate]
        
        contrast <- dba.contrast(results_block,
                                design="~Condition+Replicate",
                                contrast=c("Condition", comparison[1], comparison[2]))
        contrast <- dba.analyze(contrast,
                               bBlacklist=TRUE,
                               bGreylist=FALSE,
                               method=DBA_ALL_METHODS)

        filename <- paste0("../OUTPUT/diffbind/contrast_", factor, "_", comparison_name, ".rds")
        saveRDS(contrast, filename)

        for(diff_method in c("DBA_EDGER","DBA_DESEQ2")){

            # Plot MA plot
            pdf(paste0("../OUTPUT/diffbind/MAplot_", diff_method, "_", factor, "_", comparison_name, ".pdf"))
            if (diff_method == "DBA_EDGER"){
                dba.plotMA(contrast, th=0.1, yrange=c(-10,10), method=DBA_EDGER)

            }else if(diff_method == "DBA_DESEQ2"){
                dba.plotMA(contrast, th=0.1, yrange=c(-10,10), method=DBA_DESEQ2)
            }
            dev.off()
            
            png(paste0("../OUTPUT/diffbind/MAplot_", diff_method, "_", factor, "_", comparison_name, ".png"))
            if (diff_method == "DBA_EDGER"){
                dba.plotMA(contrast, th=0.1, yrange=c(-10,10), method=DBA_EDGER)
                
            }else if(diff_method == "DBA_DESEQ2"){
                dba.plotMA(contrast, th=0.1, yrange=c(-10,10), method=DBA_DESEQ2)
            }
            dev.off()

            # Generate report
            if(diff_method == "DBA_EDGER"){
                report <- dba.report(contrast, contrast=1, th=0.1, method=DBA_EDGER)
            }else if(diff_method =="DBA_DESEQ2"){
                report <- dba.report(contrast, contrast=1, th=0.1, method=DBA_DESEQ2)
            }

            if(!is.null(report)){
                report_gain <- report[report$Fold>0]
                report_lost <- report[report$Fold<0]
            }

            report_all <- dba.report(contrast, contrast=1, th=1)
            report_all$name <- paste0("peak_", 1:length(report_all))

            if(exists("report_gain")){
                if(length(report_gain) > 0){
                    message("exporting gained peaks")
                    report_gain$name <- paste0("peak_", 1:length(report_gain))
                    export.bed(report_gain, paste0("../OUTPUT/diffbind/diffbind_gain_", factor, "_", diff_method, "_", comparison_name, ".bed")
                    )
                }
            }

            if(exists("report_lost")) {
                if(length(report_lost) > 0 ){
                    message("exporting lost peaks")
                    report_lost$name <- paste0("peak_", 1:length(report_lost))
                    export.bed(report_lost, paste0("../OUTPUT/diffbind/diffbind_lost_", factor, "_", diff_method, "_", comparison_name, ".bed")
                    )
                }
            }

            if(exists("report_all")) {
                message("exporting all peaks")
                head(report_all)
                export.bed(report_all, paste0("../OUTPUT/diffbind/diffbind_all_", factor, "_", diff_method, "_", comparison_name, ".bed")
                )
            }
        }
    }
}

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Nearest genes to peaks
# chipenrich.           
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(chipenrich)
library(rtracklayer)
my_gains <- list.files("../OUTPUT/diffbind/", pattern="gain.*EDGER")
my_losses <- list.files("../OUTPUT/diffbind/", pattern="lost.*EDGER")

my_gains <- my_gains[grep("TEAD1|FOSL1|YAP1", my_gains)]
my_losses <- my_losses[grep("TEAD1|FOSL1|YAP1", my_losses)]

my_gains_losses <- c(my_gains, my_losses)

df_nearest_genes <- NULL
for(my_peaks_file in my_gains_losses){
    prefix <- gsub("^diffbind_.*DBA_EDGER_", "", my_peaks_file)
    prefix <- gsub(".bed", "", prefix)
    tf <- gsub("^.*(TEAD1|FOSL1|YAP1).*", "\\1", my_peaks_file)
    status <- gsub("^.*(gain|lost).*", "\\1", my_peaks_file)
    if(status == "lost"){
        status <- "loss"
    }
    
    my_peaks <- import(paste0("../OUTPUT/diffbind/", my_peaks_file))
    my_peaks <- data.frame(seqnames=seqnames(my_peaks),
                          start=start(my_peaks)-1,
                          end=end(my_peaks))
    res <- chipenrich(peaks=my_peaks, genome="hg38", genesets="hallmark", locusdef="nearest_tss", qc_plots=FALSE, out_name=NULL, n_cores=1)
    my_nearest_genes <- res$peaks
    my_nearest_genes$comparison <- prefix
    my_nearest_genes$direction <- status
    my_nearest_genes$chip_tf <- tf
    
    if(is.null(df_nearest_genes)){
        df_nearest_genes <- my_nearest_genes
    }else{
        df_nearest_genes <- rbind(df_nearest_genes, my_nearest_genes)
    }
}

head(df_nearest_genes)
tail(df_nearest_genes)
write.table(df_nearest_genes, "../OUTPUT/nearest_genes_NGS5325_NGS5516_NGS5722.tsv", sep="\t", quote=FALSE, row.names=FALSE)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Venn diagrams                                                                #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# proportional venn diagrams are nice when circles do not overlap too much. However, when they do, its better to use something like ggVennDiagrams
library(ggVennDiagram)

# Helper function that transform makeVennDiagram object into list compatible with ggVennDiagram
buildListForVennDiagram <- function(makeVennDiagram_object){
    venn_counts <- data.matrix(makeVennDiagram_object$vennCounts)
    idx <- match("Counts", colnames(venn_counts))
    my_list <- list()
    curr_value <- 1
    
    for(i in 1:nrow(venn_counts)){
        counts <- venn_counts[i,idx]
        if(counts == 0){
            next
        }
        
        for(j in 1:ncol(venn_counts[,1:(idx - 1)])){
            my_colname <- colnames(venn_counts)[j]
            flag <- venn_counts[i,j]
            if(flag[[1]] != 0){
                my_list[[my_colname]] <- c(my_list[[my_colname]], seq(from=curr_value, to=(curr_value + counts[[1]] - 1), by=1))
            }
        }
        curr_value <- (curr_value + counts[[1]] + 0)
    }
    return(my_list)
}

peaks_list <- list()
for(i in 1:nrow(sampleSheet)){
    factor <- sampleSheet[i,]$Factor
    condition <- sampleSheet[i,]$Condition
    replicate <- as.character(sampleSheet[i,]$Replicate)
    my_gr <- tryCatch({toGRanges(sampleSheet[i,]$Peaks, format="MACS2", header=FALSE)}, error = function(e) NULL)
    if(!is.null(my_gr)){
        peaks_list[[factor]][[condition]][[replicate]] <- my_gr
    }
}

my_peaks_list_3el <- list(
    "H3K27Ac_P_DMSO_reps" =                    list(c("H3K27Ac", "P_DMSO", 1),         c("H3K27Ac", "P_DMSO",  3),        c("H3K27Ac", "P_DMSO", 4)),
    "H3K27Ac_7883RCL3_G7883_vs_P_G7883_rep1" = list(c("H3K27Ac", "7883RCL3_G7883", 1), c("H3K27Ac", "P_G7883", 1),        c("H3K27Ac", "P_DMSO", 1)),
    "H3K27Ac_7883RCL3_G7883_vs_P_G7883_rep3" = list(c("H3K27Ac", "7883RCL3_G7883", 3), c("H3K27Ac", "P_G7883", 3),        c("H3K27Ac", "P_DMSO", 3)),
    
    "FOSL1_P_DMSO_reps" =                    list(c("FOSL1", "P_DMSO", 1),         c("FOSL1", "P_DMSO",  3),        c("FOSL1", "P_DMSO", 4)),
    "FOSL1_P_G7883_reps" =                   list(c("FOSL1", "P_G7883", 1),        c("FOSL1", "P_G7883", 3),        c("FOSL1", "P_G7883", 4)), 
    "FOSL1_7883RCL3_G7883_reps" =            list(c("FOSL1", "7883RCL3_G7883", 1), c("FOSL1", "7883RCL3_G7883", 3), c("FOSL1", "7883RCL3_G7883", 4)),
    "FOSL1_7883RCL3_G7883_vs_P_G7883_rep1" = list(c("FOSL1", "7883RCL3_G7883", 1), c("FOSL1", "P_G7883", 1),        c("FOSL1", "P_DMSO", 1)),
    "FOSL1_7883RCL3_G7883_vs_P_G7883_rep3" = list(c("FOSL1", "7883RCL3_G7883", 3), c("FOSL1", "P_G7883", 3),        c("FOSL1", "P_DMSO", 3)),
    "FOSL1_7883RCL3_G7883_vs_P_G7883_rep4" = list(c("FOSL1", "7883RCL3_G7883", 4), c("FOSL1", "P_G7883", 4),        c("FOSL1", "P_DMSO", 4))
)

my_peaks_list_2el <- list(
    "H3K27Ac_P_G7883_reps" =                   list(c("H3K27Ac", "P_G7883", 1),        c("H3K27Ac", "P_G7883", 3)), 
    "H3K27Ac_7883RCL3_G7883_reps" =            list(c("H3K27Ac", "P_G7883", 1),        c("H3K27Ac", "P_G7883", 3)), 
    "H3K27Ac_7883RCL3_G7883_vs_P_G7883_rep4" = list(c("H3K27Ac", "7883RCL3_G7883", 1),  c("H3K27Ac", "P_DMSO", 3))
)

# 3 elements
for(curr_name in names(my_peaks_list_3el)){
    print(curr_name)
    i <- my_peaks_list_3el[[curr_name]]
    sub_names <- c(paste(i[[1]][[1]], i[[1]][[2]], i[[1]][[3]], sep="_"),
                  paste(i[[2]][[1]], i[[2]][[2]], i[[2]][[3]], sep="_"),
                  paste(i[[3]][[1]], i[[3]][[2]], i[[3]][[3]], sep="_"))
    
    ol <- findOverlapsOfPeaks(peaks_list[[ i[[1]][[1]]  ]][[ i[[1]][[2]] ]][[ i[[1]][[3]] ]],
                             peaks_list[[ i[[2]][[1]]  ]][[ i[[2]][[2]] ]][[ i[[2]][[3]] ]],
                             peaks_list[[ i[[3]][[1]]  ]][[ i[[3]][[2]] ]][[ i[[3]][[3]] ]])
    vd <- makeVennDiagram(ol)
    colnames(vd$vennCounts)[1:length(sub_names)] <- sub_names
    colnames(vd$vennCounts)[(length(sub_names) + 2):length(colnames(vd$vennCounts))] <- sub_names
    
    p <- ggVennDiagram(buildListForVennDiagram(vd), category.names=gsub("_", "\n", names(buildListForVennDiagram(vd))), set_size=3, label_alpha=0) + ggplot2::scale_fill_distiller(palette="RdBu") + labs(title="Overlapping of peaks for: ", subtitle=paste(colnames(vd$vennCounts)[1:3], collapse="; "))  #scale_fill_gradient(low="#3884D4", high="red")
    ggsave(paste0("../OUTPUT/diffbind/vd_", curr_name, ".pdf"), plot=p, device="pdf", units="in", width=9, height=7.3, dpi=300)
    ggsave(paste0("../OUTPUT/diffbind/vd_", curr_name, ".png"), plot=p, device="png", units="in", width=9, height=7.3, dpi=300)
}

# 2 elements
for(curr_name in names(my_peaks_list_2el)){
    print(curr_name)
    i <- my_peaks_list_2el[[curr_name]]
    sub_names <- c(paste(i[[1]][[1]], i[[1]][[2]], i[[1]][[3]], sep="_"),
                  paste(i[[2]][[1]], i[[2]][[2]], i[[2]][[3]], sep="_"))
    
    ol <- findOverlapsOfPeaks(peaks_list[[ i[[1]][[1]]  ]][[ i[[1]][[2]] ]][[ i[[1]][[3]] ]],
                             peaks_list[[ i[[2]][[1]]  ]][[ i[[2]][[2]] ]][[ i[[2]][[3]] ]])
    vd <- makeVennDiagram(ol)
    colnames(vd$vennCounts)[1:length(sub_names)] <- sub_names
    colnames(vd$vennCounts)[(length(sub_names) + 2):length(colnames(vd$vennCounts))] <- sub_names
    
    p <- ggVennDiagram(buildListForVennDiagram(vd), category.names=gsub("_", "\n", names(buildListForVennDiagram(vd))), set_size=3, label_alpha=0) + ggplot2::scale_fill_distiller(palette="RdBu") + labs(title="Overlapping of peaks for: ", subtitle=paste(colnames(vd$vennCounts)[1:3], collapse="; "))  #scale_fill_gradient(low="#3884D4", high="red")
    ggsave(paste0("../OUTPUT/diffbind/vd_", curr_name, ".pdf"), plot=p, device="pdf", units="in", width=9, height=7.3, dpi=300)
    ggsave(paste0("../OUTPUT/diffbind/vd_", curr_name, ".png"), plot=p, device="png", units="in", width=9, height=7.3, dpi=300)
}

# END
