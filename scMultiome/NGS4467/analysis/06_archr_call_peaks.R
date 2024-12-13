#!/usr/bin/env Rscript

# Workflow that generate peaks
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

archr_call_peaks <- function(indir=NULL, num_threads=num_threads) {

    #setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
    #indir = "./output/archr/"
    #num_threads = 4
    #num_threads = as.numeric(num_threads)
    
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
    library(MASS)
    library(viridis)
    library(scales)
    library(BSgenome.Hsapiens.Genentech.GRCh38)
    library(ggpubr)
    library(chromVARmotifs)
    source("./tools/R/lib/utils.R")
    
    archr_proj = ArchR::loadArchRProject(indir)
    
    # Loading metadata
    mapping = getCellColData(archr_proj, select = c("hash_assignment2", "CNAME", "TEST_ARTICLE", "Treatment", "library", "HTO", "Treatment_HTO", "Treatment_library"))
    mapping$Treatment_library_HTO = paste0(mapping$Treatment_library, "_", mapping$HTO)
    archr_proj$Treatment_library_HTO = mapping$Treatment_library_HTO
    mapping = unique(mapping)
    mapping = data.frame(mapping)
    
    message("Executing addGroupCoverages()")
    archr_proj = addGroupCoverages(archr_proj, groupBy='Clusters_TileMatrix', threads=num_threads)
    message("Executing addRepdoduciblePeakSet()")
    archr_proj = addReproduciblePeakSet(archr_proj, groupBy="Clusters_TileMatrix", method="q", cutOff=0.05,
                                        pathToMacs2="/gstore/home/tremblj2/software/macs2/macs2_venv/bin/macs2", 
                                        excludeChr=c('chrMT','chrY'), genomeSize=2.7e9, threads=num_threads)
  
    peaks = getPeakSet(archr_proj)
    peaks = peaks@metadata[["PeakCallSummary"]]
    p1 = ggplot(peaks, aes(x=Group, y=Freq, fill=Var1)) +
        geom_bar(stat="identity") +
        facet_grid(. ~ Group, space="free_x", scale="free_x") +
        scale_fill_manual(values=vColors) +
        theme_minimal() + theme(panel.border=element_rect(fill=NA, linetype="solid", colour = "black", linewidth=0.5), axis.text.x=element_text(angle=90, hjust=1), strip.text.y=element_text(angle=0)) +
        ggtitle("Peak counts (MACS2)") + 
        ylab("Number of peaks (x10^3)") +
        labs(caption="Barplot showing the number of peaks per sample. Peaks were called with MACS2 using ArchR's addReproduciblePeakSet() function.")
    ggsave(paste0(indir, "/PeakCalls/tsne_peaks_Clusters_TileMatrix.png"), plot=p1, device="png", height=5.31, width=11.4, units="cm")
    ggsave(paste0(indir, "/PeakCalls/tsne_peaks_Clusters_TileMatrix.pdf"), plot=p1, device="pdf", height=5.31, width=11.4, units="cm")
    
    message("Executing addPeakMatrix()")
    archr_proj = addPeakMatrix(archr_proj, binarize=FALSE, threads=num_threads, force=TRUE)
    
    # TF annotation
    message("Annotating TFs")
    peaks_anno = genomitory::getFeatures('GMTY162:hg38_motif_bed_granges@REVISION-3')
    archr_proj = addPeakAnnotations(archr_proj, regions=peaks_anno, name='ENCODE_and_cistromeDB_TF_peaks', force=TRUE)
    archr_proj = addDeviationsMatrix(archr_proj, peakAnnotation='ENCODE_and_cistromeDB_TF_peaks', 
                                     matrixName='TFPeaksDeviationsMatrix', threads=num_threads, force=TRUE)
    
    # motif annotation
    message("Annotating motifs")
    archr_proj = addMotifAnnotations(archr_proj, motifSet='cisbp', annoName='Motif', species='Homo sapiens', force=TRUE)
    archr_proj = addDeviationsMatrix(archr_proj, peakAnnotation='Motif', threads=num_threads, force=TRUE)
  
    # TF binding matrix 
    message("TF_binind matrix")
    archr_proj = addPeakAnnotations(archr_proj, regions=epiregulon::getTFMotifInfo(genome="hg38"), name="TF_binding", force=T, logFile="addPeakAnnotations")
    archr_proj = addDeviationsMatrix(archr_proj, peakAnnotation="TF_binding", threads=num_threads)
    
    # Create bigwig files
    getGroupBW(archr_proj, groupBy = "Sample", normMethod="ReadsInTSS", tileSize=100, maxCells=1000, ceiling=4, verbose=TRUE, threads=4)
    
    
    saveArchRProject(archr_proj, overwrite=TRUE, logFile=createLogFile("saveArchRProject"), threads=num_threads)
    message("Completed archr_call_peaks.R!")
}

usage=function(errM) {
    cat("\nUsage : Rscript my_script_name.R [option] <Value>\n")
    cat("       -i        : indir containing the archR data structure.")
    cat("       -n        : num_threads.")
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
    }
}
archr_call_peaks(indir=indir, num_threads=num_threads)
