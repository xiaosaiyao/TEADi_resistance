# This code is solely intented provide a reference on the processing of the RNAseq data in reference to the results of Figure 2.

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
setwd("/path/to/my/manuscript/RNAseq/NGS4372/analysis/")
se <- getDatasetAsSE("DS000012520")
se$CELL_LINE <- sapply(strsplit(se$SAMPLE_LABEL, " "), "[", 2)
se$TREATMENT_NAME <- stringr::str_trim(gsub("SP_EXP1 ", "", se$TREATMENT_NAME))
se$TREATMENT <- paste0(se$CELL_LINE, "_", se$TREATMENT_NAME)
outdir <- "OUTPUT/"

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# run voom
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom.output=list()
for (cell in unique(se$CNAME)) {
  #set up comparisons
  contrast.list <- list(
                       c("DMSO(83)_G7883", "DMSO(83)_DMSO"),
                       c("G7883(CL3)_G7883", "G7883(CL3)_DMSO"),
                       c("H226_G7883", "H226_DMSO"),
                       c("G7883(CL3)_G7883","DMSO(83)_G7883"),
                       c("G7883(CL3)_DMSO","DMSO(83)_DMSO"),
                       c("G7883(CL3)_G7883","H226_G7883"),
                       c("G7883(CL3)_DMSO","H226_DMSO"),
                       c("G7883(CL3)_G7883", "H226_DMSO")
                   )


  voom.output[[cell]] <- runVoom(
    se,
    subset.factor="CNAME",
    subset.levels=cell,
    group="TREATMENT",
    comparisons=contrast.list
  )

  # write table
  filename <- make.names(names(voom.output[[cell]]))
  filename <- gsub("Difference.between.","", filename)

  for (n in 1:length(voom.output[[cell]])){
    df <- voom.output[[cell]][[n]]
    df$symbol <- rowData(se)$symbol
    write.table(df, col.names=NA, sep="\t",file=paste0(outdir, make.names(cell), filename[n],".txt"))
  }
}

saveRDS(voom.output, paste0(outdir,"voom.output.rds"))

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# perform GSEA and generate plots
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
msigdb_results <- list()
msigplot <- list()
for (cells in unique(se$CNAME)) {

  #Hallmark
  category=c("H","C2","C6","C8")
  subcategory=c(NA,"CGP","CGP", "CGP")

  for (j in 1:length(category)) {
    if (j ==1){
      m_t2g <- msigdbr(species="Homo sapiens", category="H") %>%
        dplyr::select(gs_name, gene_symbol)} else{
          m_t2g <- msigdbr(species="Homo sapiens", category=category[[j]]) %>%
            dplyr::select(gs_name, gene_symbol)
        }




    for (i in 1:length(voom.output[[cells]])) {
      data1 <- voom.output[[cells]][[i]]$LogFC
      names(data1) <- rowData(se)$symbol
      data1 <- sort(data1, decreasing <- TRUE)
      data1_clean <- data1[which(!is.na(names(data1)))]
      message(paste("analyzing", category[j],  names(voom.output[[cells]][i]), cells))


      msigdb <- GSEA(data1_clean, TERM2GENE=m_t2g)
      msigdb_results[[cells]][[category[j]]][[i]] <- msigdb@result[order(msigdb@result$NES), ]
      msigdb_results[[cells]][[category[j]]][[i]]  <-
        rbind(head(msigdb_results[[cells]][[category[j]]][[i]] , 10), tail(msigdb_results[[cells]][[category[j]]][[i]] , 10))
      msigdb_results[[cells]][[category[j]]][[i]]$ID <- factor(as.character(msigdb_results[[cells]][[category[j]]][[i]]$ID),
                                                              levels=unique(as.character(msigdb_results[[cells]][[category[j]]][[i]]$ID)))
      msigdb_results[[cells]][[category[j]]][[i]]$direction <- ifelse(msigdb_results[[cells]][[category[j]]][[i]]$NES > 0, "UP", "DOWN")
      msigdb_results[[cells]][[category[j]]][[i]]$direction <- factor(as.character(msigdb_results[[cells]][[category[j]]][[i]] $direction),
                                                                     levels=c("UP", "DOWN"))

      #plot results
      msigplot[[paste0(cells,".",category[j],".",names(voom.output[[cells]][i]))]]=plotMsigplot(msigdb_results[[cells]][[category[j]]][[i]] ,
                                                                                                paste(cells, category[j], names(voom.output[[cells]][i])), manual_color=c("red", "blue")) +scale_x_discrete(label=function(x) stringr::str_trunc(x, 50))
    }
  }
}
heatmaplist.m=marrangeGrob(grobs=msigplot, nrow=2, ncol=2)
ggsave(paste0(outdir,"GSEA.pdf"), heatmaplist.m, width=12, height=9)
saveRDS(msigdb_results, paste0(outdir, "msigdb_results.rds"))


