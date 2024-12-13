# This code is solely intented provide a reference on the processing of the Quantseq data.

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
setwd("/path/to/my/manuscript/Quantseq/DS000015488/")
dir.create("./OUTPUT")

se <- getDatasetAsSE("DS000015488")
se <- se[,se$cell_name == "MSTO-211H"]
colData(se)

se <- se[,grep("DMSO|7883", se$SAM_COMS)]

se$TEST_ARTICLE = "unknown"
se$TEST_ARTICLE = ifelse(grepl("DMSO", se$SAM_COMS), "DMSO", se$TEST_ARTICLE)
se$TEST_ARTICLE = ifelse(grepl("7883", se$SAM_COMS), "7883", se$TEST_ARTICLE)
se$CELLTYPE = ifelse(grepl("par", se$SAM_COMS), "parental", se$SAM_COMS)
se$CELLTYPE = ifelse(grepl("res", se$SAM_COMS), "resistant", se$CELLTYPE)
se$TREATMENT = paste0(se$CELLTYPE, "_", se$TEST_ARTICLE)
colData(se)
outdir="OUTPUT/"

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# run voom
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom.output=list()
for (cell in unique(se$CELLTYPE)) {
  #set up comparisons
  contrast.list = list(
      c("7883", "DMSO")
  )

  voom.output[[cell]] <- runVoom(
    se,
    subset.factor = "CELLTYPE",
    subset.levels = cell,
    group = "TEST_ARTICLE",
    comparisons = contrast.list,
    commit = "always"
  )

  # write table
  filename=make.names(names(voom.output[[cell]]))
  filename=gsub("Difference.between.","", filename)

  for (n in 1:length(voom.output[[cell]])){
    df=voom.output[[cell]][[n]]
    df$symbol=rowData(se)$symbol
    write.table(df, col.names = NA, sep="\t",file = paste0(outdir, make.names(cell), filename[n],".txt"))
  }
}

saveRDS(voom.output, paste0(outdir,"voom.output.rds"))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# perform GSEA and generate plots
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
voom.output = readRDS(file=paste0(outdir,"voom.output_JT.rds"))

msigdb_results = list()
msigplot = list()
for (cells in unique(se$CELLTYPE)) {

    #Hallmark
    category=c("H","C2","C6","C8")
    subcategory=c(NA,"CGP","CGP", "CGP")

    for (j in 1:length(category)) {
        if(j ==1){
            m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
            dplyr::select(gs_name, gene_symbol)} else{
              m_t2g <- msigdbr(species = "Homo sapiens", category = category[[j]]) %>%
                dplyr::select(gs_name, gene_symbol)
        }
    
        for (i in 1:length(voom.output[[cells]])) {
            comp = names(voom.output[[cells]])
            cells2 = paste0(cells, "_", comp)
            data1 = voom.output[[cells]][[i]]$LogFC
            names(data1) = rowData(se)$symbol
            data1 = sort(data1, decreasing = TRUE)
            data1_clean = data1[which(!is.na(names(data1)))]
            message(paste("analyzing", category[j],  names(voom.output[[cells]][i]), cells))
    
            msigdb <- GSEA(data1_clean, TERM2GENE = m_t2g)
            msigdb_results[[cells2]][[category[j]]][[i]] <- msigdb@result[order(msigdb@result$NES), ]
            msigdb_results[[cells2]][[category[j]]][[i]]  <-
                rbind(head(msigdb_results[[cells2]][[category[j]]][[i]] , 10), tail(msigdb_results[[cells2]][[category[j]]][[i]] , 10))
            msigdb_results[[cells2]][[category[j]]][[i]]$ID = factor(as.character(msigdb_results[[cells2]][[category[j]]][[i]]$ID),
                                                                  levels = unique(as.character(msigdb_results[[cells2]][[category[j]]][[i]]$ID)))
            msigdb_results[[cells2]][[category[j]]][[i]]$direction = ifelse(msigdb_results[[cells2]][[category[j]]][[i]]$NES > 0, "UP", "DOWN")
            msigdb_results[[cells2]][[category[j]]][[i]]$direction = factor(as.character(msigdb_results[[cells2]][[category[j]]][[i]] $direction),
                                                                         levels = c("UP", "DOWN"))
    
            #plot results
            msigplot[[paste0(cells2,".",category[j],".",names(voom.output[[cells2]][i]))]]=plotMsigplot(msigdb_results[[cells2]][[category[j]]][[i]] ,
                                                                                                    paste(cells2, category[j], names(voom.output[[cells2]][i])), 
                                                                                                    manual_color = c("red", "blue")) + 
                                                                                                    scale_x_discrete(label = function(x) stringr::str_trunc(x, 50))
        }
    }
}
heatmaplist.m=marrangeGrob(grobs = msigplot, nrow=2, ncol=2)
ggsave(paste0(outdir, "GSEA.pdf"), heatmaplist.m, width=12, height=9)
saveRDS(msigdb_results, paste0(outdir, "msigdb_results.rds"))
