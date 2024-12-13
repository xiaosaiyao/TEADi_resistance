#!/usr/bin/env Rscript

# Function that merges seRNA objects with their corresponding archr object.
# Function should be called after 'demultiplex_hto_and_scRNA.R'
# In this function, we can remove HTOs and/or SAMs (usually after having run a first pass of the data).
# Author: Julien Tremblay - julien.tremblay@contractors.roche.com
# Roche/Genentech

merge_sc_and_archr <- function(hto_frs_id=NULL, htos_to_remove=NULL, sams_to_remove=NULL) {

    #hto_frs_id = "FRS17654"
    #htos_to_remove = c("HTO-9", "HTO-12")
    #sams_to_remove=c("SAM24416357")
    
    if(is.null(hto_frs_id)){
        stop("a hto_frs_id value has to be included")
    }
    
    library(ArchR)
    library(SingleCellExperiment)
    library(data.table)
    
    hto_info = maw.utils::getFireDBResourceSetInfo(frs.id = hto_frs_id)$file
    samples = unique(hto_info$sampleName)
    
    seRNA_objs = list()
    for(sample in samples){
        curr_seRNA_obj = readRDS(paste0("./output/", sample, "/rds/", sample, "_demultiplex_hto_and_scRNA.rds"))
        seRNA_objs[[sample]] = curr_seRNA_obj
    }
    
    # Here load the multiple seRNA objects from last step and combine them.
    seRNA_final = do.call(cbind, seRNA_objs)
    
    seRNA_final$hash_assignment2 = paste0(seRNA_final$library, "_", seRNA_final$hash_assignment)
    if(is.null(htos_to_remove)){
        message("No HTOs will be removed.")
    }else{
        found_htos_to_remove = seRNA_final$hash_assignment[seRNA_final$hash_assignment %in% htos_to_remove]
        if(length(found_htos_to_remove) == 0){
            message("Did not find any of the following HTOs to remove from the scRNA_final object... ")
            message(htos_to_remove)
        }else{
            message("Removing ", length(found_htos_to_remove), " cells associated with HTOs to remove from seRNA_final object.")
            idxs_to_keep = which(!seRNA_final$hash_assignment %in% htos_to_remove)
            seRNA_final = seRNA_final[, idxs_to_keep]
        }
    }
    
    if(is.null(sams_to_remove)){
        message("No SAMs will be removed.")
    }else{
        found_sams_to_remove = seRNA_final$library[seRNA_final$library %in% sams_to_remove]
        if(length(found_sams_to_remove) == 0){
            message("Did not find any of the following HTOs to remove from the scRNA_final object... ")
            message(sams_to_remove)
        }else{
            message("Removing ", length(found_sams_to_remove), " cells associated with SAMs to remove from seRNA_final object.")
            idxs_to_keep = which(!seRNA_final$library %in% sams_to_remove)
            seRNA_final = seRNA_final[, idxs_to_keep]
        }
    }
    
    proj = ArchR::loadArchRProject("data/ArchRProject/")
    common = intersect(proj$cellNames, colnames(seRNA_final))
    proj = proj[common,]
    
    # Here at this point, we have the archr object and seRNA obj witht the same cells.
    # add HTO information
    for(curr_cell_attribute in colnames(colData(seRNA_final))){
        proj = addCellColData(
            ArchRProj = proj,
            data = colData(seRNA_final)[common, curr_cell_attribute],
            cells = common,
            name = curr_cell_attribute,
            force = TRUE
        )
    }
    
    #filter out doublets and non-confident calls
    proj = proj[which(proj$Confident == TRUE & proj$Doublet == FALSE), ]
    
    # save
    saveArchRProject(proj, outputDirectory="./output/archr/", overwrite=TRUE)
    message("merge_sc_and_archr.R completed!")
}

usage=function(errM) {
    cat("\nUsage : Rscript my_script_name.R [option] <Value>\n")
    cat("       -b        : FRS[0-9]+ id of the hto object. Will solely be used to fetch sample IDs (SAM[0-9]+)\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 1) {
    usage("missing arguments")
    stop("missing -b <string> argument.")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-b") {
        hto_frs_id=ARG[i+1]
    }
}
merge_sc_and_archr(hto_frs_id)
