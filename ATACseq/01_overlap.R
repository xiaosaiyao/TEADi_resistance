library(ChIPpeakAnno)
library(rtracklayer)
library(GenomicRanges)

G6915_down <- "/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/G6915/down.bed"
G7883_down <- "/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/G7883/down.bed"
G6915_up <- "/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/G6915/up.bed"
G7883_up <- "/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/G7883/up.bed"

grl <- list()
for (bed in c("G6915_down", "G6915_up" ,"G7883_down", "G7883_up")){

    bedfiles <- read.delim(get(bed), header = FALSE)
    colnames(bedfiles) <- c("chr", "start","end","peak")
    grl[[bed]] <- makeGRangesFromDataFrame(bedfiles)
}

pdf("OUTPUT/overlap.6915.7883.pdf")
makeVennDiagram(list(grl$G6915_down, grl$G7883_down),
                                NameOfPeaks=c("G6915 down", "G7883 down"))

makeVennDiagram(list(grl$G6915_up, grl$G7883_up),
                           NameOfPeaks=c("G6915 up", "G7883 up"))
dev.off()

down_overlap <- subsetByOverlaps(grl$G6915_down, grl$G7883_down)
up_overlap <- subsetByOverlaps(grl$G6915_up, grl$G7883_up)
G6915_unique <- subsetByOverlaps(grl$G6915_down, grl$G7883_down, invert = TRUE)
G7883_unique <- subsetByOverlaps(grl$G7883_down, grl$G6915_down, invert = TRUE)
G6915_unique_up <- subsetByOverlaps(grl$G6915_up, grl$G7883_up, invert = TRUE)
G7883_unique_up <- subsetByOverlaps(grl$G7883_up, grl$G6915_up, invert = TRUE)

export(G6915_unique, "OUTPUT/beds/G6915_down.unique.bed", "BED")
export(G7883_unique, "OUTPUT/beds/G7883_down.unique.bed", "BED")
export(G6915_unique_up, "OUTPUT/beds/G6915_up.unique.bed", "BED")
export(G7883_unique_up, "OUTPUT/beds/G7883_up.unique.bed", "BED")
export(down_overlap, "OUTPUT/beds/G7883.G6915.down_overlap.bed", "BED")
export(up_overlap, "OUTPUT/beds/G7883.G6915.up_overlap.bed", "BED")
