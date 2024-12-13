library(ggplot2)
library(ggpubr)
library(grid)

# Merge figures
setwd("/gstore/project/tead/scMultiome/NGS4467_pilot_JT/analysis/")
curr_outdir = "../OUTPUT/archr/Plots/"

list.files("../OUTPUT/", pattern=".rds")

p_samples = readRDS("../OUTPUT/reference_ordinations.rds")
p_gex = readRDS("../OUTPUT/gex_ordinations_YAP_FOSL_TEAD.rds")
p_chromvar = readRDS("../OUTPUT/chromvar_ordinations_YAP_FOSL_TEAD.rds")
p_activity = readRDS("../OUTPUT/epiregulon_activity_ordinations.rds")
p_hippo = readRDS("../OUTPUT/hippo_signatures_gex_scores_ordinations.rds")

my_theme = theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
)

my_labels = p_samples[["my_labels"]]
my_labels2 = p_samples[["my_labels2"]]
my_labels3 = p_samples[["my_labels3"]]
p_samples[[1]]$layers[[1]]$aes_params$size = 0.2
p_samples[[2]]$layers[[1]]$aes_params$size = 0.2
p_samples[[3]]$layers[[1]]$aes_params$size = 0.2
figure <- ggarrange(
    p_samples[[1]] + ggtitle("samples") + theme(legend.text=element_text(size=6), legend.key.size=unit(0.05, "cm")) + my_theme, 
    p_samples[[2]] + ggtitle("clusters") + theme(legend.text=element_text(size=6), legend.key.size=unit(0.05, "cm")) + my_theme, 
    p_samples[[3]] + ggtitle("named_clusters") + theme(legend.text=element_text(size=6), legend.key.size=unit(0.05, "cm")) + my_theme, 
    labels = NULL,
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

figure2 <- ggarrange(
    p_gex[[1]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,  
    p_gex[[2]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,  
    p_gex[[3]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme, 
    p_gex[[4]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,
    labels = NULL,
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))
figure3 <- ggarrange(
    p_chromvar[[1]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme, 
    p_chromvar[[2]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme, 
    p_chromvar[[3]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme, 
    p_chromvar[[4]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,
    labels = NULL,
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

p_activity[[1]]$layers[[2]] = NULL #$aes_params$fill=NA
p_activity[[2]]$layers[[2]] = NULL
p_activity[[3]]$layers[[2]] = NULL
p_activity[[4]]$layers[[2]] = NULL
figure4 <- ggarrange(
    p_activity[[1]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,
    p_activity[[2]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,
    p_activity[[3]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,
    p_activity[[4]] + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")) + my_theme,
    labels = NULL,
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

figure5 <- ggarrange(
    p_hippo[[1]] + my_theme, p_hippo[[2]] + my_theme, p_hippo[[3]] + my_theme, p_hippo[[4]] + my_theme, 
    p_hippo[[5]] + my_theme, p_hippo[[6]] + my_theme, p_hippo[[7]] + my_theme, p_hippo[[8]] + my_theme, 
    p_hippo[[9]] + my_theme,
    labels = NULL,
    ncol=4, nrow=3, widths=c(0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

### Try to generate one row by row
figure5a <- ggarrange(
    p_hippo[[1]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")), 
    p_hippo[[2]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[3]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[4]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

figure5b <- ggarrange(
    p_hippo[[5]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")), 
    p_hippo[[6]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[7]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[8]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

figure5c <- ggarrange(
    p_hippo[[9]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[10]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[11]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    p_hippo[[12]] + my_theme + theme(legend.text=element_text(size=6, angle=90, vjust=0.5), legend.key.size=unit(0.4, "cm")),
    ncol=4, nrow=1, widths=c(0.5, 0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

final_figure = ggarrange(
    figure, figure2, figure3, figure4, figure5a, figure5b, figure5c,
    ncol=1, nrow=7, 
    heights=c(0.47, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    common.legend=F, legend="none",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))

ggsave(paste0(curr_outdir, "/Figure_umaps_for_paper_v4", ".pdf"), plot=final_figure, device="pdf", height=23, width=11.5, units="in", limitsize=FALSE)

#################
# Violin plots
list.files("../OUTPUT/", pattern=".rds")
p_violin_gex = readRDS("../OUTPUT/gex_violin_YAP_FOSL_TEAD.rds")
p_violin_chromvar = readRDS("../OUTPUT/chromvar_violin_YAP_FOSL_TEAD.rds")
p_violin_activity = readRDS("../OUTPUT/epiregulon_activity_violin.rds")
p_violin_activity = readRDS(file="/gstore/project/tead/scMultiome/NGS4467_pilot_JT/OUTPUT/archr//regulon/epiregulon_activity_violin.rds")

figure6 <- ggarrange(
    p_violin_gex[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title=element_text(size=9)) + ylab("Gene expression\n(log2 normalized counts)"),
    p_violin_gex[[2]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()), 
    p_violin_gex[[3]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()),
    p_violin_gex[[4]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()), 
    p_violin_chromvar[[1]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title=element_text(size=9)) + ylab("chromVAR (z-score)"),
    p_violin_chromvar[[2]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()), 
    p_violin_chromvar[[3]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()),
    p_violin_chromvar[[4]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank()), 
    p_violin_activity[[1]] + ylab("Activity") + theme(axis.title.x=element_blank(), axis.title=element_text(size=9)),
    p_violin_activity[[2]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
    p_violin_activity[[3]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
    p_violin_activity[[4]] + theme(axis.title.y=element_blank(), axis.title.x=element_blank()),
    labels = NULL,
    ncol=4, nrow=3, widths=c(0.5, 0.5, 0.5, 0.5),
    common.legend=T, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))
ggsave(paste0(curr_outdir, "/Figure_violin_plots_for_paper_v3", ".pdf"), plot=figure6, device="pdf", height=12, width=10, units="in", limitsize=FALSE)

p_violin_signatures = readRDS("../OUTPUT/hippo_signatures_gex_scores_violin.rds")
figure7 <- ggarrange(
    p_violin_signatures[[1]] + theme(axis.text.x=element_blank(), axis.title=element_text(size=9)) + ylab("Score"),
    p_violin_signatures[[2]] + theme(axis.title.y=element_blank(), axis.text.x=element_blank()), 
    p_violin_signatures[[3]] + theme(axis.title.y=element_blank(), axis.text.x=element_blank()),
    p_violin_signatures[[4]] + theme(axis.text.x=element_blank(), axis.title=element_text(size=9)) + ylab("Score"),
    p_violin_signatures[[5]] + theme(axis.title.y=element_blank(), axis.text.x=element_blank()),
    p_violin_signatures[[6]] + theme(axis.title.y=element_blank(), axis.text.x=element_blank()), 
    p_violin_signatures[[7]] + theme(axis.title=element_text(size=9)) + ylab("Score"),
    p_violin_signatures[[8]] + theme(axis.title.y=element_blank()),
    p_violin_signatures[[9]] + theme(axis.title.y=element_blank()),
    p_violin_signatures[[10]] + theme(axis.title=element_text(size=9)) + ylab("Score"),
    p_violin_signatures[[11]] + theme(axis.title.y=element_blank()),
    p_violin_signatures[[12]] + theme(axis.title.y=element_blank()), 
    labels = NULL,
    ncol=3, nrow=4, widths=c(0.5, 0.5, 0.5),
    common.legend=T, legend="bottom",
    align="hv", 
    font.label=list(size=10, color="black", face="bold", family=NULL, position="top"))
ggsave(paste0(curr_outdir, "/Figure_violin_plots_gene_signature_for_paper_v4", ".pdf"), plot=figure7, device="pdf", height=11.5, width=9, units="in", limitsize=FALSE)
