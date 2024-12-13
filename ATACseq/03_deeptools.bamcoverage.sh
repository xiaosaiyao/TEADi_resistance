#!/bin/bash
#! -n 4
#! -N 1
#! --mem=60G
#! --qos=long

module add deeptools
dir=/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT
cd $dir
TEAD1_chip_bed1=/gstore/data/genomics/congee_rest_runs/63bf627a1673b2cbd51378df/SAM24425030/croo_output/peak/rep1/1_0E0S_01TSGenen_DMSO-1_TEAD1_hs_i02.R1.nodup_x_0_0D5Q_01TSGenen_Pooled_Input_hs_i78.R1.nodup.300K.regionPeak.gz
TEAD1_chip_bed2=/gstore/data/genomics/congee_rest_runs/63bf627a1673b2cbd51378df/SAM24425031/croo_output/peak/rep1/2_0E0T_01TSGenen_DMSO-2_TEAD1_hs_i04.R1.nodup_x_0_0D5Q_01TSGenen_Pooled_Input_hs_i78.R1.nodup.300K.regionPeak.gz
TEAD1_merged_bed=/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/merged.TEAD1.chipseq.bed

panTEAD_chip_bed1=/gstore/data/genomics/congee_rest_runs/63cb03581673b2cbd51602cb/SAM24425012/croo_output/peak/rep1/SAM24425012_H226_DMSO_48h_antiTEAD_CST_rep1_R1.nodup_x_SAM24425020_H226_DMSO_48h_Input_rep1_R1.nodup.300K.regionPeak.gz
panTEAD_chip_bed2=/gstore/data/genomics/congee_rest_runs/63cb06a61673b2cbd51603ad/SAM24425013/croo_output/peak/rep1/SAM24425013_H226_DMSO_48h_antiTEAD_CST_rep2_R1.nodup_x_SAM24425021_H226_DMSO_48h_Input_rep2_R1.nodup.300K.regionPeak.gz
panTEAD_merged_bed=/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/merged.panTEAD.chipseq.bed

cp $TEAD1_chip_bed1 beds/DMSO-1_TEAD1.bed.gz
cp $TEAD1_chip_bed2 beds/DMSO-2_TEAD1.bed.gz

cp $panTEAD_chip_bed1 beds/DMSO-1_panTEAD.bed.gz
cp $panTEAD_chip_bed2 beds/DMSO-2_panTEAD.bed.gz

gunzip beds/DMSO-1_TEAD1.bed.gz
gunzip beds/DMSO-2_TEAD1.bed.gz

gunzip beds/DMSO-1_panTEAD.bed.gz
gunzip beds/DMSO-2_panTEAD.bed.gz


#filter on qval (column 9)
for i in beds/DMSO-1_TEAD1.bed beds/DMSO-2_TEAD1.bed beds/DMSO-1_panTEAD.bed beds/DMSO-2_panTEAD.bed
do
awk '{ if ($9 >2 ) { print } }' $i > $i.filt.bed
done

cat  beds/DMSO-1_TEAD1.bed.filt.bed beds/DMSO-2_TEAD1.bed.filt.bed | sortBed | mergeBed > $TEAD1_merged_bed
cat  beds/DMSO-1_panTEAD.bed.filt.bed beds/DMSO-2_panTEAD.bed.filt.bed | sortBed | mergeBed > $panTEAD_merged_bed

bigwig_dir=/gstore/data/genomics/congee_rest_runs/61314afb17c8a42bb220e277/

DMSO=/DMSO/croo_output/signal/rep1/01_0A96_01E8Genen_DMSO-1_ATAC_hs_i201_r1.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig
G7883=/G7883/croo_output/signal/rep1/13_0A9H_01E8Genen_G7883-1_ATAC_hs_i212_r1.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig
G9886=/G9886/croo_output/signal/rep1/04_0A99_01E8Genen_G9886-1_ATAC_hs_i204_r1.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig
G9734=/G9734/croo_output/signal/rep1/10_09HE_01E8Genen_G9734-1_ATAC_hs_i204_r1.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig
G6915=/G6915/croo_output/signal/rep1/07_0A9C_01E8Genen_G6915-1_ATAC_hs_i207_r1.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig

computeMatrix reference-point --referencePoint center -S $bigwig_dir/*/croo_output/signal/rep*/*fc.signal.bigwig  -b 500 -a 500 -R $merged_bed  --skipZeros -o H226.TEADi.gz -p max/2
plotHeatmap -m H226.TEADi.gz -out H226.TEADi.pdf --colorMap plasma  --zMin -2 --zMax 2

computeMatrix reference-point --referencePoint center -S $bigwig_dir/$DMSO $bigwig_dir/$G9886 $bigwig_dir/$G6915  $bigwig_dir/$G9734 -b 500 -a 500 -R $merged_bed  --skipZeros -o H226.TEADi.TEAD1chip.gz -p max/2
plotHeatmap -m H226.TEADi.TEAD1chip.gz -out H226.TEADi.TEAD1chip.pdf --colorMap Greys Greens Blues Reds  --zMin 0 --zMax 40 --samplesLabel DMSO G9886 G6915 G9734
computeMatrix reference-point --referencePoint center -S $bigwig_dir/$DMSO $bigwig_dir/$G9886 $bigwig_dir/$G6915  $bigwig_dir/$G9734 -b 500 -a 500 -R $panTEAD_merged_bed  --skipZeros -o H226.TEADi.panTEADchip.gz -p max/2
plotHeatmap -m H226.TEADi.panTEADchip.gz -out H226.TEADi.panTEADchip.pdf --colorMap Greys Greens Blues Reds  --zMin 0 --zMax 40 --samplesLabel DMSO G9886 G6915 G9734 --kmeans 3 --outFileSortedRegions H226.k3.bed

up_7883=/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/G7883/up.bed
down_7883=/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/beds/G7883/down.bed
computeMatrix reference-point --referencePoint center -S $bigwig_dir/$DMSO $bigwig_dir/$G7883 $bigwig_dir/$G6915  -b 500 -a 500 -R $down_7883  --skipZeros -o H226.TEADi.7883.gz -p max/2
plotHeatmap -m H226.TEADi.7883.gz -out H226.TEADi.panTEADchip.pdf --colorMap Blues --samplesLabel DMSO G7883 G6915


