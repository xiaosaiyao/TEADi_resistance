#!/bin/bash
#! -n 4
#! -N 1
#! --mem=16G
#! --qos=medium

output_dir=/gstore/project/tead/Atacseq/H226_TEADi/OUTPUT/


cd $output_dir

export(G6915_unique_up, "OUTPUT/beds/G6915_up.unique.bed", "BED")
export(G7883_unique_up, "OUTPUT/beds/G7883_up.unique.bed", "BED")

for i in beds/G7883.G6915.down_overlap.bed beds/G6915_up.unique.bed G7883_up.unique.bed beds/G7883.G6915.up_overlap.bed # beds/G6915_down.unique.bed beds/G7883_down.unique.bed
do
findMotifsGenome.pl $output_dir/$i hg38 $output_dir/homer/`basename $i`.bg -size given -bg
done

