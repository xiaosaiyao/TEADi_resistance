ml deeptools
ml BEDTools/2.30.0-GCCcore-6.3.0
ml spaces/oncbfx
ml oncbfx/tools/dev
ml oncbfx/wiggletools
ml singularity/3.5.3-foss-2017a

wiggletools.py -s ../OUTPUT/samplesheet_KOs.csv -d ../OUTPUT/tornado_plots/ -w 8 -f Factor --statistics mean
