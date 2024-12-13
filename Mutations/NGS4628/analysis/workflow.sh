# This code is solely intented provide a reference on how we obtained the mutations annotations referenced in the manuscript.
# Mutations/variants were called as described in the manuscript. 
# To run transvar, mutation were extracted from all the vcf files of each variant caller. In our case, we had mutation calls from MuTect2, Strelka and LoFreq. The input file has to be in the following format:

#>head ../OUTPUT/mutations_all_callers.txt
#chr10:100063634G>A
#chr10:100069832T>T
#chr10:100069841C>A
#chr10:100076084C>A
#chr10:100143841G>T
#...

mkdir -p ../OUTPUT/

# run transvar
module load oncbfx/transvar/2.5.10.20211024

transvar ganno -l ../OUTPUT/mutations_all_callers.txt --refversion hg38 --ensembl > ../OUTPUiT/mutations_all_callers_transvar.tsv
