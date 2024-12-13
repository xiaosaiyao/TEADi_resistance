#!/usr/bin/env python
 
"""A wrapper script to run wiggle tools for all reps belonging to the same samples/libs.
tested with modules:
ml deeptools
ml BEDTools/2.30.0-GCCcore-6.3.0
ml spaces/oncbfx
ml oncbfx/wiggletools
ml singularity/3.5.3-foss-2017a

Julien Tremblay - julien.tremblay@contractors.roche.com
"""
 
import os
import sys
import argparse
import re
import multiprocessing
import errno
import csv
import pprint
from collections import defaultdict 

def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--samplesheet', required=False, help='samplesheet in format required by dba. Either -i and -o args are given OR the -s/--samplesheet argument.', type=argparse.FileType('r'))
    parser.add_argument('-d', '--outdir', help="Output directory. To be use in combination with -s/--samplesheet arg.")
    parser.add_argument('-p', '--peaks-identifier', help="column name corresponding to the peaks in the sample sheet.", default="Peaks")
    parser.add_argument('-b', '--bigwigs-identifier', help="column name corresponding to the bigwigs in the sample sheet.", default="bigwigs")
    parser.add_argument('-f', '--factors-identifier', help="column name corresponding to the factor in the sample sheet.", default="Factor")
    parser.add_argument('-x', '--statistics', help="Choose between mean and median.", default="mean", choices=["mean", "median"])
    parser.add_argument('-t', '--num-threads', help="Number of threads", type=int, default=1)
    parser.add_argument('-w', '--width', help="May need to be adjusted depending on samples labels lengths", type=int, default=4)
    parser.add_argument('--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('--debug', default=False, action='store_true', help='Debug mode, print commands only.')
    
    exe = which("computeMatrix")
    wiggletools_sif = which("wiggletools.sif")
    if(exe == None):
        print("computeMatrix executable is not found in your environment.", file=sys.stderr)
        exit(1)
    if(wiggletools_sif == None):
        print("computeMatrix executable is not found in your environment.", file=sys.stderr)
        exit(1)

    args = parser.parse_args(arguments)
    num_threads = args.num_threads
    outdir = args.outdir
    factor_column_id = args.factors_identifier

    infiles_list = []
    outdirs_list = []
    out = defaultdict(dict)
    if args.samplesheet:
        samplesheet = os.path.abspath(args.samplesheet.name)
        readset_csv = csv.DictReader(open(samplesheet, 'r'), delimiter='\t')
        for line in readset_csv:
            if 'Replicate' in line:
                rep = line['Replicate']
            if 'Rep' in line:
                rep = line['Rep']
            name = line['Treatment'] + "_" + rep
            peaks_file = line[args.peaks_identifier]
            bigwig_file = line[args.bigwigs_identifier]
            factor = line[factor_column_id]
            color = line['Color']
            treatment = line['Treatment']

            if treatment not in out[factor]:
                out[factor][treatment] = defaultdict(dict)

            if "colors" not in out[factor][treatment]:
                out[factor][treatment]["colors"] = []
            if "peaks" not in out[factor][treatment]:
                out[factor][treatment]["peaks"] = []
            if "bigwigs" not in out[factor][treatment]:
                out[factor][treatment]["bigwigs"] = []
            if "names" not in out[factor][treatment]:
                out[factor][treatment]["names"] = []
    
            out[factor][treatment]["peaks"].append(peaks_file)
            out[factor][treatment]["bigwigs"].append(bigwig_file)
            out[factor][treatment]["names"].append(name)
            out[factor][treatment]["colors"].append(color)
            
            if(args.verbose):
                print(name + "\t" + rep + "\t" + peaks_file + "\t" + bigwig_file, file=sys.stderr)
        
        pp = pprint.PrettyPrinter(depth=4)
        pp.pprint(out)

        for factor in out:
            print(factor)
            colors = ""

            wiggletools_names = []
            wiggletools_bigwigs = []
            peaks_to_merge = []
            my_colors = []
            for treatment in out[factor]:
                print(treatment)
                my_colors.append(out[factor][treatment]["colors"][1])

                print("#1 ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", file=sys.stderr)
                wiggletools_names.append(factor + "_" + treatment)
                wiggletools_bigwigs.append(os.path.join(outdir, factor + "_" + treatment + "_" + args.statistics + ".bigwig"))
                for my_peak_file in out[factor][treatment]["peaks"]:
                    peaks_to_merge.append(my_peak_file)

                cmd = ""
                try:
                    # wiggletools workflow
                    cmd = "singularity exec --bind /local --bind /gstore/data/genomics/pipeline_resources/cutnrun/HiChIP/HiC-Pro/hg38/GRCh38_noalt_as/ --bind /gstore/data/genomics/congee_rest_runs/ --bind $(pwd)/../ -W $(pwd) " + wiggletools_sif + " \\\n"
                    cmd += "wiggletools write_bg " + os.path.join(outdir, factor + "_" + treatment + "_" + args.statistics + ".bedGraph") + " " + args.statistics  + " " + " ".join(out[factor][treatment]["bigwigs"])
                    cmd += " && \\\nbedGraphToBigWig " + os.path.join(outdir, factor + "_" + treatment + "_" + args.statistics + ".bedGraph") + " /gstore/data/genomics/pipeline_resources/cutnrun/HiChIP/HiC-Pro/hg38/GRCh38_noalt_as/GRCh38_noalt_as.chrom.sizes " + os.path.join(outdir, factor + "_" + treatment + "_" + args.statistics + ".bigwig")
                    print(cmd)
                    
                    if not args.debug:
                        if os.system(cmd) != 0:
                            raise Exception(cmd + ' command failed...')
                except:
                    print(cmd + ' command failed...', file=sys.stderr)

            print("#2 ==> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::", file=sys.stderr)
            cmd2 = ""
            try:
            # merge peaks to be plotted on the same figure.
                cmd2 = "cat " + " ".join(peaks_to_merge) + "| sort -k1,1 -k2,2n | bedtools merge -i - > " + os.path.join(outdir, "merged_" + factor + ".bed")
                cmd2 += " && \\\ncomputeMatrix reference-point --referencePoint center -R " + os.path.join(outdir, "merged_" + factor + ".bed") + " -S " + " ".join(wiggletools_bigwigs) + " -b 500 -a 500 --skipZeros -p max/2 -o " + os.path.join(outdir, factor + "_" + args.statistics + ".gz")
                cmd2 += " && \\\nplotHeatmap -m " + os.path.join(outdir, factor + "_" + args.statistics + ".gz") + " -out " + os.path.join(outdir, factor + "_" + args.statistics + ".png") + " --samplesLabel " +  " ".join(wiggletools_names) + " --regionsLabel 'peaks' --colorMap " + " ".join(my_colors) + " --missingDataColor 1  --heatmapWidth " + str(args.width)
                cmd2 += " && \\\nplotHeatmap -m " + os.path.join(outdir, factor + "_" + args.statistics + ".gz") + " -out " + os.path.join(outdir, factor + "_" + args.statistics + ".pdf") + " --samplesLabel " +  " ".join(wiggletools_names) + " --regionsLabel 'peaks' --colorMap " + " ".join(my_colors) + " --missingDataColor 1  --heatmapWidth " + str(args.width)
                print(cmd2, file=sys.stderr)
            
                if not args.debug:
                    if os.system(cmd2) != 0:
                        raise Exception(cmd2 + ' command failed...')
            except:
                print(cmd2 + ' command failed...', file=sys.stderr)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

