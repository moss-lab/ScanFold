#!/usr/local/bin/python3.6
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

ScanFold-Scan
Contact: Ryan Andrews - randrews@iastate.edu

This program takes a fasta input file and uses a scanning window approach to
calculate thermodynamic z-scores for individual windows.

Usage:
ScanFold-Scan_spinoff.py -i fasta.fa --scan_out_path ./out.tsv --zscore_wig_file_path ./zscore.bw --mfe_wig_file_path ./mfe.bw --ed_wig_file_path ./ed.bw --pvalue_wig_file_path ./pvalue.bw --fasta_file_path ./fasta.fasta.fa --fasta_index ./fasta.fai --nodeid "/scholar" --callbackurl "https://www.google.com"
"""

import time
import sys
import argparse
import string
import math
# not needed anymore
# import pyBigWig
import re
import numpy as np
sys.path.append('/home/randrews/ViennaRNA/lib/python3.6/site-packages/')
sys.path.append('/usr/local/lib/python3.6/site-packages')
#import RNA
import random

from ScanFoldSharedIGV import *

from Bio import SeqIO
# from progressbar import *               # just a simple progress bar


def write_fasta(sequence, outputfilename, name):
    w = open(outputfilename, 'w')
    fasta_sequence = sequence

    w.write(">"+name+"\n")
    w.write(str(fasta_sequence))

def randomizer(frag):
    result = ''.join(random.sample(frag,len(frag)))
    return result;

######Function to create X number of scrambled RNAs in list #################
#test
def scramble(text, randomizations, type):
    frag = str(text)
    frag_seqs = []
    if type == "di":
        for _ in range(randomizations):
            result = dinuclShuffle(frag)
            frag_seqs.append(result)
    elif type == "mono":
        frag_seqs = multiprocessing(randomizer, [frag for i in range(randomizations)], 4)

        # for _ in range(int(randomizations)):
        #     result = ''.join(random.sample(frag,len(frag)))
        #     frag_seqs.append(result)
    else:
        print("Shuffle type not properly designated; please input \"di\" or \"mono\"")

    return frag_seqs;

def write_wig(metric_list, step, name, outputfilename, strand):
    w = open(outputfilename, 'w')
    if strand == "forward":
        w.write("%s %s %s %s %s\n" % ("fixedStep", "chrom="+name, "start="+str(input_start_coordinate), "step="+str(step), "span="+str(step)))
        for metric in metric_list:
            #print(metric)
            try:
                if metric == "#DIV/0!":
                    w.write("%s\n" % (metric))
                else:
                    metric = float(metric)
                    w.write("%f\n" % (metric))
            except:
                print("Error writing WIG file.", metric)

    if strand == "reverse":
        seqend = input_start_coordinate + length
        num_windows = math.floor((length-window_size)/step)
        remainder = length-num_windows*step
        revstart = seqend - ((num_windows+2)*(step))
        w.write("%s %s %s %s %s\n" % ("fixedStep", "chrom="+name, "start="+str(revstart), "step="+str(step), "span="+str(step)))
        for metric in reversed(metric_list):
            #print(metric)
            try:
                if metric == "#DIV/0!":
                    w.write("%s\n" % (metric))
                else:
                    metric = float(metric)
                    w.write("%f\n" % (metric))
            except:
                print("Error writing WIG file.", metric)

def transcribe(seq):
    #Function to covert T nucleotides to U nucleotides
    for ch in seq:
        rna_seq = seq.replace('T', 'U')
        return(rna_seq)

if __name__ == "__main__":
    #### Parsing arguments ####
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input filename')
    parser.add_argument('-s', type=int, default=10,
                        help='step size')
    parser.add_argument('-w', type=int, default=120,
                        help='window size')
    parser.add_argument('-r', type=int, default=30,
                        help='randomizations')
    parser.add_argument('-t', '--temp', type=int, default=37,
                        help='Folding temperature')
    parser.add_argument('-type', type=str, default='mono',
                        help='randomization type')
    parser.add_argument('-p', type=str, default='off',
                        help='print to screen option (default off:1)')
    parser.add_argument('--print_random', type=str, default='off',
                        help='print to screen option (default off)')
    parser.add_argument('--name', type=str, default = "UserInput",
                        help='name of data being analyzed')
    parser.add_argument('--split', type=str, default = "off",
                        help='name of data being analyzed')
    parser.add_argument('-d', type=str, default = "forward",
                        help='strand of genome (forward or reverse; default forward)')

    ###Required arguments for webserver:
    parser.add_argument('--scan_out_path', type=str,
                        help='ScanFold-Scan output path')
    parser.add_argument('--zscore_wig_file_path', type=str,
                        help='zscore_wig_file_path')
    parser.add_argument('--mfe_wig_file_path', type=str,
                        help='mfe_wig_file_path')
    parser.add_argument('--ed_wig_file_path', type=str,
                        help='ed_wig_file_path')
    parser.add_argument('--pvalue_wig_file_path', type=str,
                        help='pvalue_wig_file_path')
    parser.add_argument('--fasta_file_path', type=str,
                        help='fasta_file path')
    parser.add_argument('--fasta_index', type=str,
                        help='fasta index file path')


    ### Arguments for IGV
    parser.add_argument('--start', type=int,
                        help='Input start coordinate')


    ### input parms ###



    args = parser.parse_args()

    myfasta = args.input
    step_size = int(args.s)
    window_size = int(args.w)
    randomizations = int(args.r)
    temperature = int(args.temp)
    type = str(args.type)
    print_to_screen = str(args.p)
    print_random = str(args.print_random)

    scan_out_path = args.scan_out_path

    mfe_wig_file_path = args.mfe_wig_file_path
    zscore_wig_file_path = args.zscore_wig_file_path
    pvalue_wig_file_path = args.pvalue_wig_file_path
    ed_wig_file_path = args.ed_wig_file_path

    name = args.name
    fasta_file_path = args.fasta_file_path
    fasta_index = args.fasta_index

    split = args.split

    strand = args.d

    ### IGV Args

    input_start_coordinate = args.start+1

    ##################### Main Script #########################################

    ### Initialize the bigwig files
    # MFE_wig = pyBigWig.open(mfe_wig_file_path, 'w')
    # zscore_wig = pyBigWig.open(zscore_wig_file_path, 'w')
    # pvalue_wig = pyBigWig.open(pvalue_wig_file_path, 'w')
    # ED_wig = pyBigWig.open(ed_wig_file_path, 'w')



    ### Create output file:
    w = open(scan_out_path, 'w')

    # create progress bar
    # widgets = ['Test: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),
    #            ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

    ### Create list for metrics to be written to bw via pyBigWig
    MFE_list = []
    zscore_list = []
    pscore_list = []
    ED_list = []


    ##### Establish empty lists to capture calculated metrics per window ######
    zscore_total = []
    numerical_z = []
    pscore_total = []
    numerical_p = []
    MFE_total = []
    ED_total = []


    seq = ""
    do_upper_transcribe = False
    do_write_header = False

    # with open(myfasta, 'r') as forward_fasta:
    with open (myfasta, 'r') as forward_fasta:
        if '>' in forward_fasta.readlines()[0]:
            #print("Sequence in FASTA format...")
            with open (myfasta, 'r') as forward_fasta:
                for cur_record in SeqIO.parse(forward_fasta, "fasta"):
                    ### get info about fasta files like name and sequence
                    #read_name = cur_record.name #this reads fasta header (not reliable)
                    seq = cur_record.seq
                    if "-" in seq:
                        seq.ungap("-")
                        print("Gaps found in sequence.\n Removing Gaps...")
                    do_upper_transcribe = True
                    do_write_header = True

        else:
            with open (myfasta, 'r') as forward_fasta:
                #print("Sequence not in FASTA format, attempting to convert. You can also resubmit sequence with a proper FASTA header (i.e. >SequenceName)")
                raw_sequence = forward_fasta.read()
                seq = re.sub("\n", "", raw_sequence.strip())
                #print(seq)

    seq = seq.upper()
    seq = transcribe(seq)
    length = len(seq)
    record_name = name

    print("Sequence Length: "+str(len(seq))+"nt")
    number_windows = ((length+1)-window_size)/step_size
    #print(number_windows)
    if number_windows < float(1.0):
        number_windows = ((length+1)-window_size)/step_size
    else:
        number_windows = math.floor(((length+1)-window_size)/step_size)

    #print(length, step_size, window_size)
    print("Approximately "+str(int(number_windows))+" windows will be generated.")
    print("Sequence being scanned...")
    if len(seq) > 40000 :
        #print(str(len(seq)))
        raise SystemExit('Input sequence is longer than 40000 nt; in order to scan longer sequences consider using the stand alone programs (avaiable here: https://github.com/moss-lab/ScanFold)')

    ### Need to calculate reverse strand start coordinate (flipping output)
    # print("Remove LINES about REVERSE STRAND")
    # strand = "reverse"
    if strand == "reverse":
        seqend = input_start_coordinate + (length-1)
        remainder = length-(number_windows*step_size)
        revstart = seqend - ((number_windows)*(step_size))
        #print(f"Reverse Strand. Nucleotide remainders={remainder}  Sequence end={seqend} Reverse start={revstart}")
    else:
        remainder = length-(number_windows*step_size)
        seqend = input_start_coordinate + (length-1)

    # ### Add headers for bigwig files
    # """ Headers need to have the name and length of fasta sequence
    # Should be set up like XYZ_wig.addHeader("chromosomeName", length)
    # """
    # MFE_wig.addHeader([(str(record_name), int(length))])
    # zscore_wig.addHeader([(str(record_name), int(length))])
    # pvalue_wig.addHeader([(str(record_name), int(length))])
    # ED_wig.addHeader([(str(record_name), int(length))])

    ### Write the header for output:
    if do_write_header:
        w.write("i\tj\tTemperature\tNative_dG\tZ-score\tP-score\tEnsembleDiversity\tSequence\tStructure\tCentroid\tWindowSize="+str(window_size)+" StepSize="+str(step_size)+" RandomizationCount="+str(randomizations)+" ShuffleType="+type+" Name="+record_name+"\n")
    i = 0
##### Main routine using defined functions: ##########################################
### Figure out length progress bar
    # pbar = ProgressBar(widgets=widgets, max_value=int(length))
    # print(int(length/step_size))
    # pbar.start()
    while i == 0 or i <= (length - window_size):
        #print(i)
        # pbar.update(i)
        start_nucleotide = i + input_start_coordinate # This will just define the start nucleotide coordinate value
        frag = seq[i:i+int(window_size)] # This breaks up sequence into fragments
        #print(frag)
        #print(str(len(frag)))
        start_nucleotide = i + input_start_coordinate
        end_nucleotide = start_nucleotide + window_size-1
        if -1 == 0:
            print("Magic")
        # if 'N' in frag:
        #     w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str(frag)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\n")
        #     i += step_size #this ensures that the next iteration increases by "step size" length
        else:
            #print(start_nucleotide)
            #print(end_nucleotide)
            if do_upper_transcribe:
                frag = frag.upper()
                frag = frag.transcribe()
            if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                MFE = int(0.0)
                zscore = "00.00"
                ED = int(0.0)
                pscore = int(0.0)
                structure = "........................................................................................................................"
                centroid = "........................................................................................................................"
            else:
                #print(frag)
                frag = frag.upper()
                frag = transcribe(frag)
                structure, centroid, MFE, ED = rna_fold(frag, temperature)

                MFE_total.append(float(MFE))
                ED_total.append(float(ED))
                seqlist = [] # creates the list we will be filling with sequence fragments
                seqlist.append(frag) # adds the native fragment to list
                scrambled_sequences = scramble(frag, randomizations, type)
                seqlist.extend(scrambled_sequences)
                energy_list = energies(seqlist, temperature)
                if print_random == "on":
                    print(energy_list)
                try:
                    zscore = round(zscore_function(energy_list, randomizations), 2)
                except:
                    zscore = zscore_function(energy_list, randomizations)
                zscore_total.append(zscore)

                #print(zscore)
                pscore = round(pscore_function(energy_list, randomizations), 2)
                #print(pscore)
                pscore_total.append(pscore)

### Append metrics to list ###
                MFE_list.append(MFE)
                zscore_list.append(zscore)
                pscore_list.append(pscore)
                ED_list.append(ED)

            if print_to_screen == 'on':
                print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
            w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
            #gff3file.write()
            #pscore_wig.write()
            #zscore_wig.write()
            #ED_wig.write()
            #MFE_wig.write()

            i += step_size #this ensures that the next iteration increases by "step size" length
    ### Add a final window if the full sequence was not covered ###
    #print("ELSE", remainder, seqend)
    #print(str(int(seqend+1)-window_size), str(seqend))
### scan a final window no matter what:
    if strand == "forward":
        start_nucleotide = str((seqend+1)-window_size)
        end_nucleotide = str(seqend)
        frag = seq[(int(length))-window_size:int(length)] # This breaks up sequence into fragments

    if strand == "reverse":
        start_nucleotide = str((seqend+1)-window_size)
        end_nucleotide = str(seqend)
        frag = seq[(int(length))-window_size:int(length)] # This breaks up sequence into fragments

    #print(frag)
    #print(str(len(frag)))
    if -1 == 0:
        print("Magic")
    # if 'N' in frag:
    #     w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str(frag)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\n")
    #     i += step_size #this ensures that the next iteration increases by "step size" length
    else:
        #print(start_nucleotide)
        #print(end_nucleotide)
        if do_upper_transcribe:
            frag = frag.upper()
            frag = frag.transcribe()
        if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
            MFE = int(0.0)
            zscore = "00.00"
            ED = int(0.0)
            pscore = int(0.0)
            structure = "........................................................................................................................"
            centroid = "........................................................................................................................"
        else:
            #print(frag)
            structure, centroid, MFE, ED = rna_fold(frag, temperature)

            MFE_total.append(float(MFE))
            ED_total.append(float(ED))
            seqlist = [] # creates the list we will be filling with sequence fragments
            seqlist.append(frag) # adds the native fragment to list
            scrambled_sequences = scramble(frag, randomizations, type)
            seqlist.extend(scrambled_sequences)
            energy_list = energies(seqlist, temperature)
            if print_random == "on":
                print(energy_list)
            try:
                zscore = round(zscore_function(energy_list, randomizations), 2)
            except:
                zscore = zscore_function(energy_list, randomizations)
            zscore_total.append(zscore)

            #print(zscore)
            pscore = round(pscore_function(energy_list, randomizations), 2)
            #print(pscore)
            pscore_total.append(pscore)

### Append metrics to list ###
            MFE_list.append(MFE)
            zscore_list.append(zscore)
            pscore_list.append(pscore)
            ED_list.append(ED)

        if print_to_screen == 'on':
            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
        w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
        #gff3file.write()
        #pscore_wig.write()
        #zscore_wig.write()
        #ED_wig.write()
        #MFE_wig.write()


    for z in zscore_total:
        try:
            numerical_z.append(float(z))
        except ValueError:
            continue

    for p in pscore_total:
        try:
            numerical_p.append(float(p))
        except ValueError:
            continue

    window_count = len(zscore_total)

    mean_pscore = round(np.mean(numerical_p), 2)
    mean_zscore = round(np.mean(numerical_z), 2)
    mean_MFE = round(np.mean(MFE_total), 2)
    mean_ED = round(np.mean(ED_total), 2)
    write_wig(MFE_list, step_size, name, mfe_wig_file_path, strand)
    write_wig(zscore_list, step_size, name, zscore_wig_file_path, strand)
    write_wig(pscore_list, step_size, name, pvalue_wig_file_path, strand)
    write_wig(ED_list, step_size, name, ed_wig_file_path, strand)

    write_fasta(seq, fasta_file_path, name)
    write_fai(seq, fasta_index, name)

    if split == "on":
        print("Completed ScanFold-Scan, intial results shown below.\n ScanFold-Fold is running now")

    print("Mean MFE = "+str(mean_MFE))
    print("Mean Z-score = "+str(mean_zscore))
    print("Mean P-value = "+str(mean_pscore))
    print("Mean Ensemble Diversity = "+str(mean_ED))
    print("ScanFold-Scan complete, find output files below")
