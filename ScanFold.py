#!/usr/bin/python3
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

ScanFold
Contact: Ryan Andrews - randrews@iastate.edu
This program takes a fasta input file and uses a scanning window approach to
calculate thermodynamic z-scores for individual windows.

Usage:
$ python3.6 ScanFold.py fasta_filename [options]

"""

# ###For drawing images
# import matplotlib.pyplot as plt
# import forgi.visual.mplotlib as fvm
# import forgi

from ScanFoldFunctions import *

import os
import sys
import argparse
import string
import re
import math
import statistics

### MUST APPEND PATH FOR "import RNA" PYTHON BINDINGS !!!NOT FUNCTIONAL ON WINDOWS!!! ###
sys.path.append('/work/LAS/wmoss-lab/ViennaRNA/lib/python3.6/site-packages/')
sys.path.append('/home/randrews/ViennaRNA/lib/python3.6/site-packages')
import RNA

# sys.path.append('/usr/local/RNAstructure')
#import RNAstructure

import random
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from Bio import SeqIO
import time
from datetime import datetime
from datetime import timedelta
import uuid

start_time = time.time()

#### Parsing arguments ####
parser = argparse.ArgumentParser()
parser.add_argument('filename',  type=str,
                    help='input filename')

### New soft-constraint reactivity arguments ###
parser.add_argument('--react', type=str,
                    help='input SHAPE reactivity file')
parser.add_argument('-m', type=float, default=0.8,
                    help='SHAPE slope value')
parser.add_argument('-b', type=float, default=-0.2,
                    help='SHAPE intercept value')
parser.add_argument('--shapeD', action='store_true',
                    help='Use Deigan SHAPE algorithm')
parser.add_argument('--shapeZ', action='store_true',
                    help='Use Zarringhalam SHAPE algorithm')

### Core arguments
parser.add_argument('--name', type=str, default = "UserInput",
                    help='Name or ID of sequence being analyzed. Default "UserInput"')
parser.add_argument('--fold_only', action='store_true',
                    help='Utilizing output from a ScanFold-Scan, only run the ScanFold-Fold algorithm')
parser.add_argument('--dont_scan', action='store_true',
                    help='dont scan oprtion')
parser.add_argument('--fold', action='store_true', default=True,
                    help='Run ScanFold-Scan and Fold; on by default')
parser.add_argument('-f', type=int, default=-2,
                    help='z-score filter value; legacy value. Default = -2')
parser.add_argument('-c', type=int, default=1,
                    help='Competition on or off (1 or 0). Default = 1 (on)')
parser.add_argument('--out_name',  type=str,
                    help='Name of output folder (defaults to header name or date/time)', default=None)
parser.add_argument('-s', type=int, default=1,
                    help='Step size; default = 1')
parser.add_argument('-w', type=int, default=120,
                    help='Window size; default = 120')
parser.add_argument('-r', type=int, default=100,
                    help='Number of randomizations for background shuffling; default = 100')
parser.add_argument('-t', type=int, default=37,
                    help='Folding temperature in celsius; default = 37C')
parser.add_argument('--type', type=str, default='mono',
                    help='Randomization type')
parser.add_argument('--print', action='store_true',
                    help='Print to screen option (default off)')
parser.add_argument('--print_random', action='store_true',
                    help='Print all random MFEs to screen option (default off)')
parser.add_argument('--algo', type=str, default='rnafold',
                    help='Select RNA folding algorithm')
parser.add_argument('--constraints', type=str,
                    help='optional | input constraint file (as a dot-bracket format see https://www.tbi.univie.ac.at/RNA/RNAfold.1.html#heading6)')
parser.add_argument('--span', type=int,
                    help='optional | set a max bp span (default off)')
parser.add_argument('--global_refold', action='store_true',
                    help='Global refold option. Refold full sequence using Zavg <-1 and <-2 base pairs')


### New feature testing  ###
parser.add_argument('--lri', action='store_true',
                    help='optional | check for long range interactions')
parser.add_argument('--kmer', type=int, default=20,
                    help='size of k-mer for lri scan')
parser.add_argument('--kmer_step_size', type=int, default=1,
                    help='steps size of k-mer scan for lri ')
parser.add_argument('--lri_cutoff', type=int, default=-25,
                    help='steps size of k-mer scan for lri ')
parser.add_argument('--by_ed', action='store_true',
                    help='Select best bp based on lowest normalized ED value')


### Developer options ###
parser.add_argument('--out1', type=str, default = "./ScanFold.NoFilter",
                    help='out1 path (no filter)')
parser.add_argument('--out2', type=str, default = "./ScanFold.-1Filter",
                    help='out3 path (-1 filter)')
parser.add_argument('--out3', type=str, default = "./ScanFold.-2Filter",
                    help='out3 path (-2 filter)')
parser.add_argument('--out4', type=str, default = "./ScanFold.Log.txt",
                    help='log_file path (log file)')
parser.add_argument('--out5', type=str, default = "./ScanFold.FinalPartners.txt",
                    help='final_partner_file path')
parser.add_argument('--out6', type=str, default = "./IGV_BP_Track",
                    help='bp_track_file path')
# parser.add_argument('--out7', type=str, default = "./UserInput",
#                     help='fasta_file path')
parser.add_argument('--fasta_index', type=str, default = "./user_input.fai",
                    help='fasta index file path')
parser.add_argument('--dbn_file_path', type=str, default = "AllDBN-global_refold.txt",
                    help='dbn_file_path (all dbn refold)')
parser.add_argument('--dbn_file_path1', type=str, default = "Zavg_NoFilter",
                    help='dbn_file_path1 Zavg NoFilter')
parser.add_argument('--dbn_file_path2', type=str, default = "Zavg_-1_pairs",
                    help='dbn_file_path2, Zavg <-1')
parser.add_argument('--dbn_file_path3', type=str, default = "Zavg_-2_pairs",
                    help='dbn_file_path3, Zavg <-2')
parser.add_argument('--dbn_file_path4', type=str, default = "AllDBN.txt",
                    help='dbn_file_path4 (AllDBNs)')
parser.add_argument('--structure_extract_file', type=str, default = "ExtractedStructures.gff3",
                    help='structure_extract_file path')
parser.add_argument('--final_partners_wig', type=str, default = "./IGV_BP_Zavg_metrics",
                    help='final partners wig file path')

#scan args
args = parser.parse_args()
myfasta = args.filename
filename = args.filename
step_size = int(args.s)
window_size = int(args.w)
randomizations = int(args.r)
temperature = int(args.t)
type = str(args.type)
print_to_screen = args.print
print_random = args.print_random
algo = str(args.algo)
constraints = args.constraints
react = args.react
max_span = args.span
lri = args.lri
kmer = args.kmer
kmer_step_size = args.kmer_step_size
lri_cutoff = args.lri_cutoff
dont_scan = args.dont_scan
fold = args.fold
fold_only = args.fold_only
by_ed = args.by_ed
out_name = args.out_name
slope = args.m
intercept = args.b
shapeD = args.shapeD
shapeZ = args.shapeZ

#fold args
filter = int(args.f)
competition = int(args.c)
out1 = args.out1
out2 = args.out2
out3 = args.out3
out4 = args.out4
out5 = args.out5
out6 = args.out6
#out7 = args.out7
name = args.name
dbn_file_path = args.dbn_file_path
dbn_file_path1 = args.dbn_file_path1
dbn_file_path2 = args.dbn_file_path2
dbn_file_path3 = args.dbn_file_path3
dbn_file_path4 = args.dbn_file_path4

global_refold = args.global_refold

structure_extract_file = args.structure_extract_file

final_partners_wig = args.final_partners_wig

fasta_index_path = args.fasta_index
# nodeid = args.nodeid
# callbackurl = args.callbackurl


### Can hard code model options here ###
md = RNA.md()
md.temperature = int(temperature)
if max_span != None:
    md.max_bp_span = int(max_span)

### Define any funcitons that wont work if put into ScanFoldFunctions module
def getShapeDataFromFile(filepath):
    """Function adapted from Ronny Lorenz github page (Institute for Theoretical
    Chemistry University of Vienna)https://github.com/ViennaRNA/ViennaRNA/blob/master/tests/python/test-RNA-pf_window.py """
    retVec = []
    retVec.append(-999.0);  # data list is 1-based, so we add smth. at pos 0
    count=1
    with open(str(original_directory)+"/"+filepath, 'r') as f:
        lines = f.read().splitlines()
        #print(len(lines[0]))
        if len(lines[0].split('\t')) == 3:
            for line in lines:
                pos = int(line.split('\t')[0])
                nt = str(line.split('\t')[1])
                value = (line.split('\t')[2])
                if value == "NA":
                    value = -999

                if(pos==count):
                    retVec.append(float(value))
                else:
                    for i in range(pos-count):
                        retVec.append(-999.0)
                    retVec.append(float(value))
                    count=pos
                count+=1

        elif len(lines[0].split('\t')) == 2:
            print("no nucleotides detected in reactivity file")
            for line in lines:
                pos = int(line.split('\t')[0])
                value = (line.split('\t')[1])
                if value == "NA":
                    value = -999

                if(pos==count):
                    retVec.append(float(value))
                else:
                    for i in range(pos-count):
                        retVec.append(-999.0)
                    retVec.append(float(value))
                    count=pos
                count+=1
        else:
            raise("Trouble parsiging reactivity data")
    return retVec


##################### Main Script #########################################
original_directory = str(os.getcwd())
print(original_directory)
step_size = step_size
window_size = window_size
randomizations = randomizations
with open(myfasta, 'r') as forward_fasta:
    for cur_record in SeqIO.parse(forward_fasta, "fasta"):
        #must reinitialize dictionaries/values/directory for every fasta record
        #### Defining  variables ###############
        final_partners = {}
        best_bps ={}
        bp_dict = {}
        nuc_dict= {}
        nuc_dict_refold = {}
        cur_record_length = len(cur_record.seq)
        seq = cur_record.seq
        seq = seq.transcribe()
        if "-" in seq:
            raise("Gaps found in sequence. Please submit a complete sequence to ScanFold ")
        #print(str(seq))
        nuc_dict = {}
        nuc_i = 1
        read_name = cur_record.name
        cwd = os.getcwd()
        full_header = read_name
    ##### Establish empty lists to capture calculated metrics per window ######
        zscore_total = []
        numerical_z = []
        pvalue_total = []
        numerical_p = []
        MFE_total = []
        ED_total = []
        constraint_dict = {}
        constraint_list = []

        ###Try to split header
        if "|" in str(read_name):
            #print("Split")
            split_read_name = re.split('\|', str(read_name))
            read_name = split_read_name[0]
            full_header = str(read_name)
            print(read_name)

        ## this code makes an output directory for each fasta entry
        try:
            #Try to make a folder based on fasta name
            if out_name != None:
                try:
                    print("Making output folder named:"+out_name)
                    os.mkdir(cwd+"/"+out_name)
                    os.chdir(cwd+"/"+out_name)
                    folder_name = str(out_name)
                except:
                    now = datetime.now()
                    date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                    os.mkdir(cwd+"/"+out_name+"_"+date_time)
                    os.chdir(cwd+"/"+out_name+"_"+date_time)
                    folder_name = str(out_name+"_"+date_time)
                    print("Making output folder named:"+out_name+"_"+date_time)

            else:
                print("Making output folder named:"+read_name)
                os.mkdir(cwd+"/"+read_name)
                os.chdir(cwd+"/"+read_name)
                folder_name = str(read_name)
        except:
            #if folder already exists, try some other names
            try:
                #Try to make a folder based on "FULL" fasta name
                print("Folder already exists, attempting to extend name if more detailed header is available:"+full_header)
                if out_name != None:
                    now = datetime.now() # current date and time
                    #unique_filename = str(uuid.uuid4())
                    # date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                    date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                    date_time = date_time + "_" + str(random_with_N_digits(3))
                    folder_name = str(out_name+"_"+date_time)
                    print("Making output folder named:"+folder_name)
                    os.mkdir(cwd+"/"+folder_name)
                    os.chdir(cwd+"/"+folder_name)
                    folder_name = str(folder_name)
                else:
                    now = datetime.now() # current date and time
                    #unique_filename = str(uuid.uuid4())
                    # date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                    date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                    #date_time = date_time + "_" + str(random_with_N_digits(3))
                    folder_name = str("SF_Results_"+date_time)
                    print("Making output folder named:"+folder_name)
                    os.mkdir(cwd+"/"+folder_name)
                    os.chdir(cwd+"/"+folder_name)
                    folder_name = str(folder_name)
        #print(step_size, window_size, randomizations)
                #Make a folder based on date and time
            except:
                now = datetime.now() # current date and time
                #unique_filename = str(uuid.uuid4())
                # date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                date_time = now.strftime("%m-%d-%Y-%H.%M.%S")
                date_time = date_time + "_" + str(random_with_N_digits(3))
                folder_name = str("ScanFold-Results_"+date_time)
                print("Folder already exists, creating output name based on time:"+folder_name)
                os.mkdir(cwd+"/"+folder_name)
                os.chdir(cwd+"/"+folder_name)
        output = str(read_name+"."+str(filename)+".ScanFold.")
        outname = read_name+".win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)
        print("Output name="+str(outname))
        if len(cur_record) < window_size:
            print(cur_record.name+" sequence is less than window size. Moving on to next entry.")
            os.chdir(str(original_directory))
            del final_partners
            continue

        ### Added functionality to generate ScanFold models based on ED ###
        if by_ed == False:
            log_total = open(outname+".ScanFold.log", 'w')
            log_win = open(outname+".ScanFold.FinalPartners.txt", 'w')

        if by_ed == True:
            log_total = open(outname+".ScanFold.ED-weighted.log", 'w')
            log_win = open(outname+".ScanFold.ED-weighted.FinalPartners.txt", 'w')
        sirna_log = open(outname+".ntPairCounts.log", 'w')
        sirna_log.write(f"i\tnuc\twindows\tbps\n")

        w = open(outname+".out", 'w')
        if lri == True:
            print("WARNING! Using experimental LRI fuction. This has not been extensively tested, you may experience errors.")
            lri_file = open(outname+".LRI.out", 'w')

        for nuc in seq:
            #print("nuc dict", nuc, nuc_i, nuc_i)
            x = NucZscore(nuc,(nuc_i))
            nuc_dict[x.coordinate] = x
            nuc_i += 1

        if args.constraints != None:
            print("Considering constraint input")
            constraint_file = open(args.constraints, "r")
            constraints = constraint_file.readlines()[2]

            i = 1
            for nt in constraints:
                constraint_list.append(nt)
                constraint_dict[i] = nt
                i += 1

        if args.react != None:
            print("Considering SHAPE reactivity input")
            reactivities = getShapeDataFromFile(react)

        w.write("i\tj\tTemperature\tNative_dG\tZ-score\tP-score\tEnsembleDiversity\tSequence\tStructure\tCentroid\t"+read_name+"\n")


    ##### Main routine using defined functions: ##########################################
        i = 0
        if dont_scan == False and lri == False:
            ### Create list for metrics to be written to bw via pyBigWig
            MFE_list = []
            zscore_list = []
            pvalue_list = []
            ED_list = []

            print("Scanning input sequence:", read_name)
            while i == 0 or i <= (cur_record_length - window_size):
                #print("test")
                start_nucleotide = i + 1 # This will just define the start nucleotide coordinate value
                frag = seq[i:i+int(window_size)] # This breaks up sequence into fragments
                gc_content = get_gc_content(frag)
                #print(frag)
                #print(str(len(frag)))
                start_nucleotide = i + 1
                end_nucleotide = i + window_size
                #print(start_nucleotide)
                if algo == "rnastructure":
                    #print(frag)
                    #print("algo checked")
                    # print(start_nucleotide)
                    # print(end_nucleotide)
                    #frag = frag.transcribe()
                    #print(frag)
                    if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                        MFE = int(0.0)
                        zscore = "#DIV/0"
                        ED = int(0.0)
                        pvalue = int(0.0)
                        structure = "........................................................................................................................"
                        centroid = "........................................................................................................................"
                    else:
                        p = RNAstructure.RNA.fromString(str(frag))
                        # print(p)
                        fold = p.FoldSingleStrand(mfeonly=True)
                        # print(fold)
                        MFE = p.GetFreeEnergy(1)
                        #dot = p.WriteDotBracket(1)
                        structure = ''.join(get_structure(p))
                        #print(mfe)
                        pf = p.PartitionFunction()
                        #print(str(pf))
                        seqlist = [] # creates the list we will be filling with sequence fragments
                        seqlist.append(frag) # adds the native fragment to list
                        scrambled_sequences = scramble(frag, randomizations, type)
                        seqlist.extend(scrambled_sequences)
                        energy_list = energies(seqlist, temperature, algo)
                        try:
                            zscore = round(zscore_function(energy_list, randomizations), 2)
                        except:
                            zscore = zscore_function(energy_list, randomizations)
                        #print(zscore)
                        pvalue = round(pvalue_function(energy_list, randomizations), 2)
                        centroid = "........................................................................................................................"
                        ED = int(0.0)
                    if print_to_screen == True:
                        print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pvalue)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
                    w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pvalue)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
                    i += step_size #this ensures that the next iteration increases by "step size" length

                else:
                    #print(start_nucleotide)
                    #print(end_nucleotide)
                    #frag = frag.transcribe()
                    if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                        MFE = int(0.0)
                        zscore = "#DIV/0"
                        ED = int(0.0)
                        pvalue = int(0.0)
                        structure = "........................................................................................................................"
                        centroid = "........................................................................................................................"
                    else:
                        fc = RNA.fold_compound(str(frag), md) #creates "Fold Compound" object

                        if (constraints or react) == None:
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment

                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                            #fmfe = fc.pbacktrack()
                            #print(str(fmfe))
                        elif constraints != None:
                            window_constraint_list = constraint_list[start_nucleotide-1:end_nucleotide]
                            window_constraints = ''.join(window_constraint_list)
                            #print(len(window_constraints))
                            fc.hc_add_from_db(window_constraints)
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = fc.pf_fold() # calculate partition function "fold" of fragment
                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)

                        elif react != None:
                            window_react_list = reactivities[start_nucleotide:end_nucleotide+1]
                            print(window_react_list)
                            fc.pf()# performs partition function calculations
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                            ###
                            ###!!!Calculate shape MFE AFTER partition function!!!###
                            ###
                            #print(len(window_react_list))
                            if shapeD == True:
                                fc.sc_add_SHAPE_deigan(window_react_list, slope, intercept)
                            elif shapeZ == True:
                                fc.sc_add_SHAPE_zarringhalam(window_react_list)
                            else:
                                print("Using default Deigan algorthm.")
                                fc.sc_add_SHAPE_deigan(window_react_list, slope, intercept)

                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            #frag_q = fc.pf_fold() # calculate partition function "fold" of fragment
                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)

                        seqlist = [] # creates the list we will be filling with sequence fragments
                        seqlist.append(frag) # adds the native fragment to list
                        scrambled_sequences = scramble(frag, randomizations, type)
                        seqlist.extend(scrambled_sequences)
                        energy_list = energies(seqlist, temperature, algo)
                        if print_random == True:
                            print(energy_list)
                        try:
                            zscore = round(zscore_function(energy_list, randomizations), 2)
                        except:
                            zscore = zscore_function(energy_list, randomizations)
                        zscore_total.append(zscore)

                        #print(zscore)
                        pvalue = round(pvalue_function(energy_list, randomizations), 2)
                        #print(pvalue)
                        pvalue_total.append(pvalue)
                        #Iterate through dot bracket structure to determine locations of bps
                        if fold  == True:
                            #Convert sequence and structures into lists
                            sequence = list(frag)
                            structure = list(structure)

                            #Define window coordinates as string
                            #window = str(str(icoordinate)+"-"+str(jcoordinate))

                            #Determine length of window
                            length = len(sequence)

                            fold_i = 0
                            while fold_i < length:
                                #Unpaired nucleotide
                                #print(length, fold_i)
                                #print(len(structure))
                                if structure[fold_i] == '.':
                                    nucleotide = sequence[fold_i]
                                    coordinate = (fold_i + int(start_nucleotide))
                                    x = NucPair(nucleotide, coordinate, nucleotide, coordinate,
                                                zscore, MFE, ED)
                                    try:
                                        y = bp_dict[coordinate]
                                        y.append(x)
                                    except:
                                        bp_dict[coordinate] = []
                                        y = bp_dict[coordinate]
                                        y.append(x)
                                    fold_i += 1
                                #Paired nucleotide
                                else:
                                    fold_i += 1

                            #Inititate base pair tabulation variables
                            bond_order = []
                            bond_count = 0

                            #Iterate through sequence to assign nucleotides to structure type
                            m = 0
                            while  m < length:

                                if structure[m] == '(':
                                    bond_count += 1
                                    bond_order.append(bond_count)
                                    m += 1
                                elif structure[m] == ')':
                                    bond_order.append(bond_count)
                                    bond_count -= 1
                                    m += 1
                                elif structure[m] == '.':
                                    bond_order.append(0)
                                    m += 1
                                else:
                                    print("Error1")

                            #Initiate base_pair list
                            base_pairs = []

                            #Create empty variable named test
                            test = ""

                            #Iterate through bond order
                            j = 0
                            while j < len(bond_order):
                                if bond_order[j] != 0:
                                    test = bond_order[j]
                                    base_pairs.append(j+1)
                                    bond_order[j] = 0
                                    j += 1
                                    k = 0
                                    while k < len(bond_order):
                                        if bond_order[k] == test:
                                            base_pairs.append(k+1)
                                            bond_order[k] = 0
                                            k += 1
                                        else:
                                            k += 1
                                else:
                                    j += 1

                            #Iterate through "base_pairs" "to define bps
                            l = 0
                            while l < len(base_pairs):
                                lbp = base_pairs[l]
                                rbp = base_pairs[l+1]

                                lb = str(sequence[int(lbp)-1])
                                rb = str(sequence[int(rbp)-1])

                                lbp_coord = int(int(lbp)+int(start_nucleotide)-1)
                                rbp_coord = int(int(rbp)+int(start_nucleotide)-1)
                                x = NucPair(lb, lbp_coord, rb, rbp_coord, zscore, MFE, ED)
                                z = NucPair(rb, rbp_coord, lb, lbp_coord, zscore, MFE, ED)

                                #Try to append i-j pair to i-nuc for left i-nuc
                                try:
                                    y = bp_dict[lbp_coord]
                                    y.append(x)
                                #If i-nuc not defined, define it
                                except:
                                    bp_dict[lbp_coord] = []
                                    y = bp_dict[lbp_coord]
                                    y.append(x)

                                #Try to append i-j pair to i-nuc for right i-nuc
                                try:
                                    t = bp_dict[rbp_coord]
                                    t.append(z)
                                #If i-nuc not defined, define it
                                except:
                                    bp_dict[rbp_coord] = []
                                    t = bp_dict[rbp_coord]
                                    t.append(z)
                                l += 2

                    structure = ''.join(structure)
                    if print_to_screen == True:
                        if constraints != None:
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pvalue)+"\t"+str(ED)+"\n"+str(frag)+"\n"+str(window_constraints)+"\n"+str(structure)+"\n"+str(centroid)+str(gc_content)+"\n")
                        else:
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pvalue)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\t"+str(gc_content)+"\n")
                    w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pvalue)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\t"+str(gc_content)+"\n")

                    i += step_size #this ensures that the next iteration increases by "step size" length
                ### Append metrics to list ###
                MFE_list.append(MFE)
                zscore_list.append(zscore)
                pvalue_list.append(pvalue)
                ED_list.append(ED)
            ### Add a final window if the full sequence was not covered ###
            else:
                length = len(seq)
                #print("Sequence Length: "+str(len(seq))+"nt")
                number_windows = ((length+1)-window_size)/step_size
                #print(number_windows)
                if number_windows < float(1.0):
                    number_windows = ((length+1)-window_size)/step_size
                else:
                    number_windows = math.floor(((length+1)-window_size)/step_size)

                remainder = length-number_windows*step_size
                ### 1 here is the start coordinate, change if needed ###
                seqend = 1 + length-1
                #print(remainder, seqend)
                #print(str(int(seqend+1)-window_size), str(seqend))
                if remainder > 0:
                    start_nucleotide = str((seqend)-window_size)
                    end_nucleotide = str(seqend)
                    frag = seq[(int(length))-window_size:int(length)] # This breaks up sequence into fragments
                    if -1 == 0:
                        print("Magic")
                    else:
                        if frag == "N"*window_size:
                            MFE = int(0.0)
                            zscore = "00.00"
                            ED = int(0.0)
                            pscore = int(0.0)
                            structure = "........................................................................................................................"
                            centroid = "........................................................................................................................"
                        else:

                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                            MFE_total.append(float(MFE))
                            ED_total.append(float(ED))
                            seqlist = [] # creates the list we will be filling with sequence fragments
                            seqlist.append(frag) # adds the native fragment to list
                            scrambled_sequences = scramble(frag, randomizations, type)
                            seqlist.extend(scrambled_sequences)
                            energy_list = energies(seqlist, temperature, algo)
                            if print_random == "on":
                                print(energy_list)
                            try:
                                zscore = round(zscore_function(energy_list, randomizations), 2)
                            except:
                                zscore = zscore_function(energy_list, randomizations)
                            zscore_total.append(zscore)
                            pscore = round(pvalue_function(energy_list, randomizations), 2)

            ### Append metrics to list ###
                            MFE_list.append(MFE)
                            zscore_list.append(zscore)
                            pvalue_list.append(pscore)
                            ED_list.append(ED)

                        if print_to_screen == 'on':
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\t"+str(gc_content)+"\n")
                            w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\t"+str(gc_content)+"\n")

            #Define OVERALL values of metrics
            meanz = float(statistics.mean(zscore_total))
            sdz = float(statistics.stdev(zscore_total))
            minz = min(zscore_total)
            stdz = float(statistics.stdev(zscore_total))

            one_sig_below = float(meanz-stdz)
            two_sig_below = float( meanz - ( 2 * stdz) )


        if lri == True:
            #lri_file.write("%s\n" % (duplex_0_range, duplex_1_range, sequence, structure, cofold_zscore, duplex_mfe))
            lri_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("Coordinates(K-mer)", "Coordinates(Duplex)", "Sequence", "Structure", "z-score", "MFE"))
            print("Scanning for long range interactions...")
            j_win = 0
            while j_win == 0 or j_win <= (len(seq) - kmer+1):
                start_nucleotide = j_win# This will just define the start nucleotide coordinate value
                print("\nScanning for k-mer: "+str(j_win+1)+" to "+str(int(j_win)+int(kmer)))
                kmer = int(kmer)
                frag = str(seq[start_nucleotide:start_nucleotide+int(kmer)])
                k_win = 0
                while k_win == 0 or k_win <= (len(cur_record.seq) - kmer):

                    #change this to control how far duplex forms
                    if ((k_win+3) < (start_nucleotide-kmer)) or (k_win > (start_nucleotide+kmer+3)):
                        dup_frag = str(cur_record.seq[k_win:k_win+int(kmer)])
                        duplex = RNA.duplexfold(frag, dup_frag)
                        duplex_mfe = duplex.energy
                        duplex_mfe = round(duplex_mfe, 2)
                        duplex_structure = duplex.structure
                        duplex_structure_split = duplex_structure.split("&")
                        duplex_structure_0 = duplex_structure_split[0]
                        duplex_structure_1 = duplex_structure_split[1]
                        duplex_length_0 = len(duplex_structure_0)
                        duplex_length_1 = len(duplex_structure_1)
                        duplex_0_range = str(str(j_win+1+((duplex.i)-(duplex_length_0)))+"-"+str(j_win+duplex.i))
                        duplex_seq_list = []
                        for coordinate in range(j_win+1+(duplex.i)-(duplex_length_0), j_win+duplex.i+1):
                            duplex_seq_list.append(coordinate)

                        duplex_seq_list.append("&")
                        duplex_1_range = str(str(k_win+(duplex.j))+"-"+str(k_win+duplex.j+duplex_length_1-1))
                        duplex_1_list = []
                        for coordinate in range(k_win+(duplex.j), k_win+duplex.j+duplex_length_1):
                            duplex_seq_list.append(coordinate)

                        duplex_sequence_0 = frag[(duplex.i)-(duplex_length_0):duplex.i]
                        duplex_sequence_1 = dup_frag[(duplex.j)-1:((duplex.j)+(duplex_length_1)-1)]
                        start_sequence_0 = duplex.i-1+start_nucleotide-duplex_length_0
                        start_sequence_1 = duplex.j-1+start_nucleotide
                        sequence = duplex_sequence_0+"&"+duplex_sequence_1

                        duplex_mfe = round(duplex_mfe, 2)
                        if duplex_mfe < lri_cutoff:
                            cofold_seqlist = [] # creates the list we will be filling with sequence fragments
                            cofold_seqlist.append(dup_frag) # adds the native fragment to list
                            scrambled_sequences = scramble(dup_frag, randomizations, type)
                            cofold_seqlist.extend(scrambled_sequences)
                            energy_list = cofold_energies(frag, cofold_seqlist)
                            pvalue_lri = round(pvalue_function(energy_list, randomizations), 2)
                            cofold_zscore = zscore_function(energy_list, randomizations)
                            cofold_zscore = round(cofold_zscore, 2)
                            zscore_total.append(cofold_zscore)
                            if cofold_zscore < 10:
                                ed = 0
                                s = 0
                                structure = duplex_structure
                                print("Hit found for: ")
                                print(duplex_0_range, duplex_1_range, duplex_mfe, sequence, structure)
                                print(cofold_zscore)
                                lri_file.write("%s\t%s\t%s\t%s\t%f\t%f\n" % (duplex_0_range, duplex_1_range, sequence, structure, cofold_zscore, duplex_mfe))

                                #Inititate base pair tabulation variables
                                bond_order = []
                                rev_bond_order = []

                                #Iterate through sequence to assign nucleotides to structure type
                                m = 0
                                bond_count = 0
                                structure = list(duplex_structure)
                                while  m < len(duplex_structure):

                                    if structure[m] == '(':
                                        bond_count += 1
                                        bond_order.append(bond_count)
                                        m += 1
                                    elif structure[m] == ')':
                                        bond_order.append(bond_count)
                                        bond_count -= 1
                                        m += 1
                                    elif structure[m] == '.':
                                        bond_order.append(0)
                                        m += 1
                                    elif structure[m] == '&':
                                        bond_order.append(0)
                                        m += 1
                                    else:
                                        print("Error2")

                                if duplex_seq_list[0] > duplex_seq_list[-1]:
                                    m = 0
                                    bond_count = 0
                                    while  m < len(duplex_structure):
                                        structure = "".join(structure)
                                        structure = flip_structure(duplex_structure)
                                        structure = list(structure)
                                        #print(structure)

                                        if structure[m] == '(':
                                            bond_count += 1
                                            rev_bond_order.append(bond_count)
                                            m += 1
                                        elif structure[m] == ')':
                                            rev_bond_order.append(bond_count)
                                            bond_count -= 1
                                            m += 1
                                        elif structure[m] == '.':
                                            rev_bond_order.append(0)
                                            m += 1
                                        elif structure[m] == '&':
                                            rev_bond_order.append(0)
                                            m += 1
                                        else:
                                            print("Error2")

                                #Initiate base_pair list
                                base_pairs = []
                                rev_base_pairs = []

                                #Create empty variable named test
                                test = ""

                                #Iterate through bond order (downstream)
                                if duplex_seq_list[0] < duplex_seq_list[-1]:
                                    j = 0
                                    while j < len(bond_order):
                                        if bond_order[j] != 0:
                                            test = bond_order[j]
                                            base_pairs.append(j+1)
                                            bond_order[j] = 0
                                            j += 1
                                            bond_k = 0
                                            while bond_k < len(bond_order):
                                                if bond_order[bond_k] == test:
                                                    base_pairs.append(bond_k+1)
                                                    bond_order[bond_k] = 0
                                                    bond_k += 1
                                                else:
                                                    bond_k += 1
                                        else:
                                            j += 1

                                #Iterate through bond order (upstream)
                                if duplex_seq_list[0] > duplex_seq_list[-1]:
                                    j = 0
                                    while j < len(rev_bond_order):
                                        if rev_bond_order[j] != 0:
                                            test = rev_bond_order[j]
                                            rev_base_pairs.append(j+1)
                                            rev_bond_order[j] = 0
                                            j += 1
                                            bond_k = 0
                                            while bond_k < len(bond_order):
                                                if rev_bond_order[bond_k] == test:
                                                    rev_base_pairs.append(bond_k+1)
                                                    rev_bond_order[bond_k] = 0
                                                    bond_k += 1
                                                else:
                                                    bond_k += 1
                                        else:
                                            j += 1

                                #Iterate through "base_pairs" "to define bps
                                l = 0
                                while l < len(base_pairs):

                                    #duplex forms downstream
                                    if duplex_seq_list[0] < duplex_seq_list[-1]:
                                        seq_length = len(sequence)
                                        lbp = base_pairs[l]
                                        rbp = base_pairs[l+1]

                                        lb = str(sequence[int(lbp)-1])
                                        rb = str(sequence[int(rbp)-1])

                                        #lbp_coord = j_win+(duplex.i)-(duplex_length_0)+base_pairs[l]
                                        lbp_coord_key = int(base_pairs[l])-1
                                        lbp_coord = duplex_seq_list[lbp_coord_key]
                                        #print("lbp =", lbp_coord, "lb =", lb)

                                        rbp_coord_key = int(base_pairs[l+1])-1
                                        rbp_coord = duplex_seq_list[rbp_coord_key]
                                        #rbp_coord = k_win+(duplex.j)+base_pairs[l+1]-1
                                        #print("rbp =", rbp_coord, "rb =", rb)
                                        x = NucPair(lb, lbp_coord, rb, rbp_coord, cofold_zscore, duplex_mfe, ed)
                                        z = NucPair(rb, rbp_coord, lb, lbp_coord, cofold_zscore, duplex_mfe, ed)

                                        #Try to append i-j pair to i-nuc for left i-nuc
                                        try:
                                            y = bp_dict[lbp_coord]
                                            y.append(x)
                                        #If i-nuc not defined, define it
                                        except:
                                            bp_dict[lbp_coord] = []
                                            y = bp_dict[lbp_coord]
                                            y.append(x)

                                        l += 2

                                    #duplex forms upstream
                                    elif duplex_seq_list[0] > duplex_seq_list[-1]:

                                        reversed_dup_seq_list = duplex_seq_list[::-1]
                                        reversed_sequence = sequence[::-1]

                                        seq_length = len(reversed_sequence)
                                        lbp = rev_base_pairs[l]
                                        rbp = rev_base_pairs[l+1]

                                        lb = str(reversed_sequence[int(lbp)-1])
                                        rb = str(reversed_sequence[int(rbp)-1])

                                        #lbp_coord = j_win+(duplex.i)-(duplex_length_0)+base_pairs[l]
                                        lbp_coord_key = int(rev_base_pairs[l])-1
                                        lbp_coord = reversed_dup_seq_list[lbp_coord_key]
                                        #print("lbp =", lbp_coord, "lb =", lb)

                                        rbp_coord_key = int(rev_base_pairs[l+1])-1
                                        rbp_coord = reversed_dup_seq_list[rbp_coord_key]
                                        #rbp_coord = k_win+(duplex.j)+base_pairs[l+1]-1
                                        # print("rbp =", rbp_coord, "rb =", rb, "lbp =", lbp_coord, "lb =", lb)
                                        x = NucPair(lb, lbp_coord, rb, rbp_coord, cofold_zscore, duplex_mfe, ed)
                                        z = NucPair(rb, rbp_coord, lb, lbp_coord, cofold_zscore, duplex_mfe, ed)
                                        # print(x.inucleotide, x.jnucleotide, x.icoordinate, x.jcoordinate)


                                        if lbp_coord > len(cur_record.seq):
                                            print("ERROR", x.inucleotide, x.jnucleotide, x.icoordinate, x.jcoordinate)

                                        if rbp_coord > len(cur_record.seq):
                                            print("ERROR", x.inucleotide, x.jnucleotide, x.icoordinate, x.jcoordinate)


                                        #Try to append i-j pair to i-nuc for left i-nuc
                                        try:
                                            y = bp_dict[lbp_coord]
                                            y.append(x)
                                        #If i-nuc not defined, define it
                                        except:
                                            bp_dict[lbp_coord] = []
                                            y = bp_dict[lbp_coord]
                                            y.append(x)

                                        if int(lbp_coord) == int(rbp_coord):
                                            print("WARNING\n")

                                        l += 2

                                    else:
                                        print("Error", duplex_seq_list[0], duplex_seq_list[-1])
                        k_win += kmer_step_size

                    else:
                        k_win += kmer_step_size

                j_win += kmer_step_size

            #Define OVERALL values of metrics
            meanz = float(statistics.mean(zscore_total))
            sdz = float(statistics.stdev(zscore_total))
            minz = min(zscore_total)
            stdz = float(statistics.stdev(zscore_total))

            one_sig_below = float(meanz-stdz)
            two_sig_below = float( meanz - ( 2 * stdz) )
            print("LRI Scan Complete")

        best_bps = {}
        best_sum_bps = {}
        best_sum_bps_means = {}
        best_total_window_mean_bps = {}

        #Determine start and end coordinate values
        start_coordinate = str(list(nuc_dict.keys())[0])
        #print(start_coordinate)
        end_coordinate = str(list(nuc_dict.keys())[-1])
        #print(end_coordinate)

        #Iterate through initial i-nuc dictionary to determine best base pairs (round 1)
        elapsed_time = round((time.time() - start_time), 2)
        print("Elapsed time: "+str(elapsed_time)+"s")
        print("Determining best base pairs...")
        for k, v in sorted(bp_dict.items()):
            #Initiate local dictionaries to store metrics per nucleotide
            zscore_dict = {}
            pair_dict = {}
            mfe_dict = {}
            ed_dict = {}

            #Iterate through all i-j pairs per i-nucleotide to store metrics for each
            for pair in v:
                #Create a key  which contains nucleotide and coordinate info
                partner_key = str(pair.jnucleotide)+"-"+str(pair.jcoordinate)
                #print(partner_key)

                #Create a variable which contains all i-j pair info
                x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                            pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)

                #Try to append the value of each metric to metric lists per i-nuc
                try:
                    y = zscore_dict[partner_key]
                    y.append(pair.zscore)

                    m = mfe_dict[partner_key]
                    m.append(pair.mfe)

                    e = ed_dict[partner_key]
                    e.append(pair.ed)

                    z = pair_dict[partner_key]
                    z.append(x)

                #If pair not defined, define it
                except:
                    zscore_dict[partner_key] = []
                    y = zscore_dict[partner_key]
                    y.append(pair.zscore)
                    pair_dict[partner_key] = []

                    mfe_dict[partner_key] = []
                    m = mfe_dict[partner_key]
                    m.append(pair.mfe)

                    ed_dict[partner_key] = []
                    e = ed_dict[partner_key]
                    e.append(pair.ed)

                    z = pair_dict[partner_key]
                    z.append(x)


            #Caclulate and store sum of z-score per i-j pair
            sum_z = {}
            sum_z_lengths = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = sum(v1)
                test = sum_z[k1] = sum(v1)
                sum_z_lengths[k1] = len(sum_z)

            #Caclulate and store sum of ED per i-j pair
            sum_ed = {}
            sum_ed_lengths = {}
            for k1, v1 in zscore_dict.items():
                sum_ed[k1] = sum(v1)
                test = sum_ed[k1] = sum(v1)
                sum_ed_lengths[k1] = len(sum_ed)



            #Caclulate and store mean of z-score per i-j pair
            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = statistics.mean(v1)
                test = mean_z[k1] = statistics.mean(v1)

            #Caclulate and store mean MFE per i-j pair
            mean_mfe = {}
            for k1, v1 in mfe_dict.items():
                mean_mfe[k1] = statistics.mean(v1)

            #Caclulate and store mean ED per i-j pair
            mean_ed = {}
            for k1, v1 in ed_dict.items():
                mean_ed[k1] = statistics.mean(v1)

            #Caclulate and store total window counts per i-j pair
            total_windows = 0
            num_bp = 0
            for k1, v1 in zscore_dict.items():
                total_windows = total_windows + len(v1)
                #key_data = re.split("-", str(k1))
                key_i = (k1[2:])
                if int(k) == int(key_i):
                    continue
                if int(k) != int(key_i):
                    num_bp += 1

            #Print first line of log file tables (first half of log file)
            k_nuc = str((nuc_dict[k].nucleotide))
            #print(k_nuc)
            if by_ed == False:
                log_total.write("\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED"
                      +"\tSumZ\tSumZ/#TotalWindows\tBPs= "+str(num_bp)+"\n")
                log_total.write("nt-"+str(k)+"\t-\t"+str(k_nuc)+"\t"+str(total_windows)
                      +"\t-\t-\t-\t-\t-"+"\n")
            if by_ed == True:
                log_total.write("\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED"
                      +"\tSumED\tSumED/#TotalWindows\tBPs= "+str(num_bp)+"\n")
                log_total.write("nt-"+str(k)+"\t-\t"+str(k_nuc)+"\t"+str(total_windows)
                      +"\t-\t-\t-\t-\t-"+"\n")
            sirna_log.write(f"{k}\t{k_nuc}\t{total_windows}\t{num_bp}\n")
            #Print remainder of log file tables (first half of log file)
            if by_ed == False:
                total_window_mean_z = {}
                for k1, v1 in zscore_dict.items():
                    bp_window = str(len(v1))
                    #key_data = re.split("-", str(k1))
                    key_nuc = k1[:1]
                    key_i = int(k1[2:])
                    total_window_mean_z[k1] = (sum(v1))/total_windows
                    z_sum = str(round(sum(v1), 2))
                    z_avg = str(round(statistics.mean(v1), 2))
                    test = str(round(total_window_mean_z[k1], 2))
                    k1_mean_mfe = str(round(mean_mfe[k1], 2))
                    k1_mean_ed = str(round(mean_ed[k1], 2))
                    if int(k) == int(key_i):
                        #print("iNuc is "+str(key_i))
                        log_total.write(str(k)+"\tNoBP\t"+key_nuc+"\t"+bp_window+"\t"
                              +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                              +z_sum+"\t"+str(test)+"\n")
                    else:
                        #print("j is "+str(k))
                        log_total.write(str(k)+"\t"+str(key_i)+"\t"+str(key_nuc)+"\t"+bp_window+"\t"
                              +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                              +z_sum+"\t"+str(test)+"\n")

                    #Define best_bp_key based on coverage-normalized z-score
                    best_bp_key = min(total_window_mean_z, key = total_window_mean_z.get)


            if by_ed == True:
                total_window_mean_ed = {}
                for k1, v1 in ed_dict.items():
                    bp_window = str(len(v1))
                    key_data = re.split("-", str(k1))
                    key_nuc = str(key_data[0])
                    key_i = str(key_data[1])
                    total_window_mean_ed[k1] = (sum(v1))/total_windows
                    ed_sum = str(round(sum(v1), 2))
                    ed_avg = str(round(statistics.mean(v1), 2))
                    test = str(round(total_window_mean_ed[k1], 2))
                    k1_mean_mfe = str(round(mean_mfe[k1], 2))
                    k1_mean_z = str(round(mean_z[k1], 2))
                    if int(k) == int(key_i):
                        #print("iNuc is "+str(key_i))
                        log_total.write(str(k)+"\tNoBP\t"+key_nuc+"\t"+bp_window+"\t"
                              +k1_mean_mfe+"\t"+k1_mean_z+"\t"+ed_avg+"\t"
                              +ed_sum+"\t"+str(test)+"\n")
                    else:
                        #print("j is "+str(k))
                        log_total.write(str(k)+"\t"+key_i+"\t"+key_nuc+"\t"+bp_window+"\t"
                              +k1_mean_mfe+"\t"+k1_mean_z+"\t"+ed_avg+"\t"
                              +ed_sum+"\t"+str(test)+"\n")

                    #Define best_bp_key based on coverage-normalized ED
                    best_bp_key = min(total_window_mean_ed, key = total_window_mean_ed.get)

            #Access best i-j NucPairs for each metric using best_bp_key
            if by_ed == False:
                best_bp_mean_z = mean_z[best_bp_key]
                best_bp_sum_z = sum_z[best_bp_key]
                best_bp_mean_mfe = mean_mfe[best_bp_key]
                best_bp_mean_ed = mean_ed[best_bp_key]
                best_total_window_mean_z = total_window_mean_z[best_bp_key]

            if by_ed == True:
                best_bp_mean_ed = mean_ed[best_bp_key]
                best_bp_sum_ed = sum_ed[best_bp_key]
                best_bp_mean_mfe = mean_mfe[best_bp_key]
                best_bp_mean_z = mean_z[best_bp_key]
                best_total_window_mean_ed = total_window_mean_ed[best_bp_key]


            #Access best i-j pair info from key name
            best_bp_data = re.split("-", best_bp_key)
            best_nucleotide = best_bp_data[0]
            best_coordinate = best_bp_data[1]

            #Fill dictionary with coverage normalized z-score
            #print("Determining best base pair for nucleotide ", k)
            if by_ed == False:
                best_total_window_mean_bps[k] = (NucPair((nuc_dict[k]).nucleotide,
                                                 nuc_dict[k].coordinate, best_nucleotide,
                                                 best_coordinate, best_total_window_mean_z,
                                                 best_bp_mean_mfe, best_bp_mean_ed))

                #Fill dictionary with coverage average z-score
                best_bps[k] = (NucPair((nuc_dict[k]).nucleotide, (nuc_dict[k]).coordinate,
                                        best_nucleotide, best_coordinate, best_bp_mean_z,
                                        best_bp_mean_mfe, best_bp_mean_ed))

            if by_ed == True:
                best_total_window_mean_bps[k] = (NucPair((nuc_dict[k]).nucleotide,
                                                 nuc_dict[k].coordinate, best_nucleotide,
                                                 best_coordinate, best_total_window_mean_ed,
                                                 best_bp_mean_mfe, best_bp_mean_ed))

                #Fill dictionary with coverage average z-score
                best_bps[k] = (NucPair((nuc_dict[k]).nucleotide, (nuc_dict[k]).coordinate,
                                        best_nucleotide, best_coordinate, best_bp_mean_z,
                                        best_bp_mean_mfe, best_bp_mean_ed))

        ######## Detect competing partners, and select final i-j pairs #################
        final_partners = {}
        elapsed_time = round((time.time() - start_time), 2)
        print("Elapsed time: "+str(elapsed_time)+"s")

        #print header for fianl partener log file (log_win)
        log_win.write("i\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED"
            + "\t*Indicates most favorable bp has competition; bp(j) has more favorable partner or is "
            + "more likely to be unpaired"+"\n")

        #Iterate through round 1 i-j pairs
        if competition == 1:
            #print(start_coordinate, end_coordinate)
            print("Detecting competing pairs...")
            j_coord_list = []
            # for k, v in sorted(best_bps.items()):
            #     print(jcoordinate)
            #     j_coord_list.append(int(v.jcoordinate))

            for k, v in sorted(best_bps.items()):
                #print(k, v.icoordinate, v.jcoordinate)
                test_k = int(k)
                #print(sum(test_k == int(v.jcoordinate) for v in best_bps.values()))
                #print(v.jcoordinate)
                if sum(test_k == int(v.jcoordinate) for v in best_bps.values()) >= 0:
                    #print(start_coordinate, end_coordinate)
                ### Scan the entire dictionary:
                    # keys = range(int(start_coordinate), int(end_coordinate))

                ### Scan two window's length flanking nucleotide:
                    #print(len(best_bps))
                    # print(length*4)
                    length = len(cur_record.seq)
                    #print(length)
                    if (len(best_bps) < length*4):
                        #print("Scanning full dictionary")
                        #Length of input less than length of flanks
                        # keys = range(int(start_coordinate), int(end_coordinate))
                        subdict = best_total_window_mean_bps

                    elif (
                        (v.icoordinate - length*(2)) >= int(start_coordinate) and
                        (v.icoordinate + (length*2)) <= int(end_coordinate)
                        ):
                        #print(str(v.icoordinate - length*(2)))
                        # print("MIDDLE")
                        keys = range(int(v.icoordinate-(length*2)), int(v.icoordinate+(length*2)))
                        subdict = {k: best_total_window_mean_bps[k] for k in keys}

                    elif (
                        int(v.icoordinate + (length*(2))) <= int(end_coordinate)and
                        (v.icoordinate + (length*2)) <= int(end_coordinate)
                    ):
                        #print("BEGINING"+str(v.icoordinate - (length*(2)))+" "+str(end_coordinate))
                        keys = range(int(start_coordinate), int(v.icoordinate+(length*2))+1)
                        subdict = {k: best_total_window_mean_bps[k] for k in keys}

                    elif (v.icoordinate + (length*2)) >= int(end_coordinate):
                        if v.icoordinate-(length*2) > 0:
                            #print("END"+str(v.icoordinate + (length*2)))
                            keys = range(int(v.icoordinate-(length*2)), int(end_coordinate))
                        else:
                            keys =range(int(v.icoordinate-(length*2)), int(end_coordinate))
                            subdict = {k: best_total_window_mean_bps[k] for k in keys}

                    elif len(best_bps) < length:
                            subdict = best_total_window_mean_bps

                    else:
                        print("Sub-dictionary error")
                        raise ValueError("Sub-dictionary error")

                    #print(keys)
                    # # if k == 216:
                    # #     for subk, subv in subdict.items():
                    # #         print(subk, subv.icoordinate, subv.jcoordinate)
                    # #print("SubDict length for "+str(k)+"="+str(len(subdict)))
                    #
                    if len(subdict) >= 0:

                        #print("Found competing pair for "+str(k))
                        #elapsed_time = round((time.time() - start_time), 2)
                        #print(elapsed_time)
                        #print("Detecting competing pairs for nuc ", k)
                        #For each i and j in i-j pair, detect competing pairs and append to dict
                        comp_pairs_i = competing_pairs(subdict, v.icoordinate)
                        #print("i-pairs="+str(len(comp_pairs_i)))
                        comp_pairs_j = competing_pairs(subdict, v.jcoordinate)
                        #print("j-pairs="+str(len(comp_pairs_j)))
                        total_pairs = []
                        #Put pairs competing with i from i-j pair into total pair dict for i-nuc
                        for key, pair in comp_pairs_i.items():
                            #print("checking competing pairs for i")
                            # if k == 216:
                            #     print(k, pair.icoordinate, pair.jcoordinate, pair.zscore)
                            total_pairs.append(competing_pairs(subdict,
                                                               pair.jcoordinate))
                            total_pairs.append(competing_pairs(subdict,
                                                               pair.icoordinate))
                        #Put pairs competing with j from i-j pair into total pair dict for i-nuc
                        for key, pair in comp_pairs_j.items():
                            #print("checking competing pairs for j")
                            # if k == 216:
                            #     print(k, pair.icoordinate, pair.jcoordinate, pair.zscore)
                            total_pairs.append(competing_pairs(subdict,
                                                               pair.jcoordinate))
                            total_pairs.append(competing_pairs(subdict,
                                                               pair.icoordinate))
                        # print(str(k)+"nt Total comp pairs="+str(len(total_pairs)))

                        #Merge all dictionaries
                        merged_dict = {}
                        i = 0
                        for d in total_pairs:
                            #print("merging competing dictionaries "+str(i))
                            for k1, v1 in d.items():
                                # if k == 216:
                                #     print(k, k1, v1.icoordinate, v1.jcoordinate, v1.zscore)
                                merged_dict[i] = v1
                                i += 1

                        # #print("MergedDict length for "+str(k)+"="+str(len(merged_dict)))
                        # #initiate best_basepair fucntion, return best_bp based on sum
                        # if len(merged_dict) > 2:
                        #     bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                        #     #print(str(len(merged_dict))+"__11111")
                        # else:
                        #     #print("Nucleotide "+str(k))
                        #     bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                        #     #print(str(len(merged_dict))+"____222222")
                        #     #bp = best_total_window_mean_bps[k]
                        # #Check if best basepair was connected to i-nucleotide (i.e., "k")
                        # print(len(merged_dict))
                        if len(merged_dict) > 0:
                            bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")
                        else:
                            bp = NucPair(v.inucleotide, v.icoordinate,
                                                        v.inucleotide, v.icoordinate,
                                                        v.zscore,
                                                        v.mfe,
                                                        v.ed)

                        if (int(k) != bp.icoordinate) and (int(k) != int(bp.jcoordinate)):
                            # print("1 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                            #if there was a competing i-j pair print it to log file instead:
                            log_win.write("nt-"+str(k)+"*:\t"+str(v.icoordinate)+"\t"+v.jcoordinate+"\t"
                                  +str(round(v.mfe, 2))
                                  +"\t"+str(round(v.zscore, 2))
                                  +"\t"+str(round(v.ed, 2))+"\n")
                            final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                                        v.inucleotide, v.icoordinate,
                                                        best_bps[bp.icoordinate].zscore,
                                                        bp.mfe,
                                                        bp.ed)
                        #
                        # elif ((int(v.icoordinate) == int(v.jcoordinate)) and (int(bp.icoordinate) != int(bp.jcoordinate))):
                        #     #Check for instance where competing base pair
                        #     print("!!!!!!!2 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                        #     log_win.write("nt-"+str(k)+"*:\t"+str(bp.icoordinate)+"\t"+bp.jcoordinate+"\t"
                        #           +str(round(bp.mfe, 2))
                        #           +"\t"+str(round(bp.zscore, 2))
                        #           +"\t"+str(round(bp.ed, 2))+"\n")
                        #
                        #     final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                        #                                 bp.jnucleotide, bp.jcoordinate,
                        #                                 best_bps[bp.icoordinate].zscore,
                        #                                 best_bps[bp.icoordinate].mfe,
                        #                                 best_bps[bp.icoordinate].ed)
                        #
                        #
                        else:
                            #print("3 = "+str(v.icoordinate)+"_"+str(v.jcoordinate)+" AND "+str(bp.icoordinate)+"_"+str(bp.jcoordinate))
                            log_win.write("nt-"+str(k)+":\t"+str(bp.icoordinate)+"\t"+str(bp.jcoordinate)+"\t"
                                  + str(round(best_bps[k].mfe, 2))+"\t"
                                  + str(round(best_bps[k].zscore, 2))
                                  + "\t"+str(round(best_bps[k].ed, 2))+"\n")
                            final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                                                        bp.jnucleotide, bp.jcoordinate,
                                                        best_bps[bp.icoordinate].zscore,
                                                        best_bps[bp.icoordinate].mfe,
                                                        best_bps[bp.icoordinate].ed)

                    else:
                        continue
                else:
                    final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                                v.jnucleotide, v.jcoordinate,
                                                best_bps[k].zscore,
                                                best_bps[k].mfe,
                                                best_bps[k].ed)
                    #print("No competing pair found for ", k)
                    continue
        if competition == 0:
            elapsed_time = str(round((time.time() - start_time), 2))+"s"
            print("Elapsed time: "+elapsed_time)
            print("Writing DP files, can not write CT files...")
            write_bp(best_bps, out6+"."+outname+cur_record.name+".ALL.bp", start_coordinate, name, minz)
            if filter != None:
                write_dp(best_bps, output+str(filter)+cur_record.name+".dp", filter, minz)
            write_dp(best_bps, out1+"."+outname+".dp", float(10), minz)
            write_dp(best_bps, out2+"."+outname+".dp", float(-1), minz)
            write_dp(best_bps, out3+"."+outname+".dp", float(-2), minz)
            write_dp(best_bps, output+"."+outname+"mean_"+str(round(meanz, 2))+".dp", meanz, minz)
            write_dp(best_bps, output+"."+outname+"below_mean_"+str(round(one_sig_below, 2))+".dp", one_sig_below, minz)
            print("ScanFold-Fold complete, find results in...")

        #Write CT files
        if competition == 1:
            print("Trying to write CT files with -c option")
            elapsed_time = str(round((time.time() - start_time), 2))+"s"
            print("Elapsed time: "+elapsed_time)
            print("Writing CT files")
            #print(len(final_partners))

            # if filter != None or filter != -2:
            #     (final_partners, output+str(filter)+".ct", filter, strand, name)
            # write_ct(final_partners, output+"no_filter.ct", float(10), strand, name, start_coordinate)
            # write_ct(final_partners, output+"-1.ct", float(-1), strand, name, start_coordinate)
            # write_ct(final_partners, output+"-2.ct", float(-2), strand, name, start_coordinate)
            # if filter != None or filter != -2:
            #     write_ct(final_partners, output+str(filter)+".ct", filter, strand, name, start_coordinate)
            strand = 1
            write_ct(final_partners, str(dbn_file_path1)+".ct", float(10), strand, name, start_coordinate)
            write_ct(final_partners, str(dbn_file_path2)+".ct", float(-1), strand, name, start_coordinate)
            write_ct(final_partners, str(dbn_file_path3)+".ct", float(-2), strand, name, start_coordinate)
            makedbn(dbn_file_path1, "NoFilter")
            makedbn(dbn_file_path2, "Zavg_-1")
            makedbn(dbn_file_path3, "Zavg_-2")
            write_bp(final_partners, out6+"."+outname+".bp", start_coordinate, name, minz)
            write_wig_dict(final_partners, final_partners_wig+"."+outname+".wig", name, step_size)
            write_bp(best_bps, out6+"."+outname+"."+cur_record.name+".ALL.bp", start_coordinate, name, minz)

        write_fasta(nuc_dict, name+"."+outname+".fa", name)

        if lri == False:
            write_wig(MFE_list, step_size, name, outname+".scan-MFE.wig")
            write_wig(zscore_list, step_size, name, outname+".scan-zscores.wig")
            write_wig(pvalue_list, step_size, name, outname+".scan-pvalue.wig")
            write_wig(ED_list, step_size, name, outname+".scan-ED.wig")
        full_fasta_sequence = seq

        #If specified, refold and generate global fold of input sequence
        #Create a dbn file for forna
        # os.system(str("ct2dot "+str(dbn_file_path1)+".ct 1 "+str(dbn_file_path1)))
        # os.system(str("ct2dot "+str(dbn_file_path2)+".ct 1 "+str(dbn_file_path2)))
        # os.system(str("ct2dot "+str(dbn_file_path3)+".ct 1 "+str(dbn_file_path3)))

        if global_refold == True:
            #create a separate DBN file
            dbn_log_file = open(str(dbn_file_path), "w+")

            #try:
            #fold the full fasta input as a fold compound (full_fc) using model params (md)
            print("Refolding full sequence using ScanFold results as constraints...")
            elapsed_time = round((time.time() - start_time), 2)
            print("Elapsed time: "+str(elapsed_time)+"s")
            md = RNA.md()
            md.temperature = int(temperature)

            #refold from -1 constraints
            fc = RNA.fold_compound(str(full_fasta_sequence), md)
            dbn_file_filter1 = open(str(dbn_file_path2+".dbn"), "r")
            lines = dbn_file_filter1.readlines()
            filter1constraints = str(lines[2])
            refolded1filter = fc.hc_add_from_db(filter1constraints)
            (refolded_filter1_structure, refolded_filter1_MFE) = fc.mfe()

            #refold from -2 constraints
            fc = RNA.fold_compound(str(full_fasta_sequence), md)
            dbn_file_filter2 = open(str(dbn_file_path3+".dbn"), "r")
            lines = dbn_file_filter2.readlines()
            filter2constraints = str(lines[2])
            refolded2filter = fc.hc_add_from_db(filter2constraints)
            (refolded_filter2_structure, refolded_filter2_MFE) = fc.mfe()

            #extract the structure
            full_fc = RNA.fold_compound(str(full_fasta_sequence), md)
            (full_structure, full_MFE) = full_fc.mfe()

            dbn_log_file.write(">"+str(name)+"\tGlobal Full MFE="+str(full_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(full_structure)+"\n")
            dbn_log_file.write(">"+str(name)+"\tRefolded with -1 constraints MFE="+str(refolded_filter1_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(refolded_filter1_structure)+"\n")
            dbn_log_file.write(">"+str(name)+"\tRefolded with -2 constraints MFE="+str(refolded_filter2_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(refolded_filter2_structure)+"\n")
            dbn_log_file.close()
            os.system(str("cat "+str(dbn_file_path1+".dbn")+" "+str(dbn_file_path2+".dbn")+" "+str(dbn_file_path3+".dbn")+" > "+str(dbn_file_path4)))
            dbn_log_file.close()
            # except:
            #     print("Automatic refold for "+cur_record.name+" failed. Run manually")
            #     os.chdir(str(original_directory))
            #     del final_partners
            #     continue

        if global_refold == False:
            os.system(str("cat "+str(dbn_file_path1+".dbn")+" "+str(dbn_file_path2+".dbn")+" "+str(dbn_file_path3+".dbn")+" > "+str(dbn_file_path4)))


        #############
        #Begin the structure extract process
        #try: #allow a bypass is the process fails.
        #Set flanking nucleotides to be folded
        flanking = 0

        #Set number of randomizations and shuffle type ("mono" or "di")
        sub_randomizations = 100
        sub_shuffle = "mono"


        #Inititate base pair tabulation variables
        bond_order = []
        bond_count = 0

        #Read the structure of -2 filter2constraints
        if global_refold == False:
            #refold from -2 constraints
            # try:
            dbn_file_filter2 = open(dbn_file_path3+".dbn", "r")
            lines = dbn_file_filter2.readlines()
            #print(lines)
            filter2constraints = str(lines[2])
            full_fasta_sequence = str(lines[1])

        structure_raw = filter2constraints

        sequence = list(str(seq))
        #print(sequence)
        structure = list(structure_raw)
        #print(structure)
        length = len(cur_record.seq)
        #print("sequence length "+str(length)+"nt")
        length_st = len(structure)
        #print(length_st)


        #Iterate through sequence to assign nucleotides to structure type
        m = 0
        nuc_dict_refold = {}
        while  m < len(structure)-1:
            #print(m)
            if structure[m] == '(':
                #print(m, structure[m])
                bond_count += 1
                bond_order.append(bond_count)
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif structure[m] == ')':
                # print(m, structure[m])
                bond_order.append(bond_count)
                bond_count -= 1
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1

            elif str(structure[m]) == ( '.' ):
                # print(m, structure[m])
                bond_order.append(0)
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '<' ):
                # print(m, structure[m])
                bond_order.append(0)
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '>' ):
                # print(m, structure[m])
                bond_order.append(0)
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '{' ):
                # print(m, structure[m])
                bond_order.append(0)
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            elif str(structure[m]) == ( '}' ):
                # print(m, structure[m])
                bond_order.append(0)
                nuc_dict_refold[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
                m += 1
            else:
                # print("Error", bond_count, (m+1), sequence[m], structure[m])
                m += 1
                continue
                # print("no")

        #Initiate base_pair list
        base_pairs = []

        #Create empty variable named test
        test = ""

        #Iterate through bond order
        j = 0
        structure_count = 0
        structure_end = []
        structure_start = []
        while j < len(structure)-1:
            #print(nuc_dict[j].structure)
            try:
                if (nuc_dict_refold[j].bond_order == 1) and (nuc_dict_refold[j].structure == '('):
                    structure_count += 1
                    #print(nuc_dict[j].structure)
                    structure_start.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                    j += 1

                elif (nuc_dict_refold[j].bond_order == 0) and (nuc_dict_refold[j].structure == ')'):
                    structure_count += 1
                    #print(nuc_dict[j].structure)
                    structure_end.append(NucStructure(structure_count, nuc_dict_refold[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                    j += 1
                elif (nuc_dict_refold[j].bond_order == 1) and (nuc_dict_refold[j].structure == '<'):
                    structure_count += 1
                    #print(nuc_dict[j].structure)
                    structure_start.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                    j += 1

                elif (nuc_dict_refold[j].bond_order == 0) and (nuc_dict_refold[j].structure == '>'):
                    structure_count += 1
                    #print(nuc_dict[j].structure)
                    structure_end.append(NucStructure(structure_count, nuc_dict_refold[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                    j += 1
                elif (nuc_dict_refold[j].bond_order == 1) and (nuc_dict_refold[j].structure == '{'):
                    structure_count += 1
                    #print(nuc_dict[j].structure)
                    structure_start.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                    j += 1

                elif (nuc_dict_refold[j].bond_order == 0) and (nuc_dict_refold[j].structure == '}'):
                    structure_count += 1
                    #print(nuc_dict[j].structure)
                    structure_end.append(NucStructure(structure_count, nuc_dict_refold[j].coordinate, nuc_dict_refold[j].nucleotide, nuc_dict_refold[j].structure))
                    j += 1
                else:
                    j += 1
            except:
                print(j, "EXCEPT")
                j += 1
                continue

        l = 0
        extracted_structure_list = []
        #print(structure_start)
        #print(structure_end)
        while l < int(len(structure_start)):
            #print(len(structure_start), len(structure_end))
            offset = flanking = 0
            s = structure_start_coordinate =  int((structure_start[l].coordinate)-offset)
            e = structure_end_coordinate = int((structure_end[l].coordinate)+offset-1)

            seq = ""
            se_fold = ""
            for k, v in nuc_dict_refold.items():
                if s <= k <= e:
                    seq += str(v.nucleotide)
                    se_fold += str(v.structure)

            extracted_structure_list.append(ExtractedStructure(l, seq, se_fold, s, e))
            l += 1


        zscore_total = []
        numerical_z = []
        pvalue_total = []
        numerical_p = []
        MFE_total = []
        ED_total = []

        se = open(structure_extract_file, "w")
        #se.write("ScanFold predicted structures which contain at least one base pair with Zavg < -1 have been extracted from "+str(name)+" results (sequence length "+str(length)+"nt) and have been refolded using RNAfold to determine their individual MFE, structure, z-score (using 100X randomizations), and ensemble diversity score.\n")
        motif_num = 1
        for es in extracted_structure_list[:]:
        #try:
            #print(str(i))
            frag = es.sequence
            fc = RNA.fold_compound(str(frag)) #creates "Fold Compound" object
            fc.hc_add_from_db(str(es.structure))
            fc.pf() # performs partition function calculations
            frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
            (MFE_structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
            MFE = round(MFE, 2)
            MFE_total.append(MFE)
            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
            ED_total.append(ED)
            #print(structure)
            #fmfe = fc.pbacktrack()
            #print(str(fmfe))
            seqlist = [] # creates the list we will be filling with sequence fragments
            seqlist.append(frag) # adds the native fragment to list
            scrambled_sequences = scramble(frag, 100, type)
            seqlist.extend(scrambled_sequences)
            energy_list = energies(seqlist, temperature, algo)
            try:
                zscore = round(zscore_function(energy_list, 100), 2)
            except:
                zscore = zscore_function(energy_list, 100)
            zscore_total.append(zscore)

            pvalue = round(pvalue_function(energy_list, 100), 2)
            #print(pvalue)
            pvalue_total.append(pvalue)
            accession = str(name)
            ps_title = f"motif_{motif_num} coordinates {es.i} - {es.j}"
            #print(macro)

            with open(f"{name}_motif_{motif_num}.dbn", 'w') as es_dbn:
                es_dbn.write(f">{name}_motif_{motif_num}_coordinates:{es.i}-{es.j}\n{frag}\n{MFE_structure}")
            dbn2ct(f"{name}_motif_{motif_num}.dbn")
            ###Create figure using forgi
            # with open('tmp.fa', 'w') as file:
            #     file.write('>tmp\n%s\n%s' % (frag, MFE_structure))
            # cg = forgi.load_rna("tmp.fa", allow_many=False)
            # fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7, backbone_kwargs={"linewidth":3})
            # plt.savefig(ps_title+".png")
            # plt.clf()


            RNA.PS_rna_plot_a(frag, MFE_structure, "motif_"+str(motif_num)+".ps", '', '')
            gff_attributes = f'motif_{motif_num};sequence={es.sequence};structure={str(es.structure)};refoldedMFE={str(MFE_structure)};MFE(kcal/mol)={str(MFE)};z-score={str(zscore)};ED={str(ED)}'
            se.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str("."), str("RNA_sequence_secondary_structure"), str(int((es.i)+1)), str(int((es.j)+1)), str("."), str("."), str("."), gff_attributes))
            motif_num += 1

        se.close()
        # except:
        #     print("Structure Extract failed for "+folder_name+", must extract manually.")
        #     continue
        #print("TEST")
        elapsed_time = round((time.time() - start_time), 2)
        print("Total runtime: "+str(elapsed_time)+"s")
        print("ScanFold-Fold analysis complete! Output found in folder named: "+folder_name)

        os.chdir(str(original_directory))
        del final_partners
