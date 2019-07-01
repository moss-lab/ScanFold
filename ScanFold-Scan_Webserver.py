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
$ ScanFold-Scan_Webserver.py -i fasta.fa --scan_out_path ./out.tsv --zscore_wig_file_path ./zscore.bw --mfe_wig_file_path ./mfe.bw --ed_wig_file_path ./ed.bw --pvalue_wig_file_path ./pvalue.bw --fasta_file_path ./fasta.fasta.fa --fasta_index ./fasta.fai --nodeid "/scholar" --callbackurl "https://www.google.com"
"""

import time
import sys
import argparse
import string
# not needed anymore
# import pyBigWig
import re
import numpy as np
sys.path.append('/home/randrews/ViennaRNA/lib/python3.6/site-packages/')
sys.path.append('/usr/local/lib/python3.6/site-packages')
import RNA
import random
import multiprocessing
import requests

# temporary fix to disable warning about ssl
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from Bio import SeqIO
# from progressbar import *               # just a simple progress bar

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

### input parms ###

parser.add_argument('--nodeid', type=str,
                    help='node id')
parser.add_argument('--callbackurl', type=str,
                    help='callbackurl')


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
nodeid = args.nodeid
callbackurl = args.callbackurl

split = args.split

#### Defining global variables ###############

### The RNA model "md" needs to be set ABOVE FUNCTIONS! ####
### DO NOT MOVE ###
md = RNA.md()
md.temperature = int(temperature)
### DO NOT MOVE ###
### The RNA model "md" needs to be set ABOVE FUNCTIONS! ####

def multiprocessing(func, args,
                    workers):
    with ProcessPoolExecutor(workers) as ex:
        res = ex.map(func, args)
    return list(res)

#### Defining Dinucleotide function #####
# Taken from
# altschulEriksonDinuclShuffle.py
# P. Clote, Oct 2003
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only nonoverlapping words, presumably in a left to
# right fashion.
def computeCountAndLists(s):
  #WARNING: Use of function count(s,'UU') returns 1 on word UUU
  #since it apparently counts only nonoverlapping words UU
  #For this reason, we work with the indices.

  #Initialize lists and mono- and dinucleotide dictionaries
  List = {} #List is a dictionary of lists
  List['A'] = []; List['C'] = [];
  List['G'] = []; List['U'] = [];
  nuclList   = ["A","C","G","U"]
  s       = s.upper()
  s       = s.replace("T","U")
  nuclCnt    = {}  #empty dictionary
  dinuclCnt  = {}  #empty dictionary
  for x in nuclList:
    nuclCnt[x]=0
    dinuclCnt[x]={}
    for y in nuclList:
      dinuclCnt[x][y]=0

  #Compute count and lists
  nuclCnt[s[0]] = 1
  nuclTotal     = 1
  dinuclTotal   = 0
  for i in range(len(s)-1):
    x = s[i]; y = s[i+1]
    List[x].append( y )
    nuclCnt[y] += 1; nuclTotal  += 1
    dinuclCnt[x][y] += 1; dinuclTotal += 1
  assert (nuclTotal==len(s))
  assert (dinuclTotal==len(s)-1)
  return nuclCnt,dinuclCnt,List

def chooseEdge(x,dinuclCnt):
  numInList = 0
  for y in ['A','C','G','U']:
    numInList += dinuclCnt[x][y]
  z = random.random()
  denom=dinuclCnt[x]['A']+dinuclCnt[x]['C']+dinuclCnt[x]['G']+dinuclCnt[x]['U']
  numerator = dinuclCnt[x]['A']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['A'] -= 1
    return 'A'
  numerator += dinuclCnt[x]['C']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['C'] -= 1
    return 'C'
  numerator += dinuclCnt[x]['G']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['G'] -= 1
    return 'G'
  dinuclCnt[x]['U'] -= 1
  return 'U'

def connectedToLast(edgeList,nuclList,lastCh):
  D = {}
  for x in nuclList: D[x]=0
  for edge in edgeList:
    a = edge[0]; b = edge[1]
    if b==lastCh: D[a]=1
  for i in range(2):
    for edge in edgeList:
      a = edge[0]; b = edge[1]
      if D[b]==1: D[a]=1
  ok = 0
  for x in nuclList:
    if x!=lastCh and D[x]==0: return 0
  return 1

def eulerian(s):
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)
  #compute nucleotides appearing in s
  nuclList = []
  for x in ["A","C","G","U"]:
    if x in s: nuclList.append(x)
  #compute numInList[x] = number of dinucleotides beginning with x
  numInList = {}
  for x in nuclList:
    numInList[x]=0
    for y in nuclList:
      numInList[x] += dinuclCnt[x][y]
  #create dinucleotide shuffle L
  firstCh = s[0]  #start with first letter of s
  lastCh  = s[-1]
  edgeList = []
  for x in nuclList:
    if x!= lastCh: edgeList.append( [x,chooseEdge(x,dinuclCnt)] )
  ok = connectedToLast(edgeList,nuclList,lastCh)
  return ok,edgeList,nuclList,lastCh

def shuffleEdgeList(L):
  n = len(L); barrier = n
  for i in range(n-1):
    z = int(random.random() * barrier)
    tmp = L[z]
    L[z]= L[barrier-1]
    L[barrier-1] = tmp
    barrier -= 1
  return L

def dinuclShuffle(s):
  ok = 0
  while not ok:
    ok,edgeList,nuclList,lastCh = eulerian(s)
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)

  #remove last edges from each vertex list, shuffle, then add back
  #the removed edges at end of vertex lists.
  for [x,y] in edgeList: List[x].remove(y)
  for x in nuclList: shuffleEdgeList(List[x])
  for [x,y] in edgeList: List[x].append(y)

  #construct the eulerian path
  L = [s[0]]; prevCh = s[0]
  for i in range(len(s)-2):
    ch = List[prevCh][0]
    L.append( ch )
    del List[prevCh][0]
    prevCh = ch
  L.append(s[-1])
 # print(L)
  t = "".join(L)
  return t

#### Defining my functions #####
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def utf8len(s):
    return len(s.encode('utf-8'))

def write_fai (sequence, filename, name):
    w = open(filename, 'w')
    name = str(name)
    length = str(len(sequence))
    offset = str(utf8len(str(">"+name+"\n")))
    linebases = str(len(sequence))
    linewidth = str(len(sequence)+1)
    w.write("%s\t%s\t%s\t%s\t%s\n" % (name, length, offset, linebases, linewidth))

def write_fasta(sequence, outputfilename, name):
    w = open(outputfilename, 'w')
    fasta_sequence = sequence

    w.write(">"+name+"\n")
    w.write(str(fasta_sequence))

###### Function to calculate ZScore on list of MFEs #################
def pscore_function(energy_list, randomizations):
    below_native = 0
    total_count = len(energy_list)
    native_mfe = float(energy_list[0])
    #scrambled_mean_mfe = np.mean(energy_list[1:randomizations])
    for MFE in energy_list:
        if float(MFE) < float(native_mfe):
            below_native += 1

    pscore = float(float(below_native) / float(total_count))

    return pscore;

###### Function to calculate ZScore on list of MFEs #################
def zscore_function(energy_list, randomizations):
    mean = np.mean(energy_list)
    sd = np.std(energy_list)
    native_mfe = energy_list[0]
    scrambled_mean_mfe = np.mean(energy_list[1:randomizations])
    #scrambled_sd = np.std(energy_list[1:randomizations])
    if sd != 0:
        zscore = (native_mfe - scrambled_mean_mfe)/sd
    if sd == 0:
        zscore = float(00.00)
    return zscore;

def rna_folder(frag):
    fc = RNA.fold_compound(str(frag), md)
    (structure, MFE) = fc.mfe()
    return MFE;

def randomizer(frag):
    result = ''.join(random.sample(frag,len(frag)))
    return result;

###### Function to calculate MFEs using RNAfold #################
def energies(seq_list):
    energy_list = []

    energy_list = multiprocessing(rna_folder, [sequence for sequence in seq_list], 12)
    # for sequence in seq_list:
    #     #fc = RNA.fold_compound(str(sequence))
    #     (structure, MFE) = RNA.fold(str(sequence)) # calculate and define variables for mfe and structure
    #     energy_list.append(MFE) # adds the native fragment to list

    return energy_list;

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
        frag_seqs = multiprocessing(randomizer, [frag for i in range(randomizations)], 12)

        # for _ in range(int(randomizations)):
        #     result = ''.join(random.sample(frag,len(frag)))
        #     frag_seqs.append(result)
    else:
        print("Shuffle type not properly designated; please input \"di\" or \"mono\"")

    return frag_seqs;

def write_wig(metric_list, step, name, outputfilename):
    w = open(outputfilename, 'w')
    w.write("%s %s %s %s %s\n" % ("fixedStep", "chrom="+name, "start=0", "step="+str(step), "span="+str(step)))
    for metric in metric_list:
        if metric == "#DIV/0!":
            w.write("%s\n" % (metric))
        else:
            w.write("%f\n" % (metric))

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

# with open(myfasta, 'r') as forward_fasta:
with open (myfasta, 'r') as forward_fasta:
    if '>' in forward_fasta.readlines()[0]:
        print("Sequence in FASTA format...")
        with open (myfasta, 'r') as forward_fasta:
            for cur_record in SeqIO.parse(forward_fasta, "fasta"):
                ### get info about fasta files like name and sequence
                #read_name = cur_record.name #this reads fasta header (not reliable)
                seq = cur_record.seq
                #print(str(seq))
                print("Sequence Length: "+str(len(seq))+"nt")
                number_windows = int((len(seq)-int(window_size))/int(step_size)
                print("Approximately "+str(int(number_windows))+" windows will be generated.")
                print("Sequence being scanned...")
                if len(seq) > 20000 :
                    print(str(len(seq)))
                    raise SystemExit('Input sequence is longer than 20000 nt; in order to scan longer sequences consider using the stand alone programs (avaiable here: https://github.com/moss-lab/ScanFold)')
                length = len(seq)
                record_name = name

                # ### Add headers for bigwig files
                # """ Headers need to have the name and length of fasta sequence
                # Should be set up like XYZ_wig.addHeader("chromosomeName", length)
                # """
                # MFE_wig.addHeader([(str(record_name), int(length))])
                # zscore_wig.addHeader([(str(record_name), int(length))])
                # pvalue_wig.addHeader([(str(record_name), int(length))])
                # ED_wig.addHeader([(str(record_name), int(length))])

                ### Write the header for output:
                w.write("i\tj\tTemperature\tNative_dG\tZ-score\tP-score\tEnsembleDiversity\tSequence\tStructure\tCentroid\t"+record_name+"\n")
                i = 0
            ##### Main routine using defined functions: ##########################################
            ### Figure out length progress bar
                # pbar = ProgressBar(widgets=widgets, max_value=int(length))
                # print(int(length/step_size))
                # pbar.start()
                while i == 0 or i <= (length - window_size):
                    #print(i)
                    # pbar.update(i)
                    start_nucleotide = i + 1 # This will just define the start nucleotide coordinate value
                    frag = seq[i:i+int(window_size)] # This breaks up sequence into fragments
                    #print(frag)
                    #print(str(len(frag)))
                    start_nucleotide = i + 1
                    end_nucleotide = i + window_size
                    if -1 == 0:
                        print("Magic")
                    # if 'N' in frag:
                    #     w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str(frag)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\n")
                    #     i += step_size #this ensures that the next iteration increases by "step size" length
                    else:
                        #print(start_nucleotide)
                        #print(end_nucleotide)
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
                            fc = RNA.fold_compound(str(frag), md) #creates "Fold Compound" object
                            fc.pf() # performs partition function calculations
                            #frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                            #fmfe = fc.pbacktrack()
                            #print(str(fmfe))
                            seqlist = [] # creates the list we will be filling with sequence fragments
                            seqlist.append(frag) # adds the native fragment to list
                            scrambled_sequences = scramble(frag, randomizations, type)
                            seqlist.extend(scrambled_sequences)
                            energy_list = energies(seqlist)
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

    else:
        with open (myfasta, 'r') as forward_fasta:
            print("Sequence not in FASTA format, attempting to convert. You can also resubmit sequence with a proper FASTA header (i.e. >SequenceName)")
            raw_sequence = forward_fasta.read()
            seq = re.sub("\n", "", raw_sequence.strip())
            #print(seq)
            print("Sequence Length: "+str(len(seq))+"nt")
            number_windows = int(len(seq))/int(step_size)
            print("Approximately "+str(int(number_windows))+" windows will be generated.")
            print("Sequence being scanned...")
            if len(seq) > 20000 :
                print(str(len(seq)))
                raise SystemExit('Input sequence is longer than 20000 nt; in order to scan longer sequences consider using the stand alone programs (avaiable here: https://github.com/moss-lab/ScanFold)')
            length = len(seq)
            record_name = name
            i = 0
        ##### Main routine using defined functions: ##########################################
        ### Figure out length progress bar
            # pbar = ProgressBar(widgets=widgets, max_value=int(length))
            # print(int(length/step_size))
            # pbar.start()
            while i == 0 or i <= (length - window_size):
                #print(i)
                # pbar.update(i)
                start_nucleotide = i + 1 # This will just define the start nucleotide coordinate value
                frag = seq[i:i+int(window_size)] # This breaks up sequence into fragments
                #print(frag)
                #print(str(len(frag)))
                start_nucleotide = i + 1
                end_nucleotide = i + window_size
                if -1 == 0:
                    print("Magic")
                # if 'N' in frag:
                #     w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str(frag)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\n")
                #     i += step_size #this ensures that the next iteration increases by "step size" length
                else:
                    #print(start_nucleotide)
                    #print(end_nucleotide)
                    # frag = frag.transcribe()
                    if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                        MFE = int(0.0)
                        zscore = "00.00"
                        ED = int(0.0)
                        pscore = int(0.0)
                        structure = "........................................................................................................................"
                        centroid = "........................................................................................................................"
                    else:
                        fc = RNA.fold_compound(str(frag), md) #creates "Fold Compound" object
                        fc.pf() # performs partition function calculations
                        #frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
                        (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                        MFE = round(MFE, 2)
                        MFE_total.append(MFE)
                        (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                        ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                        ED_total.append(ED)            #print(structure)
                        #fmfe = fc.pbacktrack()
                        #print(str(fmfe))
                        seqlist = [] # creates the list we will be filling with sequence fragments
                        seqlist.append(frag) # adds the native fragment to list
                        scrambled_sequences = scramble(frag, randomizations, type)
                        seqlist.extend(scrambled_sequences)
                        energy_list = energies(seqlist)
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

write_wig(MFE_list, step_size, name, mfe_wig_file_path)
write_wig(zscore_list, step_size, name, zscore_wig_file_path)
write_wig(pscore_list, step_size, name, pvalue_wig_file_path)
write_wig(ED_list, step_size, name, ed_wig_file_path)

write_fasta(seq, fasta_file_path, name)
write_fai(seq, fasta_index, name)


if split == "off":
    url = str(callbackurl+"/"+str(nodeid)+"/0")
    response = requests.get(url, verify=False)
    # print(url)

if split == "on":
    print("Completed ScanFold-Scan, intial results shown below.\n ScanFold-Fold is running now")

print("Mean MFE = "+str(mean_MFE))
print("Mean Z-score = "+str(mean_zscore))
print("Mean P-value = "+str(mean_pscore))
print("Mean Ensemble Diversity = "+str(mean_ED))
print("ScanFold-Scan complete, find output files below")
