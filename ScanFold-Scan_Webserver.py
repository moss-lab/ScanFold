#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3.6

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
$ python3.6 ScanFold-Scan.py 1-input 2-stepsize 3-window size 4-randomizations
    5-temperature 6-shuffle type

    1.
    2.
    3.
    4.
    5.
    6. should be "mono" or "di"

"""

import sys
import argparse
import string
#import sequence
#import csv
#import argparse
import re
import numpy as np
sys.path.append('/home/randrews/ViennaRNA/lib/python/site-packages/')
sys.path.append('/usr/local/lib/python3.6/site-packages')
import RNA
import random
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from Bio import SeqIO


#### Defining global variables ###############

myfasta = sys.argv[1] #input filename
step_size = int(sys.argv[2])
window_size = int(sys.argv[3])
randomizations = int(sys.argv[4])
temperature = int(sys.argv[5])
type = str(sys.argv[6])
w = open(myfasta+".forward.win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)+".txt", 'w')
s = open("result_summary.forward."+myfasta+".win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)+".txt", 'w')
s.write("ReadName\tLength\tMeanMFE\tMeanZ\tMeanP\tMeanED\n")
# r = open(myfasta+".reverse.win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)+".txt", 'w')
# rs = open("result_summary.reverse."+myfasta+".win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".txt", 'w')
# rs.write("ReadName\tLength\tMeanMFE\tMeanZ\tMeanP\tMeanED\n")

md = RNA.md()
md.temperature = int(temperature)

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
        zscore = "#DIV/0!"
    return zscore;

def rna_folder(frag):
    (structure, MFE) = RNA.fold(str(frag))
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

##################### Main Script #########################################
with open(myfasta, 'r') as forward_fasta:

    for cur_record in SeqIO.parse(forward_fasta, "fasta") :

            read_name = cur_record.name

            #### this will change based on input fasta file header format #########
            #print(read.name)
            #fasta_header = read_name.split('|')
            #print(fasta_header)
            #gene_id = fasta_header[0]
            #transcript_id = fasta_header[1]
            #chromosome = "chr"+fasta_header[2]
            #gene_start = fasta_header[3]
            #gene_end = fasta_header[4]
            #strand = fasta_header[5]
        ##### Establish empty lists to capture calculated metrics per window ######
            zscore_total = []
            numerical_z = []
            pscore_total = []
            numerical_p = []
            MFE_total = []
            ED_total = []

            #gff3file = open(read.name+'.gff3', 'w')
            #pscore_wig = open(read.name+'.pscore.wig', 'w')
            #zscore_wig = open(read.name+".zscore.wig", 'w')
            #ED_wig = open(read.name+".ED.wig", 'w')
            #MFE_wig = open(read.name+".MFE.wig", 'w')
            #print(read.name, read.sequence)
            length = len(cur_record.seq)
            seq = cur_record.seq
            #print(length)
            w.write("i\tj\tTemperature\tNative_dG\tZ-score\tP-score\tEnsembleDiversity\tSequence\tStructure\tCentroid\t"+read_name+"\n")
            i = 0


    ##### Main routine using defined functions: ##########################################

            while i == 0 or i <= (length - window_size):
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
                    frag = frag.transcribe()
                    fc = RNA.fold_compound(str(frag)) #creates "Fold Compound" object
                    fc.pf() # performs partition function calculations
                    frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
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
                    print(energy_list)
                    zscore = round(zscore_function(energy_list, randomizations), 2)
                    zscore_total.append(zscore)

                    #print(zscore)
                    pscore = round(pscore_function(energy_list, randomizations), 2)
                    #print(pscore)
                    pscore_total.append(pscore)

                    print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
                    w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(zscore)+"\t"+str(pscore)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
                    #gff3file.write()
                    #pscore_wig.write()
                    #zscore_wig.write()
                    #ED_wig.write()
                    #MFE_wig.write()

                    i += step_size #this ensures that the next iteration increases by "step size" length

            #print(len(zscore_total))
            for z in zscore_total:
                try:
                    numerical_z.append(float(z))
                except ValueError:
                    continue
            #print(len(numerical_z))

            #print(len(pscore_total))
            for p in pscore_total:
                try:
                    numerical_p.append(float(p))
                except ValueError:
                    continue
            window_count = len(zscore_total)
            #print(len(numerical_p))
            #print(window_count)
            #print(type(window_count))
            #print(step_size)
            #print(type(step_size))
            #print(length)
            #print(type(length))
            #coverage = round((float(window_count)*float(step_size))/float(length), 2)
            #print(coverage)
            #print(len(MFE_total))
            #print(len(ED_total))
            mean_pscore = round(np.mean(numerical_p), 2)
            mean_zscore = round(np.mean(numerical_z), 2)
            mean_MFE = round(np.mean(MFE_total), 2)
            mean_ED = round(np.mean(ED_total), 2)
            #w.write("---\t---\t---\t---\t---\t---\t---\t---\tSummary:\tLength\tMeanMFE\tMeanZ\tMeanPscore\tMeanED\n---\t---\t---\t---\t---\t---\t---\t---\t---\t"+str(length)+"\t"+str(mean_MFE)+"\t"+str(mean_zscore)+"\t"+str(mean_pscore)+"\t"+str(mean_ED)+"\n\n")
            s.write(str(read_name)+"\t"+str(length)+"\t"+str(mean_MFE)+"\t"+str(mean_zscore)+"\t"+str(mean_pscore)+"\t"+str(mean_ED)+"\n")
