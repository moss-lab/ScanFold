#!/usr/local/bin/python3.6
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

ScanFold-Fold
Contact: Ryan Andrews - randrews@iastate.edu

This program will take the output from the ScanFold-Scan program (which can be
found at https://github.com/moss-lab/ScanFold/ScanFold-Scan.pl) and condense
the entirity of the scanning window results into a single structure. Only the
most unusually stable base pairs will be reported.

Usage:
ScanFold-Fold_IGV.py -i out.tsv --structure_extract_file ./extracted_structures.txt --out1 ./nofilter.ct --out2 ./scanfold-1.ct --out3 ./scanfold-2.ct --out4 ./log_file --out5 ./final_partner_log --out6 ./bp_track --out7 ./fasta.fa --dbn_file_path ./scanfold.dbn  --dbn_file_path1 ./scanfold1.dbn  --dbn_file_path2 ./scanfold2.dbn --dbn_file_path3 ./scanfold3.dbn --dbn_file_path4 ./scanfold4.dbn --fasta_index ./fasta.fa.fai --final_partners_wig ./fp.wig --nodeid "/scholar" --callbackurl "https://www.google.com"

    1. Name of output file from ScanFold-Scan

    2. Cutoff value which will be used in addition to the default values of -2,
      -1 and 10.

    3. The standard output is formatted as a log file which contains 1. a list
      of all nucleotides and their predicted base pairing partners; 2. a list of
      the most favorable base pairs.

"""

import math
import itertools
import operator
from collections import Counter, defaultdict
import sys
import re
import numpy as np
import os
# sys.path.append('/Users/ryanandrews/Desktop/programs/RNAstructure/exe')
# import RNAstructure
sys.path.append('/usr/local/lib/python3.6/site-packages')
import time
import argparse
from itertools import repeat
from functools import partial
#for mono z-score
import random
from io import StringIO
import tempfile
import subprocess

from ScanFoldSharedIGV import *


class NucPair:
    #Class to define a base pair
    def __init__(self, inucleotide, icoordinate, jnucleotide, jcoordinate, zscore, mfe, ed):
        self.inucleotide = inucleotide
        self.icoordinate = icoordinate
        self.jnucleotide = jnucleotide
        self.jcoordinate = jcoordinate
        self.zscore = zscore
        self.mfe = mfe
        self.ed = ed

class NucStructure:
    #Class to define a base pair
    def __init__(self, bond_order, coordinate, nucleotide, structure):
        self.bond_order = bond_order
        self.coordinate = coordinate
        self.nucleotide = nucleotide
        self.structure = structure

class NucStructureCount:
    #Class to define a base pair
    def __init__(self, structure_count, coordinate, nucleotide, structure):
        self.structure_count = structure_count
        self.coordinate = coordinate
        self.nucleotide = nucleotide
        self.structure = structure

class ExtractedStructure:
    def __init__(self, structure_count, sequence, structure, i, j):
        self.structure_count = structure_count
        self.sequence = sequence
        self.structure = structure
        self.i = i
        self.j = j

class NucZscore:
    #Nucleotide class; defines a nucleotide with a coordinate and a A,T,G,C,U
    def __init__(self, nucleotide, coordinate):
        self.nucleotide = nucleotide
        self.coordinate = coordinate

    def add_zscore(self, zscore):
        self.zscores.append(zscore)

    def add_pair(self, pair):
        self.pair.append(pair)

def NucleotideDictionary (lines):
    """
    Function to generate nucleotide dictionary where each key is the i
    coordinate of the nucleotide of the input sequence, and each value is a
    NucZscore class object (which contains the coordinate and nucleotide
    informations)
    """
    nuc_dict = {}
    for row in lines:
        if not row.strip():
            continue
        else:
            i = 1
            data = row.split('\t')
            icoordinate = data[0]
            sequence = str(data[7])
            if "-" in sequence:
                sequence = sequence.replace("-", "N")
            if str(data[8]).find("A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u") > 0:
                #print("8"+str(data[8]))
                sequence = transcribe(str(data[8]))
                sequence = sequence.replace("-", "N")
            if str(data[7]).find("A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u") > 0:
                #print("7")
                sequence = transcribe(str(data[7]))
                sequence = sequence.replace("-", "N")
            # else:
            #     raise("Could not find sequence for window")
            #
            # except:
            #     data = row.split(',')
            #     strand = int(data[11])
            #     #print(strand)
            #     icoordinate = data[0]
            #     if ("A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u") in str(data[8]):
            #         sequence_raw = transcribe(str(data[8]))
            #     elif ("A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u") in str(data[7]):
            #         sequence_raw = transcribe(str(data[7]))
            #     else:
            #         raise("Could not find sequence for window")
            #
            #
            #     #print(sequence_raw)
            #     if strand == -1:
            #         #print("NegStrand")
            #         sequence = sequence_raw[::-1]
            #         #print(sequence)
            #     else:
            #         #print("PosStrand")
            #         sequence = sequence_raw
            #
            for nuc in sequence:
                #print(nuc)
                x = NucZscore(nuc,(int(icoordinate)+int(i)-1))
                nuc_dict[x.coordinate] = x
                i += 1

    return nuc_dict;

def competing_pairs(bp_dict, coordinate):
    #Function to determine other i-nuc which compete for the same j-nuc
    comp_pairs = {}
    i = 0
    for k, v in bp_dict.items():
        if ((int(v.jcoordinate) == int(coordinate)) or
            (int(v.icoordinate) == int(coordinate))):
            x = NucPair(v.inucleotide, v.icoordinate, v.jnucleotide,
                        v.jcoordinate, v.zscore, v.mfe, v.ed)
            comp_pairs[i] = x
            i += 1
        else:
            continue

    return comp_pairs;

def best_basepair(bp_dict, nucleotide, coordinate, type):
    #Function to define best i-j pair for i-nucleotide
    zscore_dict = {}
    pair_dict = {}
    partner_key = 0
    for k, pair in sorted(bp_dict.items()):
        if int(pair.icoordinate) < int(pair.jcoordinate):
            #print("148")
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)
            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                z = pair_dict[partner_key]
                z.append(x)

            except:
                zscore_dict[partner_key] = []
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                pair_dict[partner_key] = []
                z = pair_dict[partner_key]
                z.append(x)

            sum_z = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = np.sum(v1)
                test = sum_z[k1] = np.sum(v1)

            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = np.mean(v1)
                test = mean_z[k1] = np.mean(v1)

            partner_key += 1

        elif int(pair.icoordinate) > int(pair.jcoordinate):
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)

            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                z = pair_dict[partner_key]
                z.append(x)

            except:
                zscore_dict[partner_key] = []
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                pair_dict[partner_key] = []
                z = pair_dict[partner_key]
                z.append(x)

            sum_z = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = np.sum(v1)
                test = sum_z[k1] = np.sum(v1)

            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = np.mean(v1)
                test = mean_z[k1] = np.mean(v1)

            partner_key += 1

        elif int(pair.icoordinate) == int(pair.jcoordinate):
            #print("210")
            x = NucPair(pair.inucleotide, pair.icoordinate, pair.jnucleotide,
                        pair.jcoordinate, pair.zscore, pair.mfe, pair.ed)
            try:
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                z = pair_dict[partner_key]
                z.append(x)

            except:
                zscore_dict[partner_key] = []
                y = zscore_dict[partner_key]
                y.append(pair.zscore)
                pair_dict[partner_key] = []
                z = pair_dict[partner_key]
                z.append(x)

            sum_z = {}
            for k1, v1 in zscore_dict.items():
                sum_z[k1] = np.sum(v1)
                test = sum_z[k1] = np.sum(v1)

            mean_z = {}
            for k1, v1 in zscore_dict.items():
                mean_z[k1] = np.mean(v1)
                test = mean_z[k1] = np.mean(v1)

            partner_key += 1

        else:
            print("FAIL")
            best_bp = NucPair(pair.inucleotide, pair.icoordinate,
                              pair.jnucleotide, pair.jcoordinate, pair.zscore)

            partner_key += 1

        if type == 'sum':
            best_bp_key = min(sum_z, key = sum_z.get)
        if type == 'mean':
            best_bp_key = min(mean_z, key = mean_z.get)

    try:
        v = pair_dict[best_bp_key]
        best_bp = v[0]
    except:
        print("ERROR")
        print(k)

    return best_bp;

def write_ct(base_pair_dictionary, filename, filter, strand):
    #Function to write connectivity table files from a list of best i-j pairs
    w = open(filename, 'w')
    w.write((str(len(base_pair_dictionary))+"\t"+name+"_Filter="+str(filter)+"\n"))
    strand = 1
    if strand == 1:
        for k, v in base_pair_dictionary.items():
            #print(start_coordinate)
            #print(v.icoordinate)
            icoordinate = str(int(v.icoordinate)-int(int(start_coordinate)-1))
            #print(icoordinate)
            jcoordinate = str(int(v.jcoordinate)-int(int(start_coordinate)-1))
            #print(jcoordinate)
            key_coordinate = str(int(k)-int(start_coordinate)+1)
            #print(key_coordinate)
            if float(v.zscore) < filter:
                if ((int(icoordinate) < int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif ((int(icoordinate) > int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif (int(icoordinate) < int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif (int(icoordinate) > int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif int(icoordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                #
                # elif (int(key_coordinate) != icoordinate) and (int(key_coordinate) != int(jcoordinate)):
                #     continue
                #     #w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("Error at", int(key_coordinate), v.inucleotide, icoordinate, v.jnucleotide, int(jcoordinate), v.zscore)
            else:
                if int(key_coordinate) == int(icoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                elif int(key_coordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    raise ValueError("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key_coordinateey_coordinateey)")
                continue

    if strand == -1:
        for k, v in sorted(base_pair_dictionary.items(), key=lambda x:x[0], reverse = True):
            # print(start_coordinate)
            # print(end_coordinate)
            # print("i="+str(v.icoordinate))
            # print("j="+str(v.jcoordinate))
            # print("k="+str(k))
            icoordinate = str(int(end_coordinate)+1-(int(int(v.icoordinate))))
            # print("i_after"+str(icoordinate))
            jcoordinate = str(int(end_coordinate)+1-(int(int(v.jcoordinate))))
            # print("j_after="+str(jcoordinate))
            key_coordinate = str(int(end_coordinate)-int(k)+1)
            # print("key="+str(key_coordinate))
            if float(v.zscore) < filter:
                if ((int(icoordinate) < int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif ((int(icoordinate) > int(jcoordinate)) and (int(icoordinate) == int(key_coordinate))): #test to see if reverse bp.
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(jcoordinate), int(key_coordinate)))

                elif (int(icoordinate) < int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif (int(icoordinate) > int(jcoordinate)) and (int(key_coordinate) == int(jcoordinate)):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, int(icoordinate), int(key_coordinate)))

                elif int(icoordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                #
                # elif (int(key_coordinate) != icoordinate) and (int(key_coordinate) != int(jcoordinate)):
                #     continue
                #     #w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    print("Error at", int(key_coordinate), v.inucleotide, icoordinate, v.jnucleotide, int(jcoordinate), v.zscore)
            else:
                if int(key_coordinate) == int(icoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.inucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                elif int(key_coordinate) == int(jcoordinate):
                    w.write("%d %s %d %d %d %d\n" % (int(key_coordinate), v.jnucleotide, int(key_coordinate)-1, int(key_coordinate)+1, 0, int(key_coordinate)))
                else:
                    raise ValueError("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key_coordinateey_coordinateey)")
                continue

def write_dp(base_pair_dictionary, filename, filter):
    #this function will create a dp file for IGV
    w = open(filename, 'w')
    for k, v in base_pair_dictionary.items():
        if float(v.zscore) < filter:
            probability = (v.zscore/minz)
            if int(v.icoordinate) < int(v.jcoordinate):
                #w.write("%d\t%d\t%f\n" % (k, int(v.jcoordinate), float(-(math.log10(probability)))))
                w.write("%d\t%d\t%f\n" % (v.icoordinate, int(v.jcoordinate), float(((-1/minz)*(v.zscore)))/minz))
            elif int(v.icoordinate) > int(v.jcoordinate):
                w.write("%d\t%d\t%f\n" % (int(v.icoordinate), int(v.jcoordinate), float(((-1/minz)*(v.zscore)))/minz))
            elif int(v.icoordinate) == int(v.jcoordinate):
                w.write("%d\t%d\t%f\n" % (k, int(v.jcoordinate), float(((-1/minz)*(v.zscore)))/minz))
            else:
                print("Error at:", k)

def transcribe(seq):
    #Function to covert T nucleotides to U nucleotides
    for ch in seq:
        rna_seq = seq.replace('T', 'U')
        return(rna_seq)

def flip_structure(structure):
    #Function to reverse structure in a given window, for negative strand genes
    flip = {'(':')', ')':'(', '.':'.'}
    return ''.join([flip[pair] for pair in structure[::-1]])

def write_fasta(nucleotide_dictionary, outputfilename, name):
    w = open(outputfilename, 'w')
    fasta_sequence = str()
    for k, v in nucleotide_dictionary.items():
        nucleotide = v.nucleotide
        fasta_sequence += nucleotide

    w.write(">"+name+"\n")
    w.write(str(fasta_sequence))

def write_wig_dict(nucleotide_dictionary, outputfilename, name, strand):

    w = open(outputfilename, 'w')
    #write wig file header
    w.write("%s %s %s %s %s\n" % ("fixedStep", "chrom="+name, "start="+str(start_coordinate), "step=1", "span=1"))

    #write values of zscores
    if strand == "reverse":
        for k, v, in sorted(nucleotide_dictionary.items(), reverse=True):
            w.write("%f\n" % (v.zscore))
    else:
        for k, v in nucleotide_dictionary.items():
            w.write("%f\n" % (v.zscore))

def write_bp(base_pair_dictionary, filename, start_coordinate, bpstrand, length):
    with open(filename, 'w') as w:
            #set color for bp file (igv format)
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 55, 129, 255, str("Less than -2 "+str(minz))))
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 89, 222, 111, "-1 to -2"))
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 236, 236, 136, "0 to -1"))
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 199, 199, 199, "0"))
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 228, 228, 228, "0 to 1"))
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 243, 243, 243, "1 to 2"))
        w.write("%s\t%d\t%d\t%d\t%s\n" % (str("color:"), 247, 247, 247, str("Greater than 2")))
        reverse_list = []
        i = 0
        for k, v in base_pair_dictionary.items():
            #choose color
            if float(v.zscore) < float(-2):
                score = str(0)
                #print(k, v.zscore, score)

            elif (float(v.zscore) < int(-1)) and (float(v.zscore) >= -2):
                score = str(1)
                #print(k, v.zscore, score)

            elif (float(v.zscore) < int(0)) and (float(v.zscore) >= -1):
                score = str(2)
                #print(k, v.zscore, score)

            elif float(v.zscore) == 0 :
                score = str(3)
                #print(k, v.zscore, score)

            elif 0 < float(v.zscore) <= 1:
                score = str(4)
                #print(k, v.zscore, score)

            elif 1 < float(v.zscore) <= 2:
                score = str(5)
                #print(k, v.zscore, score)

            elif float(v.zscore) > 2:
                score = str(6)
                #print(k, v.zscore, score)

            else:
                print(k, v.zscore, score)


            score = str(score)

            # ensure coordinates to start at 1 to match with converted fasta file
            # if bpstrand == "forward":
                # sc = int(int(start_coordinate)-1)
            #print(length)

            if int(v.icoordinate) < int(v.jcoordinate):
                #w.write("%d\t%d\t%f\n" % (k, int(v.jcoordinate), float(-(math.log10(probability)))))
                # print(name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score)
                # print(name, str(int(v.icoordinate)-sc), str(int(v.icoordinate)-sc), str(int(v.jcoordinate)-sc), str(int(v.jcoordinate)-sc), score)
                if bpstrand == "forward":
                    w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score))
                reverse_list.append(list([name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score]))

            elif int(v.icoordinate) > int(v.jcoordinate):
                # print(name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score)
                # print(name, str(int(v.icoordinate)-sc), str(int(v.icoordinate)-sc), str(int(v.jcoordinate)-sc), str(int(v.jcoordinate)-sc), score)
                if bpstrand == "forward":
                    w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score))
                reverse_list.append(list([name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score]))

            elif int(v.icoordinate) == int(v.jcoordinate):
                # print(name, int(v.icoordinate), int(v.icoordinate), int(v.jcoordinate), int(v.jcoordinate), score)
                # print(name, str(int(v.icoordinate)-sc), str(int(v.icoordinate)-sc), str(int(v.jcoordinate)-sc), str(int(v.jcoordinate)-sc), score)
                if bpstrand == "forward":
                    w.write("%s\t%d\t%d\t%d\t%d\t%s\n" % (name, k, k, int(v.jcoordinate), int(v.jcoordinate), score))
                reverse_list.append(list([name, k, k, int(v.jcoordinate), int(v.jcoordinate), score]))
            else:
                continue

        if bpstrand == "reverse":
            #length should be in final row, column 2
            final_row = reverse_list[-1]
            first_row = reverse_list[0]
            ### correcting an off by one error here
            start_coordinate = first_row[1]-1
            length = final_row[1]
            #print(start_coordinate, length)
            #Print header
            #establish a reverse order of lines, modify, then print
            for row in reversed(reverse_list):
                rev0 = str(row[0])
                rev1 = str((int(length) - int(row[1])+1)+int(start_coordinate))
                rev2 = str((int(length) - int(row[2])+1)+int(start_coordinate))
                rev3 = str((int(length) - int(row[3])+1)+int(start_coordinate))
                rev4 = str((int(length) - int(row[4])+1)+int(start_coordinate))
                rev5 = str(row[5])

                # rev0 = str(row[0])
                # rev1 = str((int(length) - int(row[1])+1))
                # rev2 = str((int(length) - int(row[2])+1))
                # rev3 = str((int(length) - int(row[3])+1))
                # rev4 = str((int(length) - int(row[4])+1))
                # rev5 = str(row[5])
                #print(rev0+"\t"+rev1+"\t"+rev2+"\t"+rev3+"\t"+rev4+"\t"+rev5)
                w.write(rev0+"\t"+rev1+"\t"+rev2+"\t"+rev3+"\t"+rev4+"\t"+rev5+"\n")

def get_step_win(lines):
    row1 = lines[0].split()
    row2 = lines[1].split()
    step_size = int(row2[0])-int(row1[0])
    window_size = int(row1[1])-int(row1[0])

    return step_size, window_size

######Function to create X number of scrambled RNAs in list #################
def scramble(text, randomizations, type):
    frag = str(text)
    frag_seqs = []
    if shuffle == "di":
        for _ in range(randomizations):
            result = dinuclShuffle(frag)
            frag_seqs.append(result)
    elif shuffle == "mono":
        for _ in range(int(randomizations)):
            result = ''.join(random.sample(frag,len(frag)))
            frag_seqs.append(result)
    else:
        print("Shuffle type not properly designated; please input \"di\" or \"mono\"")

    return frag_seqs;


if __name__ == "__main__":

    start_time = time.time()


    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input filename')
    parser.add_argument('-f', type=int, default=-2,
                        help='filter value')
    parser.add_argument('-c', type=int, default=1,
                        help='Competition')
    parser.add_argument('-d', type=str, default = "forward",
                        help='strand of genome (forward or reverse; default forward)')
    parser.add_argument('-l', type=str,
                        help='length of ScanFold-Scan input sequence')

    ### Required for spinoff ###
    parser.add_argument('--out1', type=str,
                        help='out1 path')
    parser.add_argument('--out2', type=str,
                        help='out3 path')
    parser.add_argument('--out3', type=str,
                        help='out3 path')
    parser.add_argument('--out4', type=str,
                        help='log_file path')
    parser.add_argument('--out5', type=str,
                        help='final_partner_file path')
    parser.add_argument('--out6', type=str,
                        help='bp_track_file path')
    parser.add_argument('--out7', type=str,
                        help='fasta_file path')
    parser.add_argument('--dbn_file_path', type=str,
                        help='dbn_file_path')
    parser.add_argument('--dbn_file_path1', type=str,
                        help='dbn_file_path1')
    parser.add_argument('--dbn_file_path2', type=str,
                        help='dbn_file_path2')
    parser.add_argument('--dbn_file_path3', type=str,
                        help='dbn_file_path3')
    parser.add_argument('--dbn_file_path4', type=str,
                        help='dbn_file_path4')
    parser.add_argument('--structure_extract_file', type=str,
                        help='structure_extract_file path')



    parser.add_argument('--global_refold', action='store_true',
                        help='global refold oprion')


    parser.add_argument('--fasta_index', type=str,
                        help='fasta index file path')
    parser.add_argument('--name', type=str, default = "UserInput",
                        help='name of data being analyzied')
    parser.add_argument('--final_partners_wig', type=str,
                        help='final partners wig file path')
    parser.add_argument('-t', '--temp', type=int, default=37,
                        help='Folding temperature')

    args = parser.parse_args()

    temperature = args.temp

    strand = args.d
    #strand = "reverse"
    filename = args.input
    filter = int(args.f)
    competition = int(args.c)
    out1 = args.out1
    out2 = args.out2
    out3 = args.out3
    out4 = args.out4
    out5 = args.out5
    out6 = args.out6
    out7 = args.out7
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


    try:
        options = str(sys.argv[3])
    except:
        options = None

    #print(options, filter)
    output_data = re.split('\.', str(filename))
    output = str(str(filename)+".ScanFold.")


    # log_total = open(str(filename)+".ScanFold.log.txt", 'w')
    # log_win = open(str(filename)+".ScanFold.final_partners.txt", 'w')
    log_total = open(out4, 'w')
    log_win = open(out5, 'w')

    #Begin parsing file - Main Loop
    with open(filename, 'r') as f:
        #Initialize bp dictionary and z-score lists
        z_score_list = []
        bp_dict = {}

        #Read all lines from ScanFold-Scan file (exept header)
        lines = f.readlines()[0:]
        #print(lines)
        #Generate nucleotide dictionary to assign each nucleotide in sequence a key
        nuc_dict = NucleotideDictionary(lines)
        step_size, window_size = get_step_win(lines)
        print("Sequence length: "+str(len(nuc_dict))+"nt")
        # if len(nuc_dict) > 1000:
        #     raise SystemExit('Input sequence is longer than 1000 nt; in order to scan longer sequences consider using the stand alone programs (avaiable here: https://github.com/moss-lab/ScanFold)')

        #Determine start and end coordinate values
        start_coordinate = str(list(nuc_dict.keys())[0])
        #print(start_coordinate)
        end_coordinate = str(list(nuc_dict.keys())[-1])
        #print(end_coordinate)
        ### get step and window size
        step_size, window_size = get_step_win(lines)
        #print(step_size, window_size)


        #Iterate through input file, read each rows metrics, sequence, etc.
        print("Reading sequence and structures...")

        for row in lines:

            #Ignore blank lines
            if not row.strip():
                continue

            #Main loop to find all i-j pairs per i-nucleotide
            else:
                #Assign metrics to variables
                try:
                    data = row.split('\t')

                    icoordinate = data[0]
                    jcoordinate = data[1]
                    temp = data[2]
                    mfe = float(data[3])
                    zscore = float(data[4])
                    pvalue = data[5]
                    ed = float(data[6])
                    sequence_raw = transcribe(str(data[7]))
                    sequence_raw = sequence_raw.replace("-", "N")
                    #print(sequence_raw)
                    structure_raw = str(data[8])
                    #print(structure_raw)
                    if ("A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u") in str(data[8]):
                        #print("8"+str(data[8]))
                        fmfe = float(data[7])
                        sequence_raw = transcribe(str(data[8]))
                        sequence_raw = sequence_raw.replace("-", "N")
                        structure_raw = str(data[9])

                    elif ("A" or "G" or "C" or "T" or "U" or "a" or "g" or "c" or "t" or "u") in str(data[7]):
                        #print("7")
                        sequence_raw = transcribe(str(data[7]))
                        sequence_raw = sequence_raw.replace("-", "N")
                        structure_raw = str(data[8])

                    #strand = 1
                    #print("Tab "+icoordinate)
                except:
                    data = row.split(',')
                    icoordinate = data[0]
                    jcoordinate = data[1]
                    temp = data[2]
                    mfe = float(data[3])
                    zscore = float(data[4])
                    pvalue = data[5]
                    ed = float(data[6])
                    fmfe = data[7]
                    sequence_raw = transcribe(str(data[8]))
                    sequence_raw = sequence_raw.replace("-", "N")
                    structure_raw = str(data[9])
                    strand = int(data[11])
                    if strand == -1:
                        #print(icoordinate)
                        sequence_forward = sequence_raw
                        sequence_reverse = sequence_forward[::-1]
                        structure_forward = structure_raw
                        #print(structure_forward)
                        structure_reverse = flip_structure(structure_forward)
                        #print(structure_reverse)
                        # print(sequence_raw)
                        # print(sequence_reverse)
                        structure_raw = structure_reverse
                        sequence_raw = sequence_reverse
                        icoordinate = data[0]
                        jcoordinate = data[1]
                    if strand == int(1):
                        icoordinate = data[0]
                        jcoordinate = data[1]

                    #print("Comma "+icoordinate)

                #Convert sequence and structures into lists
                sequence_raw = sequence_raw.replace("-", "N")
                sequence = list(sequence_raw)
                structure = list(structure_raw)

                #Define window coordinates as string
                #window = str(str(icoordinate)+"-"+str(jcoordinate))

                #Determine length of window
                length = len(sequence)

                #Append window z-score to list (to calculate overall z-score)
                z_score_list.append(zscore)

                #Iterate through dot bracket structure to determine locations of bps
                i = 0
                while i < length:
                    #Unpaired nucleotide
                    if structure[i] == '.':
                        nucleotide = sequence[i]
                        coordinate = (i + int(icoordinate))
                        x = NucPair(nucleotide, coordinate, nucleotide, coordinate,
                                    zscore, mfe, ed)
                        try:
                            y = bp_dict[coordinate]
                            y.append(x)
                        except:
                            bp_dict[coordinate] = []
                            y = bp_dict[coordinate]
                            y.append(x)
                        i += 1
                    #Paired nucleotide
                    else:
                        i += 1

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
                        print("Error")

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

                    lbp_coord = int(int(lbp)+int(icoordinate)-1)
                    rbp_coord = int(int(rbp)+int(icoordinate)-1)
                    x = NucPair(lb, lbp_coord, rb, rbp_coord, zscore, mfe, ed)
                    z = NucPair(rb, rbp_coord, lb, lbp_coord, zscore, mfe, ed)

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
                        w = bp_dict[rbp_coord]
                        w.append(z)
                    #If i-nuc not defined, define it
                    except:
                        bp_dict[rbp_coord] = []
                        w = bp_dict[rbp_coord]
                        w.append(z)
                    l += 2

            #Define OVERALL values of metrics
            meanz = float(np.mean(z_score_list))
            sdz = float(np.std(z_score_list))
            minz = min(z_score_list)
            stdz = float(np.std(z_score_list))

            one_sig_below = float(meanz-stdz)
            two_sig_below = float( meanz - ( 2 * stdz) )

    #Initiate global dictionaries to store best base pairs
    best_bps = {}
    best_sum_bps = {}
    best_sum_bps_means = {}
    best_total_window_mean_bps = {}

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
            if "-" in pair.jnucleotide:
                #print("Gap found. Converting to N.")
                partner_key = str("N")+"-"+str(pair.jcoordinate)
            else:
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
            sum_z[k1] = np.sum(v1)
            test = sum_z[k1] = np.sum(v1)
            sum_z_lengths[k1] = len(sum_z)

        #Caclulate and store mean of z-score per i-j pair
        mean_z = {}
        for k1, v1 in zscore_dict.items():
            mean_z[k1] = np.mean(v1)
            test = mean_z[k1] = np.mean(v1)

        #Caclulate and store mean MFE per i-j pair
        mean_mfe = {}
        for k1, v1 in mfe_dict.items():
            mean_mfe[k1] = np.mean(v1)

        #Caclulate and store mean ED per i-j pair
        mean_ed = {}
        for k1, v1 in ed_dict.items():
            mean_ed[k1] = np.mean(v1)

        #Caclulate and store total window counts per i-j pair
        total_windows = 0
        num_bp = 0
        for k1, v1 in zscore_dict.items():
            total_windows = total_windows + len(v1)
            key_data = re.split("-", str(k1))
            key_i = str(key_data[1])
            #print(k, k1, key_i)
            if int(k) == int(key_i):
                continue
            if int(k) != int(key_i):
                num_bp += 1

        #Print first line of log file tables (first half of log file)
        k_nuc = str((nuc_dict[k].nucleotide))
        #print(k_nuc)
        log_total.write("\ni-nuc\tBP(j)\tNuc\t#BP_Win\tavgMFE\tavgZ\tavgED"
            +"\tSumZ\tSumZ/#TotalWindows\tBPs= "+str(num_bp)+"\n")
        log_total.write("nt-"+str(k)+"\t-\t"+str(k_nuc)+"\t"+str(total_windows)
            +"\t-\t-\t-\t-\t-"+"\n")

        #Print remainder of log file tables (first half of log file)
        total_window_mean_z = {}
        for k1, v1 in zscore_dict.items():
            bp_window = str(len(v1))
            key_data = re.split("-", str(k1))
            key_nuc = str(key_data[0])
            key_i = str(key_data[1])
            total_window_mean_z[k1] = (np.sum(v1))/total_windows
            z_sum = str(round(np.sum(v1), 2))
            z_avg = str(round(np.mean(v1), 2))
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
                log_total.write(str(k)+"\t"+key_i+"\t"+key_nuc+"\t"+bp_window+"\t"
                    +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                    +z_sum+"\t"+str(test)+"\n")

        #Define best_bp_key based on coverage-normalized z-score
        best_bp_key = min(total_window_mean_z, key = total_window_mean_z.get)

        #Access best i-j NucPairs for each metric using best_bp_key
        best_bp_mean_z = mean_z[best_bp_key]
        best_bp_sum_z = sum_z[best_bp_key]
        best_bp_mean_mfe = mean_mfe[best_bp_key]
        best_bp_mean_ed = mean_ed[best_bp_key]
        best_total_window_mean_z = total_window_mean_z[best_bp_key]

        #Access best i-j pair info from key name
        best_bp_data = re.split("-", best_bp_key)
        best_nucleotide = best_bp_data[0]
        best_coordinate = best_bp_data[1]

        #Fill dictionary with coverage normalized z-score
        #print("Determining best base pair for nucleotide ", k)
        best_total_window_mean_bps[k] = (NucPair((nuc_dict[k]).nucleotide,
                                        nuc_dict[k].coordinate, best_nucleotide,
                                        best_coordinate, best_total_window_mean_z,
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
    log_win.write("\ni\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED"
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
            if sum(test_k == int(v.jcoordinate) for v in best_bps.values()) >= 0:
                #print(start_coordinate, end_coordinate)
            ### Scan the entire dictionary:
                #keys = range(int(start_coordinate), int(end_coordinate))

            ### Scan two window's length flanking nucleotide:
                #print(len(best_bps))
                # print(length*4)
                if (len(best_bps) < length*4):
                    #print("success")
                    #Length of input less than length of flanks
                    keys = range(int(start_coordinate), int(end_coordinate))
                elif (
                    (v.icoordinate - length*(2)) >= int(start_coordinate) and
                    (v.icoordinate + (length*2)) <= int(end_coordinate)
                    ):
                    #print(str(v.icoordinate - length*(2)))
                    # print("MIDDLE")
                    keys = range(int(v.icoordinate-(length*2)), int(v.icoordinate+(length*2)))

                elif (
                    int(v.icoordinate + (length*(2))) <= int(end_coordinate)and
                    (v.icoordinate + (length*2)) <= int(end_coordinate)
                ):
                    #print("BEGINING"+str(v.icoordinate - (length*(2)))+" "+str(end_coordinate))
                    keys = range(int(start_coordinate), int(v.icoordinate+(length*2))+1)

                elif (v.icoordinate + (length*2)) >= int(end_coordinate):
                    if v.icoordinate-(length*2) > 0:
                        #print("END"+str(v.icoordinate + (length*2)))
                        keys = range(int(v.icoordinate-(length*2)), int(end_coordinate))
                    else:
                        keys =range(int(v.icoordinate-(length*2)), int(end_coordinate))

                else:
                    print("Sub-dictionary error")
                    raise ValueError("Sub-dictionary error")

                #print(keys)
                subdict = {k: best_total_window_mean_bps[k] for k in keys}
                # if k == 216:
                #     for subk, subv in subdict.items():
                #         print(subk, subv.icoordinate, subv.jcoordinate)
                #print("SubDict length for "+str(k)+"="+str(len(subdict)))

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
        if filter != None:
            write_dp(best_bps, output+str(filter)+".dp", filter)
        write_dp(best_bps, out1, float(10))
        write_dp(best_bps, out2, float(-1))
        write_dp(best_bps, out3, float(-2))
        write_dp(best_bps, output+"mean_"+str(round(meanz, 2))+".dp", meanz)
        write_dp(best_bps, output+"below_mean_"+str(round(one_sig_below, 2))+".dp", one_sig_below)
        print("ScanFold-Fold complete, find results in...")

    #Write CT files
    if competition == 1:
        print("Trying to write CT files with -c option")
        elapsed_time = str(round((time.time() - start_time), 2))+"s"
        print(elapsed_time)
        print("Writing CT files")

        # if filter != None or filter != -2:
        #     write_ct(final_partners, output+str(filter)+".ct", filter, strand)
        # write_ct(final_partners, output+"no_filter.ct", float(10), strand)
        # write_ct(final_partners, output+"-1.ct", float(-1), strand)
        # write_ct(final_partners, output+"-2.ct", float(-2), strand)
        # if filter != None or filter != -2:
        #     write_ct(final_partners, output+str(filter)+".ct", filter, strand)
        write_ct(final_partners, out1, float(10), strand)
        write_ct(final_partners, out2, float(-1), strand)
        write_ct(final_partners, out3, float(-2), strand)

        #Create a dbn file for forna
        makedbn(out1, "NoFilter", dbn_file_path1)
        makedbn(out2, "Zavg_-1", dbn_file_path2)
        makedbn(out3, "Zavg_-2", dbn_file_path3)
        # subprocess.run(['ct2dot', str(out1), '1', str(dbn_file_path1)])
        # subprocess.run(['ct2dot', str(out2), '1', str(dbn_file_path2)])
        # subprocess.run(['ct2dot', str(out3), '1', str(dbn_file_path3)])


    if competition == 1:
        write_bp(final_partners, out6, start_coordinate, strand, length)
        write_wig_dict(final_partners, final_partners_wig, name, strand)
    if competition == 0:
        write_bp(best_bps, out6, start_coordinate, strand, length)
    write_fasta(nuc_dict, out7, name)
    write_fai(nuc_dict, fasta_index_path, name)

    #create output file with all DBNs and RNAfold with constraints
    dbn_log_file = open(dbn_file_path, "w+")

    #generate full fasta sequence as a string
    full_fasta_sequence = str()
    for k, v in nuc_dict.items():
        nucleotide = v.nucleotide
        full_fasta_sequence += nucleotide

    dbn_files_to_combine = []

    #If specified, refold and generate global fold of input sequence
    if global_refold == True:
        #fold the full fasta input as a fold compound (full_fc) using model params (md)
        print("Refolding full sequence using ScanFold results as constraints...")
        elapsed_time = round((time.time() - start_time), 2)
        print("Elapsed time: "+str(elapsed_time)+"s")

        #refold from -1 constraints
        #refolded_filter1_structure, _, refolded_filter1_MFE, _ = rna_refold(full_fasta_sequence, int(temperature), dbn_file_path2)

        #refold from -2 constraints
        refolded_filter2_structure, _, refolded_filter2_MFE, _ = rna_refold(full_fasta_sequence, int(temperature), dbn_file_path3)

        #extract the structure
        #full_structure, _, full_MFE, _ = rna_fold(full_fasta_sequence, int(temperature))

        #dbn_log_file.write(">"+str(name)+"\tGlobal Full MFE="+str(full_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(full_structure)+"\n")
        #dbn_log_file.write(">"+str(name)+"\Refolded with -1 constraints MFE="+str(refolded_filter1_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(refolded_filter1_structure)+"\n")
        dbn_log_file.write(">"+str(name)+"\Refolded with -2 constraints MFE="+str(refolded_filter2_MFE)+"\n"+str(full_fasta_sequence)+"\n"+str(refolded_filter2_structure)+"\n")
        dbn_log_file.close()
        dbn_files_to_combine = [str(dbn_file_path), str(dbn_file_path1), str(dbn_file_path2), str(dbn_file_path3)]
    else:
        dbn_files_to_combine = [str(dbn_file_path1), str(dbn_file_path2), str(dbn_file_path3)]

    with open(str(dbn_file_path4), 'w') as outfile:
        for fname in dbn_files_to_combine:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    #############
    #Begin the structure extract process

    bp_dict = {}

    #Set flanking nucleotides to be folded
    flanking = 0

    #Set number of randomizations and shuffle type ("mono" or "di")
    randomizations = 100
    shuffle = "mono"


    #Inititate base pair tabulation variables
    bond_order = []
    bond_count = 0

    dbn_file_filter2 = open(dbn_file_path3, "r")
    lines = dbn_file_filter2.readlines()
    filter2constraints = str(lines[2])

    #Read the structure of -2 filter2constraints
    if global_refold == False:
        #refold from -2 constraints
        full_fasta_sequence = str(lines[1])

    structure_raw = filter2constraints
    sequence_raw = full_fasta_sequence

    sequence = list(sequence_raw)
    #print(sequence)
    structure = list(structure_raw)
    #print(structure)
    length = len(sequence)
    #print(length)
    length_st = len(structure)
    #print(length_st)


    #Iterate through sequence to assign nucleotides to structure type
    m = 0
    nuc_dict = {}
    while  m < length:
        if structure[m] == '(':
            #print(m, structure[m])
            bond_count += 1
            bond_order.append(bond_count)
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1

        elif structure[m] == ')':
        #    print(m, structure[m])
            bond_order.append(bond_count)
            bond_count -= 1
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1

        elif str(structure[m]) == ( '.' ):
        #    print(m, structure[m])
            bond_order.append(0)
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '<' ):
        #    print(m, structure[m])
            bond_order.append(0)
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '>' ):
        #    print(m, structure[m])
            bond_order.append(0)
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '{' ):
        #    print(m, structure[m])
            bond_order.append(0)
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        elif str(structure[m]) == ( '}' ):
        #    print(m, structure[m])
            bond_order.append(0)
            nuc_dict[m] = NucStructure(bond_count, (m+1), sequence[m], structure[m])
            m += 1
        else:
            #print("Error", bond_count, (m+1), sequence[m], structure[m])
            m += 1
            continue
            # print("no")


    #print(bond_order)
    #Initiate base_pair list
    base_pairs = []

    #Create empty variable named test
    test = ""

    #Iterate through bond order
    j = 0
    structure_count = 0
    structure_end = []
    structure_start = []
    while j < length:
        try:
            if (nuc_dict[j].bond_order == 1) and (nuc_dict[j].structure == '('):
                structure_count += 1
                #print(nuc_dict[j].structure)
                structure_start.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict[j].nucleotide, nuc_dict[j].structure))
                j += 1

            elif (nuc_dict[j].bond_order == 0) and (nuc_dict[j].structure == ')'):
                structure_count += 1
                #print(nuc_dict[j].structure)
                structure_end.append(NucStructure(structure_count, nuc_dict[j].coordinate, nuc_dict[j].nucleotide, nuc_dict[j].structure))
                j += 1
            else:
                j += 1
        except:
            j += 1
            continue

    #print(structure_start[0].coordinate, structure_end[0].coordinate)

    #print(len(structure_start))

    l = 0
    extracted_structure_list = []
    while l < int(len(structure_start)):
        offset = flanking
        s = structure_start_coordinate =  int((structure_start[l].coordinate)-offset-1)
        e = structure_end_coordinate = int((structure_end[l].coordinate)+offset-1)

        seq = ""
        fold = ""
        for k, v in nuc_dict.items():
            if s <= k <= e:
                seq += str(v.nucleotide)
                fold += str(v.structure)

        extracted_structure_list.append(ExtractedStructure(l, seq, fold, s, e))

        l += 1

    #print(len(extracted_structure_list))

    zscore_total = []
    numerical_z = []
    pscore_total = []
    numerical_p = []
    MFE_total = []
    ED_total = []

    with open(structure_extract_file, "w") as se:
        se.write("ScanFold predicted structures which contain at least one base pair with Zavg < -2 have been extracted from "+str(name)+" results (sequence length "+str(length)+"nt) and have been refolded using RNAfold to determine their individual MFE, structure, z-score (using 100X randomizations), and ensemble diversity score.\n")
        for i in extracted_structure_list[:]:
            frag = i.sequence
            frag = frag.upper()
            frag = transcribe(frag)
            # fc = RNA.fold_compound(str(frag)) #creates "Fold Compound" object
            # fc.pf() # performs partition function calculations
            # frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment
            # (MFE_structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
            # MFE = round(MFE, 2)
            # MFE_total.append(MFE)
            # (centroid, distance) = fc.centroid() # calculate and define variables for centroid
            # ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
            # ED_total.append(ED)            #print(structure)

            MFE_structure, centroid, MFE, ED = rna_fold(str(frag), temperature)
            MFE = round(MFE, 2)
            MFE_total.append(MFE)
            ED = round(ED, 2) # this caclulates ED based on last calculated partition funciton
            ED_total.append(ED)            #print(structure)


            #fmfe = fc.pbacktrack()
            #print(str(fmfe))
            seqlist = [] # creates the list we will be filling with sequence fragments
            seqlist.append(frag) # adds the native fragment to list
            scrambled_sequences = scramble(frag, 100, type)
            seqlist.extend(scrambled_sequences)
            energy_list = energies(seqlist, temperature)
            try:
                zscore = round(zscore_function(energy_list, 100), 2)
            except:
                zscore = zscore_function(energy_list, 100)
            zscore_total.append(zscore)

            pscore = round(pscore_function(energy_list, 100), 2)
            #print(pscore)
            pscore_total.append(pscore)


            se.write("Structure number "+str(int(i.structure_count)+1)+"\n")
            se.write("Coordinates ="+str(int((i.i)+1))+' to '+str(int((i.j)+1))+"\n")
            se.write(str(i.sequence)+"\n")
            se.write(str(i.structure)+"\n")
            se.write(MFE_structure+"\n")
            se.write('z-score= '+str(zscore)+"\n")
            se.write('ED='+str(ED)+"\n")
            se.write('MFE='+str(MFE)+"\n")
            se.write('\n')

    dbn_log_file.close()
    print("ScanFold-Fold analysis complete! Refresh page to ensure proper loading of IGV")
    #print(url)
