#!/usr/bin/python3.6
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
$ python3.6 ScanFold-Fold.py 1-input 2-filter > 3-logfile

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
sys.path.append('/Users/ryanandrews/Desktop/programs/RNAstructure/exe')
import RNAstructure

filename = sys.argv[1]
filter = float(sys.argv[2])

output_data = re.split('\.', str(filename))
output = str(str(output_data[0])+".ScanFold.")
log_total = open(output+"log", 'w')
log_win = open(output+"final_partners.log", 'w')

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
            sequence = transcribe(data[8])
            zscore = float(data[4])
            for nuc in sequence:
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

    v = pair_dict[best_bp_key]
    best_bp = v[0]

    return best_bp;

def write_ct(base_pair_dictionary, filename, filter):
    #Function to write connectivity table files from a list of best i-j pairs
    w = open(filename, 'w')
    w.write((str(len(base_pair_dictionary))+"\t"+filename+"\n"))
    for k, v in base_pair_dictionary.items():
        if float(v.zscore) < filter:
            if ((int(v.icoordinate) < int(v.jcoordinate)) and (int(v.icoordinate) == int(k))): #test to see if reverse bp.
                w.write("%d %s %d %d %d %d\n" % (int(k), v.inucleotide, int(k)-1, int(k)+1, int(v.jcoordinate), int(k)))

            elif ((int(v.icoordinate) > int(v.jcoordinate)) and (int(v.icoordinate) == int(k))): #test to see if reverse bp.
                w.write("%d %s %d %d %d %d\n" % (int(k), v.inucleotide, int(k)-1, int(k)+1, int(v.jcoordinate), int(k)))

            elif (int(v.icoordinate) < int(v.jcoordinate)) and (int(k) == int(v.jcoordinate)):
                w.write("%d %s %d %d %d %d\n" % (int(k), v.jnucleotide, int(k)-1, int(k)+1, int(v.icoordinate), int(k)))

            elif (int(v.icoordinate) > int(v.jcoordinate)) and (int(k) == int(v.jcoordinate)):
                w.write("%d %s %d %d %d %d\n" % (int(k), v.jnucleotide, int(k)-1, int(k)+1, int(v.icoordinate), int(k)))

            elif int(v.icoordinate) == int(v.jcoordinate):
                w.write("%d %s %d %d %d %d\n" % (int(k), v.inucleotide, int(k)-1, int(k)+1, 0, int(k)))
            #
            # elif (int(k) != v.icoordinate) and (int(k) != int(v.jcoordinate)):
            #     continue
            #     #w.write("%d %s %d %d %d %d\n" % (int(k), v.inucleotide, int(k)-1, int(k)+1, 0, int(k)))
            else:
                print("Error at", int(k), v.inucleotide, v.icoordinate, v.jnucleotide, int(v.jcoordinate), v.zscore)
        else:
            if int(k) == int(v.icoordinate):
                w.write("%d %s %d %d %d %d\n" % (int(k), v.inucleotide, int(k)-1, int(k)+1, 0, int(k)))
            elif int(k) == int(v.jcoordinate):
                w.write("%d %s %d %d %d %d\n" % (int(k), v.jnucleotide, int(k)-1, int(k)+1, 0, int(k)))
            else:
                raise ValueError("WriteCT function did not find a nucleotide to match coordinate (i or j coordinate does not match dictionary key)")
            continue

def transcribe(seq):
    #Function to covert T nucleotides to U nucleotides
    for ch in seq:
        rna_seq = seq.replace('T', 'U')
        return(rna_seq)

with open(filename, 'r') as f:
    #Initialize bp dictionary and z-score lists
    z_score_list = []
    bp_dict = {}

    #Read all lines from ScanFold-Scan file (exept header)
    lines = f.readlines()[1:]

    #Generate nucleotide dictionary to assign each nucleotide in sequence a key
    nuc_dict = NucleotideDictionary(lines)

    #Iterate through input file, read each rows metrics, sequence, etc.
    for row in lines:

        #Ignore blank lines
        if not row.strip():
            continue

        #Main loop to find all i-j pairs per i-nucleotide
        else:
            #Assign metrics to variables
            data = row.split('\t')
            icoordinate = data[0]
            jcoordinate = data[1]
            temp = data[2]
            mfe = float(data[3])
            zscore = float(data[4])
            pvalue = data[5]
            ed = float(data[6])
            fmfe = data[7]
            sequence_raw = transcribe(str(data[8]))
            structure_raw = data[9]
            centroid = data[10]

            #Convert sequence and structures into lists
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

#Initiate global dictionaries to store best base pairs
best_bps = {}
best_sum_bps = {}
best_sum_bps_means = {}
best_total_window_mean_bps = {}

#Iterate through initial i-nuc dictionary to determine best base pairs (round 1)
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
        if int(k) == int(key_i):
            continue
        if int(k) != int(key_i):
            num_bp += 1

    #Print first line of log file tables (first half of log file)
    k_nuc = str((nuc_dict[k].nucleotide))
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
            log_total.write(str(k)+"\tNoBP\t"+key_nuc+"\t"+bp_window+"\t"
                  +k1_mean_mfe+"\t"+z_avg+"\t"+k1_mean_ed+"\t"
                  +z_sum+"\t"+str(test)+"\n")
        else:
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

#print header for fianl partener log file (log_win)
log_win.write("\ni\tbp(i)\tbp(j)\tavgMFE\tavgZ\tavgED"
    + "\t*Indicates most favorable bp has more favorable partner or is "
    + "more likely to be unpaired (competing coordinates are reported)"+"\n")

#Iterate through round 1 i-j pairs
for k, v in sorted(best_bps.items()):

    #For each i and j in i-j pair, detect competing pairs and append to dict
    comp_pairs_i = competing_pairs(best_total_window_mean_bps, v.icoordinate)
    comp_pairs_j = competing_pairs(best_total_window_mean_bps, v.jcoordinate)
    total_pairs = []

    #Put pairs competing with i from i-j pair into total pair dict for i-nuc
    for key, pair in comp_pairs_i.items():
        total_pairs.append(competing_pairs(best_total_window_mean_bps,
                                           pair.jcoordinate))

    #Put pairs competing with j from i-j pair into total pair dict for i-nuc
    for key, pair in comp_pairs_j.items():
        total_pairs.append(competing_pairs(best_total_window_mean_bps,
                                           pair.jcoordinate))

    #Merge all dictionaries
    merged_dict = {}
    i = 0
    for d in total_pairs:
        for k1, v1 in d.items():
            merged_dict[i] = v1
            i += 1

    #initiate best_basepair fucntion, return best_bp based on sum
    bp = best_basepair(merged_dict, v.inucleotide, v.icoordinate, "sum")

    #Check if best basepair was connected to i-nucleotide (i.e., "k")
    if (int(k) != bp.icoordinate) and (int(k) != int(bp.jcoordinate)):
        #if there was a competing i-j pair print it to log file instead:
        log_win.write("nt-"+str(k)+"*:\t"+str(bp.icoordinate)+"\t"+bp.jcoordinate+"\t"
              +str(round(best_bps[k].mfe, 2))
              +"\t"+str(round(best_bps[k].zscore, 2))
              +"\t"+str(round(best_bps[k].ed, 2))+"\n")
        final_partners[k] = NucPair(v.inucleotide, v.icoordinate,
                                    v.inucleotide, v.icoordinate,
                                    best_bps[k].zscore,
                                    best_bps[k].mfe,
                                    best_bps[k].ed)
    else:
        #if there was no competing i-j pair, print to log file:
        log_win.write("nt-"+str(k)+":\t"+str(bp.icoordinate)+"\t"+bp.jcoordinate+"\t"
              + str(round(best_bps[k].mfe, 2))+"\t"
              + str(round(best_bps[k].zscore, 2))
              + "\t"+str(round(best_bps[k].ed, 2))+"\n")
        final_partners[k] = NucPair(bp.inucleotide, bp.icoordinate,
                                    bp.jnucleotide, bp.jcoordinate,
                                    best_bps[bp.icoordinate].zscore,
                                    best_bps[bp.icoordinate].mfe,
                                    best_bps[bp.icoordinate].ed)

#Write CT files
write_ct(final_partners, output+str(filter)+".ct", filter)
write_ct(final_partners, output+"no_filter.ct", float(10))
write_ct(final_partners, output+"-1.ct", float(-1))
write_ct(final_partners, output+"-2.ct", float(-2))
#write_ct(final_partners, output+"mean_"+str(round(meanz, 2))+".ct", meanz)
#write_ct(final_partners, output+"below_mean_"+str(round(one_sig_below, 2))+".ct", one_sig_below)

#Write DBN files from CT files
os.system(str("ct2dot "+output+"no_filter.ct 1 "+output+"no_filter.dbn"))
os.system(str("ct2dot "+output+"-1.ct 1 "+output+"-1.dbn"))
os.system(str("ct2dot "+output+"-2.ct 1 "+output+"-2.dbn"))
