
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import subprocess
import numpy as np
import random
import re
from multiprocessing import get_context
import os

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
      try:
        a = edge[0]; b = edge[1]
        #print(a, b)
        #print(edge)
        if D[b]==1: D[a]=1
      except:
        continue
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


def rna_fold(frag, temperature):
    args = ["RNAfold", "-p", "-T", str(temperature)]
    fc = subprocess.run(args, input=str(frag), check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = str(fc.stdout)
    test = out.splitlines()
    structure = test[1].split()[0]
    centroid = test[3].split()[0]
    MFE = test[1].split(" ", 1)[1]
    try:
        MFE = float(re.sub('[()]', '', MFE))
    except:
        print("Error parsing MFE values", test)
    ED = float(test[4].split()[-1])

    return (structure, centroid, MFE, ED)

def rna_refold(frag, temperature, constraint_file):
    args = ["RNAfold", "-p", "-T", str(temperature), '-C', constraint_file]
    fc = subprocess.run(args, input=frag, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = str(fc.stdout)
    test = out.splitlines()
    structure = test[2].split()[0]
    centroid = test[4].split()[0]
    MFE = test[2].split(" ", 1)[1]
    try:
        MFE = float(re.sub('[()]', '', MFE))
    except:
        print("Error parsing MFE values", test)
    ED = float(test[5].split()[-1])

    return (structure, centroid, MFE, ED)

def rna_folder(arg):
    (frag, temperature) = arg
    _, _, MFE, _ = rna_fold(frag, temperature)
    return MFE

###### Function to calculate MFEs using RNAfold #################
def energies(seq_list, temperature):
    energy_list = []

    energy_list = multiprocessing(rna_folder, [(sequence, temperature) for sequence in seq_list], 4)
    # for sequence in seq_list:
    #     #fc = RNA.fold_compound(str(sequence))
    #     (structure, MFE) = RNA.fold(str(sequence)) # calculate and define variables for mfe and structure
    #     energy_list.append(MFE) # adds the native fragment to list

    return energy_list;

def multiprocessing(func, args,
                    workers):

    if 'SCANFOLDMPUSETHREADS' in os.environ:
        with ThreadPoolExecutor(workers) as ex:
            res = ex.map(func, args)
    else:
        ctx = get_context()
        if 'SCANFOLDMPMETHOD' in os.environ:
            ctx = get_context(os.environ['SCANFOLDMPMETHOD'])
        with ProcessPoolExecutor(workers, ctx) as ex:
            res = ex.map(func, args)

    return list(res)

def makedbn(ctfile, name, dbnfullname):
    icoord = ()
    jcoord = ()
    kcoord = ()
    lcoord = ()
    #dbnfullname = ctfile+".dbn"
    sequence = ""
    with open(ctfile,'r') as ct1:
        dot = ""
        data = ct1.readlines()[1:]
        for line in data:
            rows = line.split()
            icoord = int(rows[0])
            jcoord = int(rows[-2])
            if len(rows[1]) > 1:
                print("skipping header")
                continue

            elif int(rows[-2]) == 0:
                dot += '.'
                sequence += rows[1]
            else:
                sequence += rows[1]
                if icoord < jcoord:
                    with open(ctfile,'r') as ct2:
                        data = ct2.readlines()
                        for line in data[icoord:]:
                            rows = line.split()
                            kcoord = int(rows[0])
                            lcoord = int(rows[-2])

                            if kcoord == jcoord:
                                #print(kcoord, icoord, "(")
                                dot += '('
                                break
                            else:
                                if lcoord == 0:
                                    pass
                                elif lcoord < icoord:
                                    print('Non-nested pair found: '+str(lcoord)+' to '+str(kcoord))
                                    dot += '<'
                                    break
                                else:
                                    pass
                elif icoord > jcoord:
                    with open(ctfile,'r') as ct2:
                        data = ct2.readlines()
                        for line in data[jcoord:]:
                            rows = line.split()
                            kcoord = int(rows[0])
                            lcoord = int(rows[-2])
                            if kcoord == icoord:
                                #print(kcoord, icoord, ")")
                                dot += ')'
                                break
                            else:
                                if lcoord == 0:
                                    pass
                                elif lcoord < jcoord:
                                    dot += '>'
                                    break
                                else:
                                    pass
                else:
                    print('Error in non-nested pair search'+'\n'+'i = '+str(icoord)+'\n'+'j = '+str(jcoord)+'\n'+'k = '+str(kcoord)+'\n'+'l = '+str(lcoord))
    #print(sequence)
    #print(dot)
    with open(dbnfullname,'w') as dbn:
        dbn.writelines(f">{name}\n")
        dbn.writelines(f"{sequence}\n")
        dbn.writelines(f"{dot}\n")
