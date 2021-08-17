# ScanFold

*ScanFold has been implemented as a WebServer.* Try it out at https://mosslabtools.bb.iastate.edu/scanfold
and read more about its uses here: https://doi.org/10.1016/j.ymeth.2019.11.001

The ScanFold pipeline is a set of scripts which scan a large RNA sequence (using ScanFold-Scan.py) and subsequently extract  structural motifs (using ScanFold-Fold.py) which have evidence of being ordered by evolution to form an unusually stable structure (potentially to serve a functional role).

## Updates
*(Update from 2/17/2021) ScanFold-Scan and ScanFold-Fold have been combined into a single script called ScanFold.py. This script can incorporate soft constraints as pseudo energies using RNAfold's implementation of the Deigan or Zarringhalam algorithms.*

*(Update from 12/10/2019) ScanFold-Scan can now consider hard constraints as input during the scanning window process. Constraints should be formatted as .dbn files, where line 1 is a header, line 2 is the sequence and line three contains the hard constraints.*

  Here's a list of allowed constraints (https://www.tbi.univie.ac.at/RNA/RNAfold.1.html#heading6):

## Dependencies 
The ScanFold-Scan.py and ScanFold-Fold.py scripts were built using Python3.6 and utilize several outside python modules.

ScanFold-Scan dependencies:
1. ViennaRNA; Must ensure that python3.6 can import "RNA" as a module. 
  See details for that here: https://www.tbi.univie.ac.at/RNA/documentation.html#install
2. NumPy (www.numpy.org/)
3. SeqIO from BioPython (https://biopython.org/wiki/SeqIO)
 
ScanFold-Fold dependencies:
1. RNAStructure (https://rna.urmc.rochester.edu/Overview/Python.html)

ScanFold.py dependencies:
1. ViennaRNA
2. SeqIO from BioPython
