# ScanFold
*UPDATE! ScanFold has been implemented as a WebServer.* Try it out at https://mosslabtools.bb.iastate.edu/scanfold

The ScanFold pipeline is a set of scripts which scan a large RNA sequence (using ScanFold-Scan.py) and subsequently extract  structural motifs (using ScanFold-Fold.py) which have evidence of being ordered by evolution to form an unusually stable structure (potentially to serve a functional role).  


The ScanFold-Scan.py and ScanFold-Fold.py scripts were built using Python3.6 and utilize several outside python modules.

ScanFold-Scan dependencies:
1. ViennaRNA; Must ensure that python3.6 can import "RNA" as a module. 
  See details for that here: https://www.tbi.univie.ac.at/RNA/documentation.html#install
2. NumPy (www.numpy.org/)
3. SeqIO from BioPython (https://biopython.org/wiki/SeqIO)
 
ScanFold-Fold dependencies:
1. RNAStructure (https://rna.urmc.rochester.edu/Overview/Python.html)

