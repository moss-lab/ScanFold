# ScanFold

The ScanFold pipeline is a set of scripts which scan (using ScanFold-Scan.py) a large RNA sequence and subsequently extract (using ScanFold-Fold.py) structural motifs which have evidence of being ordered by evolution (potentially to serve a functional role).  

The ScanFold-Scan.py and ScanFold-Fold.py scripts were built using Python3.6 and utilize several outside python modules.

Dependencies (for ScanFold-Scan):
 ViennaRNA; Must ensure that python3.6 can import "RNA" as a module. See details for that here: https://www.tbi.univie.ac.at/RNA/documentation.html#install
 NumPy (www.numpy.org/)
 SeqIO from BioPython (https://biopython.org/wiki/SeqIO)
 
Dependencies (for ScanFold-Fold):
 RNAStructure (https://rna.urmc.rochester.edu/Overview/Python.html)
