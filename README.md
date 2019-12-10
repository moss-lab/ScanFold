# ScanFold
*UPDATE! ScanFold has been implemented as a WebServer.* Try it out at https://mosslabtools.bb.iastate.edu/scanfold
and read more about its uses here: https://doi.org/10.1016/j.ymeth.2019.11.001

*UPDATE (12/10/2019)! ScanFold-Scan can now consider hard constraints as input during the scanning window process.* Constraints should be foratted as .dbn files, where line 1 is a header, line 2 is the sequence and line three contains the hard constraints.

Here's a list of allowed constraints:

"." : no constraint at all
"x" : base must not pair (forced single stranded)
"|" : paired with another base (undefined location up/downstream)
">" : base i is paired with a base j>i (downstream)
"<" : base i is paired with a base j<i (upstream)
matching brackets ( ): base i pairs with base j (defines exact base paired nt)

The ScanFold pipeline is a set of scripts which scan a large RNA sequence (using ScanFold-Scan.py) and subsequently extract  structural motifs (using ScanFold-Fold.py) which have evidence of being ordered by evolution to form an unusually stable structure (potentially to serve a functional role).  


The ScanFold-Scan.py and ScanFold-Fold.py scripts were built using Python3.6 and utilize several outside python modules.

ScanFold-Scan dependencies:
1. ViennaRNA; Must ensure that python3.6 can import "RNA" as a module. 
  See details for that here: https://www.tbi.univie.ac.at/RNA/documentation.html#install
2. NumPy (www.numpy.org/)
3. SeqIO from BioPython (https://biopython.org/wiki/SeqIO)
 
ScanFold-Fold dependencies:
1. RNAStructure (https://rna.urmc.rochester.edu/Overview/Python.html)

