Cola Sequence Aligner 
===================

Cola consists of an efficient implementation of a collection of sequence
alignment algorithms, extending the Smith-Waterman and Smith-Waterman-Gap-Affine
methods by the ability to apply a scoring function that is an arbitrary
function of the size of consecutive nucleotide matches.

If you want to integrate Cola with your code, that is simple. Look into src/runCola as
an example of how the Cola aligner object can be used.

- Getting started
  - Create a git clone by running git clone  https://github.com/nedaz/cola.git
  N.B. Make sure you use git to clone the repository and not use other methods such as svn checkout or download az zip

- Installation
  - Use gcc version 5 or higher and CMake higher than 3.5
  - Run ./configure in parent directory
  - Run make -C build -j 10
  - Binaries will be located in /bin directory


Once you have the executables created and will be able to follow the instructions below:

The main access point for a user to run cola is through runCola which accepts the following arguments:
-q  : Query sequence in FASTA format
-t  : Target sequence in FASTA format
-a  : Aligner type - Choose a number from 1 to 4 for the following modes - 1 : NSGA , 2 : NS , 3 : SWGA, 4 : SW 
-o  : Aligner gap open penalty
-m  : Aligner mismatch penalty
-e  : Aligner gap extension penalty
-all: align all to all sequences (default is false)
-p<double> : Maximum acceptable P-value (def=1)
-i<double> : Minium acceptable identity (def=0)
-b<int> : The bandwidth for banded mode, default is for unbanded (def=-1)

N.B. Only the first 3 arguments are compulsory and the rest are optional. 
If the aligner penalty parameters are not provided by user, defaults will be used. 

See following examples for more detail.

To run in one of 4 modes using the sample data provided in the sample directory:

1) Nonlinear Scoring Gap-Affine (NSGA): 
./runCola -t sample/homo.X.part.fa -q sample/canis.X.part.fa -a 1 

2) Nonlinear Scoring (NS): 
./runCola -t sample/homo.X.part.fa -q sample/canis.X.part.fa -a 2 

3) Smith-Waterman Gap-Affine (SWGA): 
./runCola -t sample/homo.X.part.fa -q sample/canis.X.part.fa -a 3 

4) Smith-Waterman (SW): 
./runCola -t sample/homo.X.part.fa -q sample/canis.X.part.fa -a 4 


See Cola manuscript for more detail on the underlying algorithm and methodology.

