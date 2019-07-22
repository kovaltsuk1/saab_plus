[![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)](https://www.python.org/download/releases/2.7/)


# Structural Annotation of BCR Repertoires
SAAB+: annotation of BCR Repertoires with structural information, using separate tools
to map non-CDR-H3 and CDR-H3 loop sequence to repersentative antibody structures.

## Getting Started
Please follow instructions below to install SAAB+ pipeline on our machine

### Prerequisites
SAAB+ comes with the antibody customized version of **FREAD** package.  
**anarci** and **scalop** are not supplied and need to be downloaded separately.

> * [anarci v1.3](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci) - is the antibody sequence numbering tool. In SAAB+ pipeline, anarci is used to filter the sequences for structural viability. Upon installation, anarci automatically downloads the latest version of IMGT germlines.
> * [scalop](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/scalop) - is the antibody canonical class annotation tool. SCALOP is downloaded with the latest version of antibody canonical classes.
> * [hmmer](http://hmmer.org/download.html) - is required to build antibody HMM profiles that anarci uses. anarci was tested with 3.1b2 hmmer version.

### Installing

To install SAAB+ as root run:

```
python setup.py install
```
For users without root access install locally using:

```
python setup.py install --user
```
## Checking installation
To check if installation was successful, simply run
```
SAAB_PLUS_DIAG
```
This script generates a file __diagnostics.log__.  
If installation was successful, the diagnostics.log should look like:
```
2019-07-18 16:21:55,780 INFO 	Writing DIAGNOSTICS log
2019-07-18 16:21:55,943 INFO  Successfully imported: anarci
2019-07-18 16:21:55,944 INFO  Successfully imported: anarci germlines
2019-07-18 16:21:55,946 INFO  Successfully imported: anarci Accept
2019-07-18 16:21:56,244 INFO  Successfully imported: scalop
2019-07-18 16:21:56,291 INFO  Successfully imported: FREAD
2019-07-18 16:21:56,291 INFO  Successfully imported: FREAD ESS table
2019-07-18 16:21:56,291 INFO  Successfully imported: prosci module
2019-07-18 16:21:56,295 INFO  Successfully imported: Common module
2019-07-18 16:21:56,295 INFO          PDB template info file is located: OK
2019-07-18 16:21:56,331 INFO   Directory with PDB frameworks is located: OK
2019-07-18 16:21:56,332 INFO  Directory with FREAD templates is located: OK
2019-07-18 16:21:56,658 INFO  Number of FREAD templates found: 3759
2019-07-18 16:21:56,658 INFO  Directory with numbered PDB frameworks: OK
```

### Running a test example
SAAB+ takes antibody sequences in the fasta format as the input. e.g.
>&gt;seq1  
>SLRLSCAASGFTFSGHWMYWVRQAPGKGLVWVARINND.....  
>&gt;seq2  
>SLRLSCAASGFTFRSYWMSWVRQAPGRGLEWIARIND......

The EXAMPLE folder contains a test fasta file.
To run SAAB+ pipeline on human BCR data across twenty CPU cores, 
```
SAAB_PLUS -f Fasta_example.fa -n 20 -o human
```
### Interpreting SAAB+ outputs
SAAB+ outputs a zipped tab-delimited text file.
```
          Protein_Seq - Full antibody amino acid sequences. Only sequences that passed anarci structural viability assessment are retained.  
                H3pdb - CDR-H3 structure that was predicted by FREAD  
                Canon - SCALOP predicted non-CDR-H3 canonical classes  
           Redundancy - Number of Proiten_Seq copies in the input fasta file  
            Framework - PDB structure that was used in CDR-H3 structure prediction  
              CDRHSeq - Sequence of CDR-H3 loop  
                  ESS - FREAD score for the predicted CDR-H3 structure  
StructurallyAnnotated - Sequnces, whose FREAD CDR-H3 structure prediction scores were above the quality threshold      
```
