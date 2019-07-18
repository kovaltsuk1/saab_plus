# Structural Annotation of BCR Receptor Repertoires
SAAB+: annotation of BCR Receptor Repertoires with structural information, using separate tools
to map non-CDR-H3 and CDR-H3 loop sequence to repersentative structures.

## Getting Started
Please follow instructions below to install SAAB+ pipeline on our machine

### Prerequisites
SAAB+ comes with **FREAD** package, however, **anarci** and **scallop** are not supplied.

> * [anarci v1.3](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci) - anarci is neccessary to number BCR Repertoire sequences and
                             filter the sequences for structural viability.
                             anarci automatically downloads the latest IMGT germlines.
> * [scalop](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/scalop) - SCALOP is antibody canonical class annotation software.
                        SCALOP is downloaded with the latest number of canonical classes

### Installing
To install SAAB+, run
```
python setup.py install --user
```
## Checking installation
To check if installation was successful, simply run
```
SAAB_PLUS_DIAG
```
This script generates a file __diagnostics.log__
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

