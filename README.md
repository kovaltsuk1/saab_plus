# Structural Annotation of BCR Receptor Repertoires
SAAB+: annotation of BCR Receptor Repertoires with structural information, using separate tools
to map non-CDR-H3 and CDR-H3 loop sequence to repersentative structures.

## Getting Started
Please follow instructions below to install SAAB+ pipeline on our machine

### Prerequisites
SAAB+ comes with **FREAD** package, however, **anarci** and **scallop** are not supplied.

> * [anarci v1.3](http://) - anarci is neccessary to number BCR Repertoire sequences and
                             filter the sequences for structural viability.
                             anarci will download the latest IMGT germlines.
> * [scalop](http://) - SCALOP is antibody canonical class annotation software.
                        SCALOP will be downloaded with the latest number of canonical classes

### Installing
To install SAAB+, run
> python setup.py install --user

## Checking installation
To chec

