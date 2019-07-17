import pandas as pd
import os, sys
from argparse import ArgumentParser
from aboss_utils.run_ANARCI_parsing import Main
import itertools
from Bio import SeqIO
from collections import Counter
from RunSAAB import RunSaab
from code.apply_ess_cutoff import ApplyESS
import logging

logging.basicConfig(filename='saab_plus.log', level=logging.INFO,
                     format="%(asctime)-15s %(levelname)s %(message)s",
                     datefmt='%Y-%m-%d %H:%M:%S')

class Filtering(object):
    """
    Anootation of Ig-seq data with structural information
    """
    def __init__(self, arguments):
        self.filename = arguments["filename"]
        self.ncores = arguments["ncores"]
        self.species = arguments["organism"]
        self.chain = arguments["chain"]
        self.output_name = arguments["output_name"]
        self.output_dir = arguments["output_dir"]
        print "SAAB plus log is found in saab_plus.log"

    def unique_aa(self):
        "reading fasta and converting to python dict"
        records = SeqIO.parse(self.filename, 'fasta')
        return Counter(str(rec.seq) for rec in records)

    def convert_fasta_to_frame(self):
        "converting fasta file to pandas dataframe"
        df = pd.DataFrame( self.unique_aa().items() , columns=["Protein_Seq", "Redundancy"])
        self.gziped = "{0}_df.txt.gz".format(self.filename)
        df.to_csv(self.gziped, sep = "\t", compression="gzip")

    def run(self):
        # Starting anarci parsing
        print "\n\t\tInitiating anarci parsing" 
        logging.info("\tStarting ANARCI parsing\n")
        anarci_parsing = Main(input_file = self.gziped, 
                             ncpu = self.ncores, 
                             chain = self.chain, 
                             species = self.species, 
                             output_dir =".", 
                             output_name = "anarci_parsed.txt") 
        anarci_parsing.run()
        os.remove(self.gziped)

        print "\n\t\tanarci parsing is done\n"
        logging.info("\tANARCI parsing is finished")
        # Starting SAAB running
        anarci_gzip_file = "anarci_parsed.txt.gz"
        if not os.path.isfile(anarci_gzip_file):
            logging.error("anarci parsed txt file is not found: {0}".format(anarci_gzip_dir))
            raise AssertionError("File is not found: ", anarci_gzip_file)

        # Running SAAB+
        print "\t\tInitiating SAAB plus"
        logging.info("\tStarting SAAB plus pipeline\n")
        saab_instance = RunSaab(anarci_gzip_file, 
                                self.output_name, 
                                self.output_dir) 

        saab_instance.initiate_SAAB(self.ncores)
        os.remove(anarci_gzip_file)

        logging.info("\tSAAB+ is finished")
        logging.info("\tApplying quality filters")
        apply_ess_instance = ApplyESS(os.path.join(self.output_dir, self.output_name+".gz") )
        apply_ess_instance.adding_pdb_length()
        apply_ess_instance.apply_ess()
        apply_ess_instance.save()
        logging.info("\touput_file is: {0}".format(os.path.join(self.output_dir, self.output_name+".gz")))
        logging.info("\n\n\t\tDONE!\n\n")
        print "\n\t\tDone!\n"

def initArgParser():
    parser = ArgumentParser(description=__doc__)
    group = parser.add_argument_group('study arguments')
    group.add_argument('-f', action='store', dest='filename', default=False,
                        help="Filename to be analyzed")
    group.add_argument('-s', action='store', dest='organism', default=False,
                        help="organism information in the study [human, mouse, rhesus, alpaca]")
    group.add_argument('-n', action='store', dest='ncores', default=16, type=int,
                        help="number of cores to use in multiprocessing")
    group.add_argument("-c", action="store", dest="chain", default="H",
                        help="Chain information H (currently only H) ")
    group.add_argument("-o", action="store", dest="output_name", default="output.txt",
                        help="output filename, Default is output_name.txt")
    group.add_argument("-d", action="store", dest="output_dir", default=".",
                        help="output filename, Default is current working directory")
    return parser

def printStatus(args):
    print "      FILE: {0}".format(args.filename)
    print "     CORES: {0}".format(args.ncores)
    print "  ORGANISM: {0}".format(args.organism)
    print "     CHAIN: {0}".format(args.chain)
    print "OutputName: {0}".format(args.output_name)
    print " OutputDir: {0}".format(args.output_dir)

    logging.info("\t\tInitiating structural annotation pipeline")
    logging.info("\t      FILE: {0}".format(args.filename))
    logging.info("\t     CORES: {0}".format(args.ncores))
    logging.info("\t  ORGANISM: {0}".format(args.organism))
    logging.info("\t     CHAIN: {0}".format(args.chain))
    logging.info("\tOutputName: {0}".format(args.output_name))
    logging.info("\t OutputDir: {0}".format(args.output_dir))

    return vars(args)

if __name__ == "__main__":

    parser = initArgParser()
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if not args.filename:
        raise AssertionError("Study was not specified, exiting!!!")

    if not args.organism:
        raise AssertionError("Organism was not selected")

    arguments = printStatus(args)
    STUDY = Filtering(arguments)
    STUDY.convert_fasta_to_frame()
    STUDY.run()
