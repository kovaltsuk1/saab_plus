import os, json
import pandas as pd
import multiprocessing as mp
from argparse import ArgumentParser
from itertools import tee
from aboss_utils.ProgressBar import returnProgress
from code.StructuralAlignment import align_single_sequence
from code.DataManagement.SAbDab import structural_reference
from code.PrepareForSAAB import oas_output_parser
from code.Common.Common import chain_id_converter
from collections import namedtuple
import logging

formatLoops = {"cdrh1":"H1",
               "cdrh2":"H2",
               "cdrh3":"H3Seq",
               "cdrl1":"L1",
               "cdrl2":"L2",
               "cdrl3":"L3Seq"} # Light chain for future work

chain_formatter = {"H": "Heavy",
                   "L": "Light"
                   }

sequence_obj = namedtuple("sequence_obj", "sequence numbering sequencemeta")

def find_redundancy(redundancy):
    """
    if redundancy is not integer
    return 1
    """
    try:
        return int(redundancy)
    except:
        return 1

def prepare_saab_data(sequence):
    """
    Processing data after anarci parsing.
    Preparing data for SAAB+
    ------------
    Parameters
        sequence - sequence object ( OAS database format )
    ------------
    Return
                 sequence.Sequence  - full (not-numbered) antibody sequence
        oas_output_parser(Numbered) - antibody sequence that is imgt numbered
                                      to comply with SAAB+ input format
                sequence_info_dict  - Dictionary that contains sequence metadata
                                      which is requeired for SAAB+ to run
    """

    cdr3sequence = sequence.CDRH3
    VGene = sequence.VGene[:5]
    Numbered = json.loads( sequence.Numbered )
    CDRs = [ loop for loop in Numbered.keys() if "cdr" in loop ]
    sequence_info_dict = { formatLoops[loop] : Numbered[loop] if "3" not in loop else cdr3sequence
                                                                                        for loop in CDRs  }
    sequence_info_dict["V"] = VGene
    sequence_info_dict["Redundancy"] = find_redundancy( sequence.Redundancy )
    return sequence_obj( sequence.Sequence, oas_output_parser(Numbered), sequence_info_dict )

class RunSaab(object):
    """
    SAAB+ main script
    ---------------
    Parameters
        input_file   - Anarci parsed and structurally filtered txt file
                       generated by run_ANARCI_parsing.py module
               ncpu  - Number of CPU
        output_name  - Name of the output file. This file will be gzipped
        output_dir   - Directory, where to save the output file
        self.strucs  - Dictionary that contains pre-numbered frameworks
                       with imgt numbering
        self.chain   - Antibody chain analysis (currently only heavy chain
                       is supported in SAAB+)
    """
    def __init__(self, input_file, ncpu,  output_name, output_dir, chain ):
        self.input_file = input_file
        self.ncpu = ncpu
        self.output_name = output_name
        self.output_dir = output_dir

        assert os.path.isfile( self.input_file ), "Input file does not exist"
        self.rows_to_skip = 0 # if anarci analysis was interrupted
        self.currently_analysed = 0
        self.chain = chain_formatter[chain]

        self.strucs = structural_reference( chain_id_converter[self.chain] )
        logging.info("\tNumbered PDBs have been loaded into memory")

    def create_output_dir(self):
        try:
            os.makedirs( self.output_dir )
            logging.info("Output Directory has been created: {0}".format( self.output_dir ))
        except:
            logging.info("Output Directory already exists")

    def initiate_SAAB(self):
        """
        Function that loads anarci parsed DataFrame.
        Then, we iterate through it in chunks of 60k sequences.
        Each chunk is distributed across self.ncpu.
        Intermediate Output is saved as a .txt file
        """
        if os.path.isfile( os.path.join(self.output_dir, self.output_name) ):
            self.rows_to_skip = self.to_skip()
            logging.info("SAAB annotated file already exists.\n \
                         Checking the number of rows to skip: {0}".format( self.rows_to_skip))
        # Creating an iterator and counting total number of entries
        iterator = pd.read_csv( self.input_file, header = None, iterator = True, chunksize = 60000, sep = ",", 
                                                                     names=[ "Sequence","Redundancy", "CDRH3", 
                                                                             "Aboss1", "Aboss2", "VGene", 
                                                                             "JGene", "Numbered", "Track"],
                                                                     skiprows = self.rows_to_skip)
        df_iterator, total_number_iterator = tee( iterator )
        # Checking the dataset size
        total_dataunit_size = sum([ chunk.shape[0]  for chunk in total_number_iterator ])
        if total_dataunit_size != 0:
            for chunk in df_iterator :
                data_for_saab = []
                for sequence in chunk.itertuples():
                    data_for_saab.append( prepare_saab_data( sequence ) )
                    self.currently_analysed += 1
                logging.info( returnProgress( float( self.currently_analysed ) / total_dataunit_size ) ) 
                pool = mp.Pool( processes = self.ncpu )
                results = [ pool.apply_async(align_single_sequence, args = ( data_for_saab[i::self.ncpu],
                                                                             self.strucs,
                                                                             self.chain )) 
                                                                             for i in xrange(self.ncpu) ]
                pool.close()
                self.writing_output(results)
                del results

        self.gzip_file()
        self.create_output_dir()
        os.rename("{0}.gz".format( self.output_name ), os.path.join( self.output_dir,"{0}.gz".format( self.output_name )))
        logging.info("{0} has been structurally annotated \n".format( os.path.join( self.output_dir, self.output_name ))) 

    def gzip_file(self):
        os.system('gzip -f '+self.output_name)

    def writing_output(self, results):
        """
        Writing SAAB+ outputs to a txt file
        """
        correct_sequences_portion = [ r.get() for r in results ]
        with open(self.output_name, "a") as txtFile:
            for portion in correct_sequences_portion:
                for seq in portion[0]:
                    seq_object = portion[0][seq]
                    txtFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(seq, 
                                                                               seq_object.CDR_H3_template, 
                                                                               seq_object.Canonical_classes, 
                                                                               seq_object.Redundancy, 
                                                                               seq_object.Framework_template,
                                                                               seq_object.CDR_H3_sequence,
                                                                               seq_object.ESS))
    def to_skip(self):
        """
        Loading already structurally annotated DataFrame
        and checking how many rows we need to skip
        """
        df = pd.read_csv(os.path.join( self.output_dir, self.output_name), header = None, sep = "\t", 
                                                                           names=[ "seq", "H3pdb", "Canon", 
                                                                                    "Redundancy", "FrameWork", 
                                                                                    "CDRHSeq", "ESS" ] )
        return len(df)

def initArgParser():

    parser = ArgumentParser(description=__doc__)
    group = parser.add_argument_group('input_file arguments')
    group.add_argument('-f', action='store', dest='filename', default=False,
                        help="ANARCI parsed file")
    group.add_argument('-n', action='store', dest='ncores', default=16, type=int,
                        help="number of cores to use in multiprocessing")
    group.add_argument("-o", action="store", dest="output_name", default="output_structure.txt",
                        help="Output output_name")
    group.add_argument("-d", action="store", dest="output_dir", default=".",
                        help="Output directory name")
    return parser

def printStatus(args):
    print "   FILENAME: {0}".format(args.filename)
    print "      CORES: {0}".format(args.ncores)
    print "Output FILE: {0}".format(args.output)
    print " Output Dir: {0}".format(args.output_dir)
    return vars(args)

if __name__ == "__main__":

    parser = initArgParser()
    args = parser.parse_args()

    input_file = args.filename
    ncpu = args.ncores
    output = args.output
    output_dir = args.output_dir
    arguments = printStatus(args)
    
    if not args.filename:
        raise AssertionError("Filename is not specified")

    saab_instance = RunSaab(input_file, ncpu, output, output_dir)
    saab_instance.initiate_SAAB()
