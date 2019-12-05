from anarci import run_anarci, all_germlines
import multiprocessing as mp
from os.path import join
import csv, sys, os
from itertools import tee
import numpy as np, pandas as pd
from saab_plus.aboss_utils.species_viability import checking_structural_viability, checking_sequence_viability,\
                                                                        list_of_frames, CDR3, FW1 
from argparse import ArgumentParser
csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saab_plus.aboss_utils.ProgressBar import returnProgress 
import logging

def check_anarci():

    # Recording ANARCI version used
    logging.info("\tchecking version of anarci used")
    if len(all_germlines["V"]["H"]["human"]) >= 203:
        logging.info("\tANARCI version 1.3 is used")
    elif len(all_germlines["V"]["H"]["human"]) == 201:
        logging.info("\tANARCI version 1.2 is used")
    else:
        logging.info("\tNot sure about anarci version")

def anarci_fun(list_of_inputs, ncores):
    """
    ANARCI running function
    ---------------
    Parameters
        list_of_inputs - List of antibody sequences
                ncores - Number of CPU cores to use
    ---------------
    Return
         output - ANARCI annotated and numbered sequences
    """
    list_for_anarci_parsing = [( "scFv"+str(n), list_of_inputs[n][0] ) for n in range( len( list_of_inputs )) ]
    try:
        output = run_anarci( list_for_anarci_parsing, scheme='imgt', assign_germline = True, ncpu = ncores)
    except:
        logging.error("\tANARCI parsing failed to parse this chunk of sequences")
        output = []
    return output

class Main(object):
    """
    Class that initiates structural filtering analysis.
    --------------
    Parameters
        input_file - Txt file where antibody amino acid sequences are
                     in the column with Protein_Seq as a header
              ncpu - Number of CPU 
             chain - H (Heavy) at the moment
           species - human as default
    --------------
    Return
        output_name - Txt file that contains anarci parsed and numbered sequences
    """
    def __init__(self, input_file, 
                       ncpu = 4,
                       chain = "H", 
                       species = "human", 
                       output_name = "test_anarciparsed.txt", 
                       output_dir = "."):

       self.input_file = input_file
       if not os.path.isfile( self.input_file ):
           logging.warning("Input file does not exist: {0}".format( self.input_file ))
           raise AssertionError("Input file does not exist: ", self.input_file )

       self.ncpu = ncpu
       self.rate_of_analysis = 200000
       self.chain = chain
       self.species = species
       self.output_dir = output_dir
       self.output_name = output_name
       self.create_temp_csv()

    def create_temp_csv(self):
        "creating temp folder for an output"
	try:
            os.makedirs(self.output_dir)
	except:
            with open(join(self.output_dir, self.output_name), "wb") as csv_file:
                writer = csv.writer(csv_file, delimiter=',')

    def write_to_tempcsv(self, to_csv_list):
        "Function that writes temporary outputs to a csv file"
        rows = []
        for entry in to_csv_list:
            rows.append([ entry.input_seq, 
                          entry.redundancy,
                          entry.CDRH3,
                          "Unknown",   # Aboss flag1 ignoring
                          "Unknown",   # Aboss flag2 ignoring
                          entry.IG_V,
                          entry.IG_J,
                          entry.numbered, 
                          entry.seqID])

        with open( join( self.output_dir, self.output_name), "ab" ) as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            writer.writerows(rows)
               
    def find_seq_number(self, df):
        total = sum([len(chunk) for chunk in df])
        logging.info("\tTotal number of sequences: {0}".format( total ))
        return total

    def remove_stars(self, in_list):
        # Removing sequences if there is a * in antibody sequence
        return [x for x in in_list if "*" not in x[0] and len(x[0]) < 152]

    def _open_dataframe(self):
        self.df = pd.read_csv( self.input_file, 
                               sep = "\t",
                               iterator = True, 
                               chunksize = self.rate_of_analysis,
                               index_col = 0)

    def split_anarci_outputs(self, output, inputs):

        # Splitting lists to run structural viability in parallel
        anarci_output = np.array_split( output[1], self.ncpu )
        species_matrix = np.array_split( output[2], self.ncpu )
        input_sequences_list = np.array_split( inputs, self.ncpu )

        return anarci_output, species_matrix, input_sequences_list

    def get_list_to_save(self, correct_sequences_portion):
        # Collecting info from ANARCI parsing to save
        # This will be used later on to run SAAB+
        temp_write_to_csv_list = []
        for portion in range(len(correct_sequences_portion)):				
            temp_write_to_csv_list += correct_sequences_portion[portion][3]
        return temp_write_to_csv_list

    def run(self):
        check_anarci()
        self._open_dataframe()
	count = 0
	count_correct = 0
        
        df_data, df_count = tee( self.df )        
        num_of_entries = self.find_seq_number( df_count )

	for chunk in df_data:
            inputs = dict( chunk[[ "Protein_Seq", "Redundancy" ]].get_values() ).items()
            count += len( chunk )
            if not inputs:
                continue
            inputs = self.remove_stars( inputs )
            output = anarci_fun( inputs, self.ncpu )

            if not output:
                continue

            if len(output[1]) < 3000:
                self.ncpu = 1

            anarci_output, species_matrix, input_sequences_list = self.split_anarci_outputs( output, inputs )
           
            pool = mp.Pool( processes = self.ncpu )
            results = [ pool.apply_async(checking_sequence_viability, args = ( anarci_output[i], 
                                                                               species_matrix[i], 
                                                                               input_sequences_list[i], 
                                                                               self.chain, 
                                                                               self.species )) for i in xrange(self.ncpu) ]
            correct_sequences_portion = [ r.get() for r in results ]
            del results
            pool.close()
            inputs = []

            temp_write_to_csv_list = self.get_list_to_save( correct_sequences_portion )

            del correct_sequences_portion
            count_correct += len(temp_write_to_csv_list)

            progress = returnProgress( float(count)/num_of_entries, round((float(count_correct))/count,3) )
            # Write outputs to csv file
            self.write_to_tempcsv(temp_write_to_csv_list)
            
            # Saving outputs and printing status
            logging.info("\t{0}".format(progress))
        
        self.gzip_output_file()
    def gzip_output_file(self):
        "gzipping output_name"
        os.system("gzip -f " + join(self.output_dir, self.output_name))

def initArgParser():

    parser = ArgumentParser(description=__doc__)
    group = parser.add_argument_group('study arguments')
    group.add_argument('-f', action='store', dest='filename', default=False,
                        help="filename to be analyzed which is located in Datasets directory")
    group.add_argument('-s', action='store', dest='organism', default=False,
                        help="organism information in the study [human, mouse, rhesus, rat]")
    group.add_argument('-n', action='store', dest='ncores', default=16, type=int,
                        help="number of cores")
    group.add_argument("-c", action="store", dest="chain", default="H",
                        help="chain information (currently only H")
    group.add_argument("-o", action="store", dest="output_name", default="output_anarciparsed.txt",
                        help="output filename")
    group.add_argument("-d", action="store", dest="output_dir", default=".",
                        help="output filename")
    return parser

def printStatus(args):
    print "      FILE: {0}".format(args.filename)
    print "     CORES: {0}".format(args.ncores)
    print "  ORGANISM: {0}".format(args.organism)
    print "     CHAIN: {0}".format(args.chain)
    print "OutputName: {0}".format(args.output_name)
    print " OutputDir: {0}".format(args.output_dir)
    return vars(args)

if __name__ == "__main__":

    parser = initArgParser()
    args = parser.parse_args()
    if not args.filename:
        raise AssertionError("Study was not specified, exiting!!!")
    if not args.organism:
        raise AssertionError("Organism was not selected")
    if not args.output_name:
        args.output_name = "{0}_anarciparsed.txt".format(os.path.splitext(args.filename)[0])
    arguments = printStatus(args)

    STUDY = Main(args.filename, 
            ncpu = args.ncores,
            chain = args.chain, 
            species = args.organism,
            output_dir = args.output_dir,
            output_name = args.output_name)
    STUDY.run() 

