from ..Common.Common import numbered_datasets_location
from os.path import join
from os import listdir
import pickle
import logging

#This can be loaded in whole.
def structural_reference():	
    source = join(numbered_datasets_location,'sabdab',"IMGT")
    all_strucs = {}
    for chunk in listdir(source):
        logging.info("\tLoading Structures chunk: {0}".format(chunk))
        d = pickle.load(open(join(source,chunk)))
        logging.info("################")
        for pdb in d:
            all_strucs[pdb] = d[pdb]
            
    return all_strucs

if __name__ == '__main__':
	strucs = fetch_sabdab_seqs()
	print "Got",len(strucs)

