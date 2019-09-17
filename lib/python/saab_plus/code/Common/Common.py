####Common functions & locaions##########
import random,os
from os.path import join
import json
#################
#Common functions
#################

#Generates a random filename
def random_filename():
	filename = str(int(1000000*random.random()))+str(int(1000000*random.random()))+str(int(1000000*random.random()))
	return filename

#################
#Common locations
#################
where_i_am = os.path.dirname(os.path.abspath(__file__))

#Where the numbered_datasets are stored
numbered_datasets_location = join(where_i_am,'../../data/numbered')

fread_db = join(where_i_am, "../../data/fread_db/db_CDRH3")

template_db = join(where_i_am, "../../data/structures/IMGT/")

#Where the structures of the PDBs are stored ( for FREAD comparisons.)
#structures_location = join(where_i_am,'../../data/structures/IMGT')

# pdb_length 
pdb_length_location = join(where_i_am, '../../data/pdb_length/pdb_lengths.csv')

numbered_germlines = join(where_i_am, '../../data/amino_sub/numbered_germline_genes.pkl')

canonical_germ = join(where_i_am, '../../data/amino_sub/germlines_can_annotated.txt')

# clustering file
clusters = join(where_i_am, '../../data/pdb_length/clusters.csv')
definitions = {		
		
		"chothia" : {"L1" : ["L24", "L25", "L26", "L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34"], 
				  "L2" : ["L50", "L51", "L52", "L53", "L54", "L55", "L56"],
				  "L3" : ["L89", "L90", "L91", "L92", "L93", "L94", "L95", "L96", "L97"], 
				  "H1" : ["H26", "H27", "H28", "H29", "H30", "H31", "H32"], 
				  "H2" : ["H52", "H53", "H54", "H55", "H56"] ,
				  "H3" : ["H95", "H96", "H97", "H98", "H99", "H100", "H101", "H102"]},
		"imgt" : {
				"H1" : ["H27", "H28", "H29", "H30", "H31", "H32", "H33", "H34", "H35", "H36","H37", "H38"], 
				"H2" : ["H56", "H57", "H58", "H59", "H60", "H61", "H62", "H63", "H64", "H65"] ,
				"H3" : ["H105", "H106", "H107", "H108", "H109", "H110", "H111", "H112", "H113", "H114", "H115", "H116", "H117"],
				"L1" : ["L27", "L28", "L29", "L30", "L31", "L32", "L33", "L34", "L35", "L36","L37", "L38"], 
				"L2" : ["L56", "L57", "L58", "L59", "L60", "L61", "L62", "L63", "L64", "L65"] ,
				"L3" : ["L105", "L106", "L107", "L108", "L109", "L110", "L111", "L112", "L113", "L114", "L115", "L116", "L117"],  }
		# need to introduce IMGT scheme here
		}


#Is this a CDR? If so tell which.
#e.g. res is "L52", do not supply insertion.
def is_CDR(res,deff='imgt'):

	if 'K' in res:
		res = res.replace('K','L')

	for CDR in definitions[deff]:
		if res in definitions[deff][CDR]:
			return CDR

	return False

####################
#Numbering by ANARCI
####################

#ANARCI-- number.
from anarci import run_anarci

#BUlk number sequences. 
#Input: dictionary from sequences ids to raw sequences
#Output: dictionary from sequence ids to numbered dictionaries of chothia ids to aa
def bulk_number(sequences):
	#print "[Common.py] Bulk Numbering ",len(sequences),'sequences...'
	numbered_sequences = {}
	for sequence_id in sequences:
		numbered_sequences[sequence_id] = number_sequence(sequences[sequence_id])
		
	return numbered_sequences

#From a mapping from chothia ids to amino acids, get the full sequence.
def get_sequence(s):
	sequence = ""
	for chid in sorted(s):
		
		sequence += s[chid][0]
	return sequence	

#Given a raw sequence, perform anarci annotation.
#Input: raw sequence
#output dictionary from chothia ids to amino acids and chothia cdr annotations


def number_sequence(query_seq):
	#Number the query sequence
    try:
	res = run_anarci([('q',query_seq)],scheme='imgt',assign_germline=True)
	query= res[1][0][0][0]
	    
	chain =  res[2][0][0]['chain_type']
	    
	#Translate the query into the format we want
	#[chid,insert,chain] = (aa,is_CDR)
	_query = {}
	    
	for sid in sorted(query):
		
		chid = (sid[0][0],sid[0][1].replace(' ',''),chain)
		if sid[1] == '-':
			continue
		_query[chid] = (sid[1],is_CDR(chain+str(sid[0][0])))
	
	query = _query

	return query,chain
    except:
        return None,None
    
###Strictly for testing purposes.	
if __name__ == '__main__':

#	print fetch_antigen_map()

#	quit()

	#Test bulk numbering
	seqs = {
				'1' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES',
				'2' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES',
				'3' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES',
				'4' :'DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES'
		   }
#	print bulk_number(seqs)
	
	
	#Test numbering
	print number_sequence('DIKMTQSPSSMYASLGERVTITCKASQDIRKYLNWYQQKPWKSPKTLIYYATSLADGVPSRFSGSGSGQDYSLTISSLES')

