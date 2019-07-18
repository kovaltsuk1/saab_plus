from FREAD.pyfread_api import run_fread
import os
from ..Common.Common import fread_db, template_db
fread_cache = {}


#For fread, one residue BEFORE the start of the loop.
loop_starts = {'L1':26,'L2':55,'L3':104,'H1':26,'H2':55,'H3':104}

#Get the CDR sequences from the numbered object.
#Get the CDR sequences from the numbered object.
def extract_cdrs(numbered):
    """
    Function to extract cdrh3 sequence from
    anarci numbering
    """
    cdrs = {}
    pos_112 = ""

    for chid in sorted(numbered):
        cdr_code = numbered[chid][1]
        if cdr_code != False:
            #means we are a CDR
            aa = numbered[chid][0]
            if cdr_code not in cdrs:
                cdrs[cdr_code] = ""
            if cdr_code != "H3":
                cdrs[cdr_code]+=aa
            else:
                if chid[0] < 112:
                    cdrs[cdr_code]+=aa
                elif chid[0]== 112:
                    pos_112 +=aa
                elif chid[0] > 112:
                    cdrs[cdr_code]+=pos_112[::-1]
                    pos_112 = ""
                    cdrs[cdr_code]+=aa
    return cdrs				

#Perform fread search for the given sequence.
#loop -- H1, H2 etc.
#template_pdb -- e.g. 12e8
#templaet_chain -- e.g. P
#sequence -- e.g. 'DPEIGD'
#return_top -- how many top results to fetch?
#Result: json-formatted best match {'seq': u'DPEIGD', 'str': u'12e8P', 'scr': 45} or None
def perform_loop_alignment(loop, template_pdb, template_chain,
                                      sequence, return_top=10):

    cache_key = (loop,template_pdb+template_chain,sequence)
    if cache_key in fread_cache:
            return fread_cache[cache_key]

    #Template location.
    template = os.path.join(template_db,"{0}{1}_no_cdrs.pdb".format(template_pdb, template_chain))

    results = run_fread(fread_db, template, loop_starts[loop], sequence, template_chain,'')

    #Get the best matches.
    matches = []
    
    for decoy in results:
            #Debugging...
            #print decoy.struc, decoy.startres, decoy.startinscode, decoy.length, decoy.score,decoy.seq
            matches.append((decoy.score,{'str':decoy.struc,'seq':decoy.seq,'scr':decoy.score,'qu':sequence, "ca":decoy.anchor_rmsd_open}))
    
    fread_results = sorted(matches,reverse=True)[0:min(return_top,len(matches))]
    #Cache the resulst
    fread_cache[cache_key] = fread_results

    return fread_results

#Strictly for debugging constituent funcionality.
if __name__ == '__main__':

	import sys
	cmd = sys.argv[1]
	if cmd == 'fetch_loop':
		loop = 'H2'
		sequence = 'DPEIGD'
		print perform_loop_alignment(loop,'12e8','P',sequence)
	if cmd == 'fetch_cdrs':
		import pickle
		query = pickle.load(open('sample_sequence.pckl'))
		print extract_cdrs(query[0])
	
