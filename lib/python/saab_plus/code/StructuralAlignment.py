from Common.Common import get_sequence, numbered_datasets_location, number_sequence
from Alignment.Align import align_sequences
from Alignment.LoopAlignment import perform_loop_alignment, extract_cdrs
from DataManagement.SAbDab import structural_reference
import multiprocessing as mp
import numpy as np
from scalop.inhouse_predict import _assign
from string import ascii_letters
import logging
from collections import namedtuple

output_tuple = namedtuple("output_tuple",["CDR_H3_template",
                                         "Canonical_classes",
                                         "Redundancy",
                                         "Framework_template",
                                         "CDR_H3_sequence",
                                         "ESS"])

formatDict= { 
              "Light": {"CDR3": "L3",
                        "CDR1": "L1",
                        "CDR2": "L2"},
              "Heavy": {"CDR3":"H3",
                       "CDR1": "H1",
                       "CDR2":"H2" } 
              }

def parsecdr( numcdrseq ):
    "parsing CDRs for SCALOP"
    bigcdrseq = []
    for pos, res in sorted( numcdrseq ):
        if pos[-1] in ascii_letters:
            _pos = (int(pos[:-1]), pos[-1])
        else:
            _pos = (int(pos), ' ')
            bigcdrseq.append([_pos,res])
    return bigcdrseq

def get_best_match(query,structures,region=None):
    """
    Finding the framework with the highest sequence identity
    """
    curr_best_sid = 0
    curr_best_pdb = 0		
    j = 0
    for struc in structures:  
        if query[1] == structures[struc][1]:
            j+=1
            """
            Here we get sids by sequence identity?	
            """
            sid = align_sequences( query[0], 
                                   structures[struc][0], 
                                   region = region )
            if sid> curr_best_sid:
                curr_best_sid = sid
                curr_best_pdb = struc
    results = {'best_sid':curr_best_sid,
               'best_pdb':curr_best_pdb}
    return results

def findmaxESSminCA(results):
    "sort predictions by ESS and Ca-Ca values"
    maxESS = max([x[0] for x in results])
    return np.argmin([x[1]["ca"] for x in results if x[0] == maxESS])

#Given anchors, employ FREAD to get the best CDR-templates.
def get_best_cdr_match( cdrs , fread_template, chain):
    """Running FREAD and sorting outputs
    Inputs:
        CDRH3 sequence in cdr["H3"]
        Framework pdb
    """
    pdb_template = fread_template[0:4]
    pdb_chain = fread_template[4]

    if not cdrs.get("H3", None):
        return {}, 0

    fread_results = perform_loop_alignment( "H3", 
                                            pdb_template,
                                            pdb_chain, 
                                            cdrs["H3"])
    if not fread_results:
        return {}, 0
    maxESSminCA = findmaxESSminCA( fread_results )
    essScore = fread_results[maxESSminCA][1]["scr"]
    results= { formatDict[chain]["CDR3"]: fread_results[maxESSminCA][1]["str"] }

    return results, essScore

def fetch_canonicals(query):
    """
    Function to structurally annotate canonical 
    classes.
    Passing numbered canonical class sequences
    """
    can = {}
    for canonical in ["H1", "H2"]:
        if query.sequencemeta.get( canonical, None ):
            try:
                _,_,_can,_ = _assign(parsecdr(query.sequencemeta[canonical].items()), 
                                     canonical,
                                     'imgt', 
                                     'imgt', 
                                     'latest')
                can[canonical] = _can
            except:
                can[canonical] = "None"
        else:
            can[canonical] = "None"
    return can

def init_fread(CDRSequences, full_results, chain):
    """
    Initiating FREAD template prediction
    -----------
    Parameters
        CDRSequences - Dictionary that contains CDR-H3 loop sequence
                       e.g. {H3:"ARFDY"}
        full_results - Dictionary that has the pdb code for the 
                       framework match. It will be used in loop grafting
               chain - Antibody chain information (e.g. H)
    ------------
    Returns
        fread_results - Best predicted PDB template
             essScore - ESS score for the selected FREAD template
    """
    try:
        fread_results, essScore = get_best_cdr_match( CDRSequences, 
                                                      full_results['best_pdb'], 
                                                      chain )
    except:
        essScore = 0
        fread_results = False
        logging.error("\tFREAD could not run: {0}".format( query.sequence ))
    return fread_results, essScore

def align_single_sequence(queries, structures, chain):
    """
    Function that annotates anarci numbered sequences
    with structural information.
    """
    output_dict = {}
    for query in queries:
        fread_results = False

        #Get best framework pdb_template
        full_results = get_best_match( query.numbering, 
                                       structures)
        CDRSequences = extract_cdrs( query.numbering[0] )
        CDR3Sequence = CDRSequences.get( "H3", None )

        # SCALOP
        can = fetch_canonicals( query ) 

        if "best_pdb" not in full_results:
            # if we cannot find framework match
            output_dict[ query.sequence ] = output_tuple( CDR_H3_template = "None",
                                                          Canonical_classes = can,
                                                          Redundancy = query.sequencemeta["Redundancy"],
                                                          Framework_template = "None", 
                                                          CDR_H3_sequence = CDR3Sequence, 
                                                          ESS = 0 )
            continue

        # FREAD
        if CDR3Sequence:
            fread_results, essScore = init_fread( CDRSequences, 
                                                  full_results, 
                                                  chain ) 
        else:
            essScore = 0
            CDR3Sequence = ""

        # Recording outputs
        try:
            if fread_results:
                output_dict[query.sequence] = output_tuple( CDR_H3_template = fread_results[formatDict[chain]["CDR3"]],
                                                            Canonical_classes = can,
                                                            Redundancy = query.sequencemeta["Redundancy"], 
                                                            Framework_template = full_results['best_pdb'],
                                                            CDR_H3_sequence = CDR3Sequence, 
                                                            ESS = essScore )
            else:
                output_dict[query.sequence] = output_tuple( CDR_H3_template = "None",
                                                            Canonical_classes = can,
                                                            Redundancy = query.sequencemeta["Redundancy"], 
                                                            Framework_template = full_results['best_pdb'],
                                                            CDR_H3_sequence = CDR3Sequence, 
                                                            ESS = essScore )
        except IndexError:
            output_dict[query.sequence] = output_tuple( CDR_H3_template = "None",
                                                        Canonical_classes = can,
                                                        Redundancy = query.sequencemeta["Redundancy"], 
                                                        Framework_template = full_results['best_pdb'],
                                                        CDR_H3_sequence = CDR3Sequence, 
                                                        ESS = essScore )
    return (output_dict, "_")

if __name__ == '__main__':
    pass
