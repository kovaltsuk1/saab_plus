from collections import defaultdict, namedtuple
from saab_plus.aboss_utils.region_definitions import Accept
from saab_plus.code.Common.Common import numbered_germlines
import cPickle as pickle
import json

list_of_frames ={
        "H":["cdrh3","fwh1", "cdrh1", "fwh2", "cdrh2", "fwh3","fwh4"],
        "L":["cdrl3","fwl1", "cdrl1", "fwl2", "cdrl2", "fwl3","fwl4"]}

good_aminos = sorted(list("QWERTYIPASDFGHKLCVNM"))

with open(numbered_germlines, "rb") as gg:
    numbered_germlines = pickle.load(gg)

output_tuple = namedtuple("output_tuple","input_seq redundancy CDRH3 IG_V IG_J numbered seqID")

CDR3 = {"H":"cdrh3","L":"cdrl3"}
FW1 = {"H":"fwh1","L":"fwl1"}
FW3 = {"H":"fwh3","L":"fwl3"}

a = Accept()
a.numbering_scheme = "imgt"
a.definition = "imgt"

def check_amino_viability(aa):
    if aa not in good_aminos and aa != "-":  # Checking here for unusual aa while skipping "-" entries at some position (e.g. cases of canonical CDRs)
        return False
    return True

def check_skip(aa):
    if aa == "-":
        return True
    return False

def rabbit_length_check(position_73, rabbit_position_84,
                                        regions, quality):
    if len(regions["fwh3"]) - position_73 + rabbit_position_84  < 38:   # Accounting for insertions at postion 73 and deletion at 84
            quality = ("Bad", "FW3 length is shorter IMGT defined"), rabbit_position_84
    if len(regions["fwh2"]) < 17:
            quality = ("Bad", "FW2 length is shorter IMGT defined")
    if len(regions["fwh4"]) < 10:
            quality = ("Bad", "FW4 length is shorter IMGT defined")
    if len(regions["cdrh3"]) > 37:
            quality = ("Bad", "CDR-H3 is chimeric")
    return quality

def standard_species_length_check(position_73, regions, quality):

    # We iterate heavy chain antibody regions, and check for their IMGT defined lengths
    for region, length in [("fwh2",17), ("fwh3",38),("fwh4",9), ("cdrh3", 37)]:
        if region == "fwh3":
            if len(regions[region]) - position_73 < length:
                quality = ("Bad", "{0} length is shorter IMGT defined".format(region))
                return quality
        elif region == "cdrh3":
            if len(regions[region]) > length:
                quality = ("Bad", "CDR-H3 is chimeric")
                return quality
        else:
            if len(regions[region]) < length:
                quality = ("Bad", "{0} length is shorter IMGT defined".format(region))
                return quality
    return quality

def light_chain_species_length_check(absent_positions,position_73,regions,quality):
        if len(regions["fwl3"]) + absent_positions - position_73 < 38: # Accounting for insertions at postion 73,81,83
                quality = ("Bad", "FW3 length is shorter IMGT defined")
        if len(regions["fwl2"]) < 17:
                quality = ("Bad", "FW2 length is shorter IMGT defined")
        if len(regions["fwl4"]) < 10:
                quality = ("Bad", "FW4 length is shorter IMGT defined")
        return quality

def standard_species_length_check_fw1(fw,all_region_aa,quality):
    """
    Here we check for indels in fw1. The problem is that for many Ig-seq datasets FW1 primers are used that start after position 15.
    """
    frame1 =  [int(x) for x in all_region_aa[fw]]
    if 10 not in frame1:
        if min(frame1) < 10:
            if max(frame1) - min(frame1) != len(frame1):
                quality = ("Bad", "FW1 length is shorter IMGT defined")
        else:
            if max(frame1) - min(frame1) + 1 != len(frame1):
                quality = ("Bad", "FW1 length is shorter IMGT defined")
    else:
        if max(frame1) - min(frame1) + 1 != len(frame1):
            quality = ("Bad", "FW1 length is shorter IMGT defined")
    return quality

def rabbit_length_check_fw1(fw,all_region_aa,rabbit_position_2,quality):
    """
    Framework 1 sometimes misses position 2 in rabbits.
    """
    frame1 =  [int(x) for x in all_region_aa[fw]]
    if 10 not in frame1:
        if min(frame1) == 1:
            if max(frame1) - min(frame1) - rabbit_position_2 != len(frame1):
                quality = ("Bad", "FW1 length is shorter IMGT defined")
        elif min(frame1) >= 3 and min(frame1) < 10:
            if max(frame1) - min(frame1) != len(frame1):
                quality = ("Bad", "FW1 length is shorter IMGT defined")
        else:
            if max(frame1) - min(frame1) + 1 != len(frame1):
                quality = ("Bad", "FW1 length is shorter IMGT defined")
    else:
        if max(frame1) - min(frame1) + 1 - rabbit_position_2 != len(frame1):
            quality = ("Bad", "FW1 length is shorter IMGT defined")
    return quality

def check_light_chain_missing_pos(residue, frame, position_meta):

    if frame == "fwl1":
        if str(residue) == "23":
            position_meta["conserved_position"] += 1
    if frame == "fwl3":
        if str(residue) == "104":
            position_meta["conserved_position"] += 1
        if str(residue) == "81" or str(residue) == "82":
            position_meta["absent_positions"] -= 1
        elif str(residue) == "73":
            position_meta["position_73"] += 1
    elif frame == "fwl4":
        if str(residue) == "118":
            position_meta["conserved_position"] += 1
def check_heavy_chain_missing_pos(residue, frame, position_meta,
                                                  species):
    if frame == "cdrh1":
        if str(residue) == "27":
            position_meta["conserved_position"] += 1
        elif str(residue) == "2" and species == "rabbit":
            position_meta["rabbit_position_2"] -=1
    if frame == "fwh3":
        if str(residue) == "104":
            position_meta["conserved_position"] += 1
        elif str(residue) == "73":
            position_meta["position_73"] += 1
        elif str(residue) == "84" and species == "rabbit":
            position_meta["rabbit_position_84"] -=1
    elif frame == "fwh4":
        if str(residue) == "118":
            position_meta["conserved_position"] += 1

def checking_structural_viability(numbering, species, chain, anarci_genes):
    """
    The main function of ANARCI parsing. Here we check for unusual indels in antibody regions across different species. We also check if the first numbered position of the amino acid sequence starts before postion 23
    -------------
    Parameters
            numbering - a list of ANARCI numbered antibody residues
    -------------
    Return:
            all_region_aa - ANARCI parsed numbered residues that comply with out knowledge of antibody folding
                  quality - Sequence quality. Good == No unusual indels, 3 CDRs and fw4 are covered.
                  regions - Separated amino acid sequence into antibody regions
    """
    all_region_aa = defaultdict(dict)
    regions = {}
    quality = "good"

    position_meta = {
            "conserved_position" : 0,
            "absent_positions": 2,
            "position_73": 0,
            "rabbit_position_2": 1,
            "rabbit_position_84": 1
                    }
    # Here we iterate over antibody regions
    # and check for structural viability
    for frame in list_of_frames[chain]:
        region = ''
        a.set_regions([frame])
        for n in range(len(numbering)):
            residue = list(numbering[n])
            if a.accept(residue[0], chain) != 1:
                continue
                    
            if residue[1] == "X":
                if residue[0][0] == 73:
                    residue[1] = "-"
                    continue

                if frame[0].lower() == "f":

                    # if X found in framework we replace it with germline amino acid
                    if frame[-1] == "4":
                        residue[1] = numbered_germlines[anarci_genes["IG_J"]][residue[0][0]]
                    else:
                        residue[1] = numbered_germlines[ anarci_genes["IG_V"]][residue[0][0]]
                else:
                    # if X found in CDR-H3 we return sequence as bad
                    if frame.lower() == "cdrh3":
                        quality = ("Bad", "X in CDR-H3 at position {0}".format(residue))
                        return None, quality, None
                    else:
                        # If X found in canonical CDR, we skip it.
                        residue[1] = "-"
                        continue
            else:
                if not check_amino_viability(residue[1]):
                    quality = ("Bad", "Weird amino acid")
                    return None, quality, None

            # Checking for missing amino acids
            if check_skip(residue[1]):
                continue
            # Here we check for indels outside canonical CDR
            if residue[0][1] != ' ':

                all_region_aa[frame][str(residue[0][0]) + residue[0][1]] = residue[1]
                if frame != CDR3[chain]:
                    quality = ("Bad", str(residue[0][0]) + residue[0][1])
                    return None,quality,None
            else:
                    all_region_aa[frame][str(residue[0][0])] = residue[1]

                    # Here we record the key positions that must be present in the antibody sequence
                    # We need to account for positions such as 73 that are often not present
                    # In light chain, positions 81 and 82 are often absent
                    # Position 118 tells us if numbering of CDR-H3 was successful
                    if chain == "L":
                        check_light_chain_missing_pos(residue[0][0], frame, position_meta)
                    elif chain == "H":
                        check_heavy_chain_missing_pos(residue[0][0], frame, position_meta,
                                                                            species)
            # Here we record amino acid sequences
            # of the antibody regions
            region += residue[1]
        regions[frame] = region

    # if we do not have the key conserved positions, we stop       
    if position_meta["conserved_position"] != 3:
            quality = ("Bad", "No Conserved Position")
            return all_region_aa, quality, regions

    # As rabbit has different antibody region lengths to other analysed species (e.g. human, mouse, rat, alpaca)
    # Also Heavy and Light chains have different region lengths
    if chain == "H" and species != "rabbit":
            quality = standard_species_length_check(position_meta["position_73"],
                                                    regions,
                                                    quality)
    elif chain == "H" and species == "rabbit":
            quality = rabbit_length_check(position_meta["position_73"],
                                          position_meta["rabbit_position_84"],
                                          regions,
                                          quality)
    elif chain == "L":
            quality = light_chain_species_length_check(position_meta["absent_positions"], 
                                                       position_meta["position_73"],
                                                       regions,
                                                       quality)
    if quality == "good":
        # Check for indels in framework 1 as some sequences do not start with position 1 
        fw = FW1[chain]
        if len(all_region_aa[fw]) == 0:
            return all_region_aa, quality, regions
        if species != "rabbit":
            quality = standard_species_length_check_fw1(fw, 
                                                        all_region_aa, 
                                                        quality)
        # Special case for rabbits
        else:
            quality = rabbit_length_check_fw1(fw,
                                              all_region_aa,
                                              position_meta["rabbit_position_2"], 
                                              quality)
    return all_region_aa, quality, regions

def check_depth(l):
    "anarci versions have different depths"
    if isinstance(l, list):
        return 1 + max(check_depth(item) for item in l)
    else:
        return 0

def check_V_and_J_gene_IMGT_alignment(species_matrix, gene, species):
    """
    Function that iterates through ANARCI species matrix to find sequence identities
    ----------
    Parameters
        species_matrix - ANARCI top species alignements with corresponding scores
                  gene - v_gene or j_gene to retrive the best alignment scores
    ----------
    Return
        IG_gene - IG gene, or None if no confident alignemnts
    """
    IG_gene = False

    depth = check_depth(species_matrix[0]["germlines"][gene])

    if depth == 1:
        if species_matrix[0]["germlines"][gene][0][0] != species:
            return IG_gene
        else:
            if species_matrix[0]["germlines"][gene][1] > 0.5:
                return species_matrix[0]["germlines"][gene][0][1] 

    if depth == 2:
        if species_matrix[0]["germlines"][gene][0][0][0] != species:
            range_genes = len(species_matrix[0]["germlines"][gene][0])
            for n in range(range_genes):
                if species_matrix[0]["germlines"][gene][0][n][0] == species:
                    if species_matrix[0]["germlines"][gene][1][n] > 0.5:
                        IG_gene = species_matrix[0]["germlines"][gene][0][n][1]
                    break
        else:
            if species_matrix[0]["germlines"][gene][1][0] > 0.50:
                IG_gene = species_matrix[0]["germlines"][gene][0][0][1]
    return IG_gene

def check_IMGT_alignment(species_matrix, species, chain):
    """
    Checking that antibody sequences align to IGHV and IGHJ genes
    ------------
    Parameters
        species_matrix - ANARCI top species alignements with corresponding scores
    ------------
    Return
                  IG_V - IGHV gene
                  IG_J - IGHJ gene
        IGH_IGJ_status - None if IGHV or IGHJ do not align to HMMs
    """
    # Here we aling sequence to HMM 
    IGH_IGJ_status = "unknown"
    IG_V = check_V_and_J_gene_IMGT_alignment(species_matrix,"v_gene", species)
    IG_J = check_V_and_J_gene_IMGT_alignment(species_matrix,"j_gene", species)
    if chain == "H":
        if not IG_J or not IG_V:
            IGH_IGJ_status = None
        elif IG_V == False or IG_J == False:
            IGH_IGJ_status = None
        elif IG_V[:3] != "IGH" or IG_J[:3] != "IGH":
            IGH_IGJ_status = None
    elif chain == "L":
        if not IG_J or not IG_V:
            IGH_IGJ_status = None
        elif IG_V[:3] == "IGH" or IG_J[:3] == "IGH":
            IGH_IGJ_status = None
        elif IG_V == False or IG_J == False:
            IGH_IGJ_status = None
    return IG_V, IG_J, IGH_IGJ_status

def viable_outputs(anarci_output):
    "some batches can be empty"
    try:
        if anarci_output == None or anarci_output[0][0] == []:
            return False
        return True
    except:
        if anarci_output.all() == None:
            return False
        return True

def checking_sequence_viability(anarci_output, species_matrix, input_sequences_list, chain, species):
    """
    Function that initiates checking for indels, HMM alignmets,chimeric sequences
    -------------
    Parameters
               anarci_output - ANARCI numbering outputs
              species_matrix - ANARCI top species alignements with corresponding scores
        input_sequences_list - Original full length antibody sequences
    -------------
    Return
                 anarci_parsed - ANARCI parsed sequences
                   V_gene_dict - Track of IG_V gene
                   J_gene_dict - Track of IG_J gene
        anarci_parsing_outputs - List of (Antibody sequence, Redundancy, IG_V, IG_J, Numbered residue/positions, Sequence ID)
    """
    anarci_parsed = {}
    V_gene_dict, J_gene_dict = {}, {}
    anarci_parsing_outputs = []
    #Checking if have parsing outputs
    if anarci_output == []:
        return anarci_parsed, V_gene_dict, J_gene_dict, anarci_parsing_outputs
    for u in range(len(anarci_output)):
        if not viable_outputs(anarci_output[u]):
            continue
        else:
            # Here we check for the correct antibody chain 
            IG_V, IG_J, IGH_IGJ_status = check_IMGT_alignment(species_matrix[u],
                                                                species, chain)
            if IGH_IGJ_status == None:
                continue

            aa, status, regions = checking_structural_viability(anarci_output[u][0][0],
                                                                species, chain, {"IG_V": IG_V,
                                                                                 "IG_J": IG_J})
            # if antibody sequence fails anarci
            if status != "good":
                continue
            # Some sequences have very long fw4 regions or duplicated fw4 regions
            if len(input_sequences_list[u][0][input_sequences_list[u][0].rfind(regions[list_of_frames[chain][-1]]):]) > 13 or input_sequences_list[u][0].count(regions[list_of_frames[chain][-1]]) > 1:
                continue
            anarci_parsed[input_sequences_list[u][0]] = json.dumps(aa)
            V_gene_dict[IG_V] = 1
            J_gene_dict[IG_J] = 1

            # Saving paring Outputs
            anarci_parsing_outputs.append( output_tuple( input_seq = input_sequences_list[u][0],
                                                        redundancy = input_sequences_list[u][1],
                                                             CDRH3 = regions[CDR3[chain]],
                                                              IG_V = IG_V,
                                                              IG_J = IG_J,
                                                          numbered = json.dumps(aa),
                                                             seqID = "Missing")
                                         )
    return anarci_parsed, V_gene_dict, J_gene_dict, anarci_parsing_outputs

if __name__ == "__main__":
    pass
