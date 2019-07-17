import pandas as pd
import cPickle as pickle
import os
import numpy as np
from Common.Common import pdb_length_location

class ApplyESS:
    
    def __init__(self, filename):
        self.filename = filename
        self.pdb_length = self.load_pdb_length()
        if not os.path.isfile(self.filename):
            raise AssertionError("File does not exist: ", self.filename)
        self.df = pd.read_csv(self.filename, sep="\t", header=None,
                                    names=["Protein_Seq", "H3pdb", "Canon", "Redundancy",
                                              "FrameWork", "CDRHSeq", "ESS"])
    def load_pdb_length(self):
        with open(pdb_length_location, "rb") as gg:
            data = pickle.load(gg)
            return dict(data[["pdb", "Length"]].get_values())

    def adding_pdb_length(self):
        self.df["PDB_length"] = self.df["H3pdb"].apply(self.add_pdb_length)
        self.df["CDRH3Length"] = self.df["CDRHSeq"].apply(len)
        self.df["PDB_length"] = np.where(self.df["PDB_length"] == 0, self.df["CDRH3Length"], self.df["PDB_length"])
        self.df.drop(columns=["CDRH3Length"], inplace=True)

    def add_pdb_length(self, line):
        # Checking CDR-H3 length of PDB templates
        if line == "None":
            return 0
        return self.pdb_length[line]
       
    def check_ess(self, essScore, loop_length):
        """Function that checks if ESS is score is above 
        the cut-off for a given CDR-H3 length
        """
        if loop_length <13 and essScore >= 25:
            return True
        elif 13 <= loop_length <15 and essScore >= 45:
            return True
        elif 15 <= loop_length <17 and essScore >= 55:
            return True
        return False
    
    def apply_ess(self):
        self.df['StructurallyAnnotated'] = self.df.apply(lambda x: self.check_ess(x['ESS'], x['PDB_length']), axis=1)
        self.df.drop(columns=["PDB_length"], inplace=True)
    
    def save(self):
        self.df.to_csv(self.filename, sep="\t", compression="gzip")

if __name__ == "__main__":

    apply_ess_instance = ApplyESS("../output_structure.txt.gz")
    apply_ess_instance.adding_pdb_length()
    apply_ess_instance.apply_ess()
    apply_ess_instance.save()

