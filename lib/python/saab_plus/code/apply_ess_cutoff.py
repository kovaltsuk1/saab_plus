import pandas as pd
import cPickle as pickle
import os
import numpy as np
from Common.Common import pdb_length_location


class ApplyESS:
    
    def __init__(self, filename):
        """
        Applying ESS score thresholds to SAAB+ outputs
        """
        self.filename = filename
        self.pdb_length = self.load_pdb_length()
        if not os.path.isfile(self.filename):
            raise AssertionError("File does not exist: ", self.filename)
        self.df = pd.read_csv(self.filename, sep="\t", header=None,
                                    names=["Protein_Seq", "CDR-H3_template", "Canonical_classes", "Redundancy",
                                              "Framework_template", "CDR-H3_sequence", "ESS"])
        # Vectorized check_ess for speed
        self.check_ess_vf = np.vectorize(self.check_ess)

    def load_pdb_length(self):
        with open(pdb_length_location, "rb") as gg:
            data = pickle.load(gg)
            return dict(data[["pdb", "Length"]].get_values())

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

    def adding_pdb_length(self):
        self.df["PDB_length"] = self.df["CDR-H3_template"].apply(self.add_pdb_length)
        self.df["CDRH3Length"] = self.df["CDR-H3_sequence"].apply(len)
        self.df["PDB_length"] = np.where(self.df["PDB_length"] == 0, self.df["CDRH3Length"], self.df["PDB_length"])
        self.df.drop(columns=["CDRH3Length"], inplace=True)

    def add_pdb_length(self, line):
        # Checking CDR-H3 length of PDB templates
        if line == "None":
            return 0
        return self.pdb_length[line]
       
    def apply_ess(self):
        self.df['Annotation'] = self.check_ess_vf(self.df["ESS"],
                                                 self.df["PDB_length"]) 
        self.df.drop(columns=["PDB_length"], 
                               inplace=True)
    def save(self):
        self.df.to_csv(self.filename, sep="\t", compression="gzip")

if __name__ == "__main__":

    apply_ess_instance = ApplyESS("../output_structure.txt.gz")
    apply_ess_instance.adding_pdb_length()
    apply_ess_instance.apply_ess()
    apply_ess_instance.save()

