import pandas as pd
import os
import numpy as np
from Common.Common import pdb_length_location, clusters

class ApplyESS:
    
    def __init__(self, filename):
        """
        Applying ESS score thresholds to SAAB+ outputs
        -----------
        Parameters:
            filename - name of the file that containg SAAB+ 
                       annotated sequences
            self.pdb_length - Dictionary of pdb_templates and
                              their corresponding CDR-H3 lengths
        """
        self.filename = filename
        self.pdb_length = self.load_pdb_length()
        if not os.path.isfile(self.filename):
            raise AssertionError("File does not exist: ", self.filename)
        self.df = pd.read_csv(self.filename, sep="\t", header=None,
                                    names=["Protein_Seq", "CDR-H3_template", 
                                           "Canonical_classes", "Redundancy",
                                           "Framework_template", "CDR-H3_sequence", 
                                           "ESS"])
        # Vectorized check_ess for speed
        self.check_ess_vf = np.vectorize(self.check_ess)

    def load_pdb_length(self):

        df = pd.read_csv(pdb_length_location, sep="\t", index_col=0)
        return dict(df[["pdb", "Length"]].get_values())

    def check_ess(self, essScore, loop_length):
        """
        Function that checks if ESS is score is above 
        the cut-off for a given CDR-H3 length
        """
        if 4 < loop_length <13 and essScore >= 25:
            return True
        elif 13 <= loop_length <15 and essScore >= 45:
            return True
        elif 15 <= loop_length <17 and essScore >= 55:
            return True
        return False

    @property
    def adding_pdb_length(self):
        self.df["PDB_length"] = self.df["CDR-H3_template"].apply(self.add_pdb_length)
        self.df["CDRH3Length"] = self.df["CDR-H3_sequence"].apply(len)
        self.df["PDB_length"] = np.where(self.df["PDB_length"] == 0, self.df["CDRH3Length"], self.df["PDB_length"])
        self.df.drop(columns=["CDRH3Length"], inplace=True)
        return self

    def add_pdb_length(self, line):
        # Checking CDR-H3 length of PDB templates
        if line == "None":
            return 0
        return self.pdb_length[line]
       
    def apply_ess(self):
        """
        Applying quality control to our CDR-H3 sequence structural annotation.
        Each CDR-H3 loop length has defined quality controls.
        If good quality CDR-H3 template was not found, we do not consider this
        prediction in our analysis
        ----------------
        Parameters:

            self.df["ESS"] - Value of your homology template match. Higher values
                             indicate better matches. Loops of different lengths
                             have defined ESS quality thresholds.
            self.df["PDB_Length"] - Length of CDR-H3 loop
        ----------------
        Returns:

            self.df["Annotation"] - Boolean. It indicates whether a good quality
                                    CDR-H3 template was found
        """
        self.df['Annotation'] = self.check_ess_vf(self.df["ESS"],
                                                 self.df["PDB_length"]) 
        self.df.drop(columns=["PDB_length"], 
                               inplace=True)
        return self

    def clustering(self):
        """
        Mapping CDR-H3 templates to CDR-H3 clusters
        CDR-H3 cluster represents most connected CDR-H3 template
        amongst clustered templates within 0.6A
        --------
        Return:

            self.df["Clusters"] - Clusters column now contains information
                                  of clustered CDR-H3 templates.
        """
        df_clusters = pd.read_csv(clusters, sep="\t", index_col=0).reset_index().rename(columns={"clusters":"Clusters"})
        self.df = self.df.merge(right=df_clusters, left_on="CDR-H3_template", 
                                                  right_on="templates", how="left").drop(columns=["templates"]) 
        return self

    def save(self):
        """
        Saving output file as a gzip file to save space
        """
        self.df.to_csv(self.filename, sep="\t", 
                                      compression="gzip")

if __name__ == "__main__":

    apply_ess_instance = ApplyESS("../output_structure.txt.gz")
    apply_ess_instance.adding_pdb_length()
    apply_ess_instance.apply_ess()
    apply_ess_instance.clustering()
    apply_ess_instance.save()

