import pandas as pd
import os
import numpy as np
from Common.Common import pdb_length_location, clusters
from ast import literal_eval


class ApplyESS:
    chain_converter = {"H": 
                         {"CDR3":"H3",
                          "CDR1":"H1",
                          "CDR2":"H2"},
                       "L":
                         {"CDR3":"L3",
                          "CDR1":"L1",
                          "CDR2":"L2"}
                         }

    def __init__(self, filename, chain ):
        """
        Applying ESS score thresholds to SAAB+ outputs
        -----------
        Parameters:
               filename - Name of the file that containg SAAB+ 
                          annotated sequences
                  chain - Chain information for BCR repertoire  
        self.pdb_length - Dictionary of pdb_templates and
                          their corresponding CDR-H3 lengths
        """
        self.filename = filename
        self.chain = self.chain_converter[chain]["CDR3"]
        self.CDR1 = self.chain_converter[chain]["CDR1"]
        self.CDR2 = self.chain_converter[chain]["CDR2"]
        self.CDR3_template = "CDR-{0}_template".format(self.chain)
        self.CDR3_sequence = "CDR-{0}_sequence".format(self.chain) 
        self.pdb_length = self.load_pdb_length()
        if not os.path.isfile( self.filename ):
            raise AssertionError("File does not exist: ", self.filename)
        self.df = pd.read_csv( self.filename, sep = "\t", header = None,
                                                          names=[ "Protein_Seq", self.CDR3_template, 
                                                                  "Canonical_classes", "Redundancy",
                                                                  "Framework_template", self.CDR3_sequence , 
                                                                   "ESS"])
        # Vectorized check_ess for speed
        self.check_ess_vf = np.vectorize( self.check_ess )

    def load_pdb_length(self):
        # Loading a mapping dictionary {"PDB: "amino acid length"}
        df = pd.read_csv( pdb_length_location[self.chain], 
                          sep = "\t", 
                          index_col = 0)
        return dict(df[["pdb", "Length"]].get_values())

    def check_ess(self, essScore, loop_length):
        """
        Function that checks whether the ESS score is above 
        the cut-off for a given CDR-H3 length
        """
        if 4 < loop_length < 13 and essScore >= 25:
            return True
        elif 13 <= loop_length < 15 and essScore >= 45:
            return True
        elif 15 <= loop_length < 17 and essScore >= 55:
            return True
        return False

    @property
    def adding_pdb_length(self):
        self.df["PDB_length"] = self.df[self.CDR3_template].apply( self.add_pdb_length )
        self.df["CDRH3Length"] = self.df[self.CDR3_sequence].apply( len )
        self.df["PDB_length"] = np.where( self.df["PDB_length"] == 0, self.df["CDRH3Length"],
                                                                      self.df["PDB_length"] )
        self.df.drop( columns = ["CDRH3Length"], inplace = True)
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
        self.df['Annotation'] = self.check_ess_vf( self.df["ESS"],
                                                   self.df["PDB_length"]) 
        self.df.drop( columns = [ "PDB_length" ], inplace = True )
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
        df_clusters = pd.read_csv( clusters[self.chain], 
                                   sep = "\t", 
                                   index_col = 0).reset_index().rename(
                                                                        columns = { "clusters": "Clusters" })
        self.df = self.df.merge( right = df_clusters, 
                                 left_on = self.CDR3_template, 
                                 right_on = "templates", 
                                 how = "left").drop( columns = [ "templates" ] ) 

        self.df["Clusters"].fillna( "None", inplace = True) # CDR-H3 below the ESS will have None as Cluster name
        return self

    def expanding_canonical_class_dict(self):
        """
        Reformating the canonical class column.
        Instead of Canonical_classes, we will have 
        two new columns: H1_Canonical_Class
                         H2_Canonical_Class
        """
        self.df = pd.concat( [ self.df, self.df.apply(lambda x: pd.Series([ literal_eval(x["Canonical_classes"])[self.CDR1], 
                                                                            literal_eval(x["Canonical_classes"])[self.CDR2] ], 
                                                                            index = [ "{0}_Canonical_Class".format(self.CDR1), 
                                                                                      "{0}_Canonical_Class".format(self.CDR2) ] ) , axis = 1 )],
                                                                                                                 axis  = 1 )
        self.df.drop( columns = [ "Canonical_classes" ], inplace = True) # Cleaning

    def save(self):
        """
        Saving output file as a gzip file to save space
        """
        self.df.to_csv( self.filename, sep = "\t", compression = "gzip" )

if __name__ == "__main__":

    apply_ess_instance = ApplyESS("../output_structure.txt.gz")
    apply_ess_instance.adding_pdb_length()
    apply_ess_instance.apply_ess()
    apply_ess_instance.clustering()
    apply_ess_instance.expanding_canonical_class_dict()
    apply_ess_instance.save()

