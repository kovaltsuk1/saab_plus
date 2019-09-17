import pandas as pd
from saab_plus.code.Common.Common import canonical_germ

class AddGermCan:

    def __init__(self, saab_output, input_dataframe):
        self.saab_output = saab_output
        self.input_dataframe = input_dataframe
        self.open_pandas()

    def open_pandas(self):
        self.saab_plus = pd.read_csv(self.saab_output, sep = "\t", 
                                                       index_col = 0)
        self.input_dataframe = pd.read_csv(self.input_dataframe, sep = "\t", 
                                                              index_col = 0, header=0).drop(columns=["Redundancy"])
        return self

    def merge_V_Gene(self):
        self.saab_plus = self.saab_plus.merge( self.input_dataframe, left_on = "Protein_Seq", 
                                                                    right_on = "Protein_Seq")
        return self

    def merge_germ_can(self):
        df = pd.read_csv( canonical_germ, sep="\t",
                                          header = 0, 
                                          index_col = 0).drop(columns = ["Productive", 
                                                                          "Sequence"]).rename(columns={"H1":"H1_germ", 
                                                                                                        "H2":"H2_germ"})
        self.saab_plus = self.saab_plus.merge(df, left_on = "V_GENE",
                                                  right_on = "V_GENE")
        return self
    def save(self):
        self.saab_plus.to_csv(self.saab_output, 
                                    sep = "\t", 
                                    compression = "gzip")
