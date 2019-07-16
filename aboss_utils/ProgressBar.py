import sys 
import numpy as np

def returnProgress(proportion, anarci_parsed=None):

    barLength = 100
    if not isinstance(proportion, float ):
        return None

    if proportion >= 1:
        proportion = 1
    block = int(round(barLength*proportion))

    if not anarci_parsed:
        return "\rProgress: [{0}] {1}%".format( "#"*block + "-"*(barLength-block), round(proportion*100,2))

    return "\rProgress: [{0}] {1}%\nPercent Of ANARCI Parsed Sequences: {2}%".format( "#"*block + "-"*(barLength-block), 
                                                                                                round(proportion*100,2), 
                                                                                                anarci_parsed*100)
if __name__ == "__main__":
    for x in np.linspace(0,1,11):
        returnProgress(x,0.7)
