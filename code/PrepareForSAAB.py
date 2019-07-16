import cPickle as pickle
import csv
import sys, os,json
import shutil
import multiprocessing as mp

def oas_output_parser(numbered):
    saab_input = dict()
    for region in numbered:
        if "cdr" not in region:
            CDR= False
            if region[2] == "h":
                    chain = "H"
            else:
                    chain = "L"
        else:
            if region[3] == "h":
                    chain = "H"
            else:
                    chain = "L"
            CDR =chain+region[-1]
        for position in numbered[region]:
            try:
                saab_input[int(position),"",chain] = (numbered[region][position], CDR)
            except:
                saab_input[int(position[:-1]),position[-1:],chain] = (numbered[region][position], CDR)
    return (saab_input,chain)

if __name__ =="__main__":
    pass

