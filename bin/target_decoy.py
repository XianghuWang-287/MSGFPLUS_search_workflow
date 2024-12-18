from urllib.request import urlretrieve
from pyopenms import *
from tqdm import tqdm
import pandas as pd
searchdb = "/home/user/LabData/XianghuData/MALDI_Data/Human_uniprotkb_proteome_UP000005640_2024_09_21.fasta"

targets = list()
decoys = list()
FASTAFile().load(searchdb, targets) # read FASTA file into a list of FASTAEntrys
decoy_generator = DecoyGenerator()
for entry in tqdm(targets):
    rev_entry = FASTAEntry(entry) # copy entry
    rev_entry.identifier = "XXX_" + rev_entry.identifier # mark as decoy
    aas = AASequence().fromString(rev_entry.sequence) # convert string into amino acid sequence
    rev_entry.sequence = decoy_generator.reverseProtein(aas).toString() # reverse
    decoys.append(rev_entry)

target_decoy_database = "../data/Human_uniprotkb_proteome_td.fasta"
FASTAFile().store(target_decoy_database, targets + decoys)