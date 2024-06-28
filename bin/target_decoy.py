from urllib.request import urlretrieve
from pyopenms import *
from tqdm import tqdm
import pandas as pd
searchdb = "/home/user/LabData/XianghuData/MS_Cluster_datasets/DB_ID/UP000006548_3702.fasta"

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

target_decoy_database = "../data/search_td_mouse.fasta"
FASTAFile().store(target_decoy_database, targets + decoys)