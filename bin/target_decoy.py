from urllib.request import urlretrieve
from pyopenms import *
import pandas as pd
searchdb = "./target-decoy-database-main/target_decoy/src/human_swiss_prot_target.fasta"

targets = list()
decoys = list()
FASTAFile().load(searchdb, targets) # read FASTA file into a list of FASTAEntrys
decoy_generator = DecoyGenerator()
for entry in targets:
    rev_entry = FASTAEntry(entry) # copy entry
    rev_entry.identifier = "XXX_" + rev_entry.identifier # mark as decoy
    aas = AASequence().fromString(rev_entry.sequence) # convert string into amino acid sequence
    rev_entry.sequence = decoy_generator.reverseProtein(aas).toString() # reverse
    decoys.append(rev_entry)

target_decoy_database = "../data/search_td.fasta"
FASTAFile().store(target_decoy_database, targets + decoys)