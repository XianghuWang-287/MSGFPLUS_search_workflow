from pyteomics import mzid,auxiliary
import csv
import os

# def is_decoy(psm, prefix="XXX_"):
#
#     # Get the list of protein accessions from the PSM
#     protein_accessions = psm.get('proteins', [])
#
#     # Check if all proteins have the decoy prefix
#     return all(protein.startswith(prefix) for protein in protein_accessions)
def is_decoy(psm, prefix=None):

    return all(pe['isDecoy'] for sii in psm['SpectrumIdentificationItem']
            for pe in sii['PeptideEvidenceRef'])
def get_msgf_evalue(psm):
    try:
        return float(psm['SpectrumIdentificationItem'][0]['MS-GF:EValue'])
    except (KeyError, ValueError):
        return float('inf')

input_mzid_file_path = '/home/user/LabData/XianghuData/MSGFPLUS_search_workflow/nf_output'

mzid_files = [os.path.join(input_mzid_file_path,f) for f in os.listdir(input_mzid_file_path) if f.lower().endswith('.mzid')]
# Filter PSMs based on the desired spectrum-level FDR threshold
filtered_psms = mzid.filter.chain.from_iterable(mzid_files,fdr = 0.01,key =get_msgf_evalue,decoy_prefix = "XXX_")


# Function to extract relevant information for each PSM
def extract_psm_info(psm):
    file_name = psm['name']
    scan = psm['scan number(s)']
    peptide = psm['SpectrumIdentificationItem'][0]['PeptideSequence']
    proteins = [pe['accession'] for sii in psm['SpectrumIdentificationItem']
            for pe in sii['PeptideEvidenceRef']]
    proteins = ','.join(proteins)
    charge = psm['SpectrumIdentificationItem'][0].get('chargeState', '')
    spec_prob = psm['SpectrumIdentificationItem'][0].get('MS-GF:SpecEValue', '')
    db_evalue = get_msgf_evalue(psm)
    Frag_method = psm['SpectrumIdentificationItem'][0].get('AssumedDissociationMethod', '')


    return {
        'MzIDFileName': file_name,
        'ScanNumber': scan,
        'FragMethod': Frag_method,
        'Peptide': peptide,
        'Protein': proteins,
        'Charge': charge,
        'SpecProb': spec_prob,
        'DB:EValue': db_evalue,
    }


# Path to the output filtered TSV file
output_filtered_tsv_file_path = '../data/PXD21518/PXD021518_filtered.tsv'

# Open the TSV file and write the header
with open(output_filtered_tsv_file_path, 'w', newline='') as tsv_file:
    fieldnames = [
        'MzIDFileName', 'ScanNumber', 'FragMethod',
        'Peptide', 'Protein', 'Charge', 'SpecProb', 'DB:EValue'
    ]
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    # Iterate over filtered PSMs and write the information to the TSV file
    for psm in filtered_psms:
        psm_info = extract_psm_info(psm)
        writer.writerow(psm_info)

print('Filtered results saved to:', output_filtered_tsv_file_path)

