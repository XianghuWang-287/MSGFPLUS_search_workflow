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

input_mzid_file_path = '/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_results/EK_Q_10_2.mzid'

psms = mzid.read(input_mzid_file_path)
# Filter PSMs based on the desired spectrum-level FDR threshold
filtered_psms = mzid.qvalues(psms,key =get_msgf_evalue,is_decoy= is_decoy,decoy_prefix = "XXX_",full_output=True)


# Function to extract relevant information for each PSM
def extract_psm_info(psm,score,qvalue):
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
    Qvalue=qvalue
    Score = score


    return {
        'MzIDFileName': file_name,
        'ScanNumber': scan,
        'FragMethod': Frag_method,
        'Peptide': peptide,
        'Protein': proteins,
        'Charge': charge,
        'SpecProb': spec_prob,
        'DB:EValue': db_evalue,
        'Qvalue': qvalue,
        'Score': score
    }


# Path to the output filtered TSV file
output_filtered_tsv_file_path = '../data/results/filtered.tsv'

# Open the TSV file and write the header
with open(output_filtered_tsv_file_path, 'w', newline='') as tsv_file:
    fieldnames = [
        'MzIDFileName', 'ScanNumber', 'FragMethod',
        'Peptide', 'Protein', 'Charge', 'SpecProb', 'DB:EValue', 'Qvalue','Score'
    ]
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    # Iterate over filtered PSMs and write the information to the TSV file
    for psm in filtered_psms:
        if psm['q']<=0.01:
            psm_info = extract_psm_info(psm['psm'],psm['score'],psm['q'])
            writer.writerow(psm_info)

print('Filtered results saved to:', output_filtered_tsv_file_path)

