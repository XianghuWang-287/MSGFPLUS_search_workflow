import matplotlib.pyplot as plt
from pyteomics import mass
import pandas as pd
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
# Function to read the spectrum data from a CSV file
def read_spectrum(csv_file):
    data = pd.read_csv(csv_file)
    return data['mz'].values, data['intensity'].values

def find_closest_mass(mass_list, target_mass, tolerance):
    for mass in mass_list:
        if abs(mass - target_mass) <= tolerance:
            return mass
    return None

# Read the spectra for peptides A and B
mz1, intensities1 = read_spectrum('../data/Peptide_spec_demo/EK_Q_07.mzML_scan_3926.csv')
mz2, intensities2 = read_spectrum('../data/Peptide_spec_demo/EK_Q_07.mzML_scan_6535.csv')

# Define peptide sequence for peptide A

# Create Spectrum objects
spectrum1 = sus.MsmsSpectrum("spectrum1", 430.9007,3, mz1, intensities1)
spectrum2 = sus.MsmsSpectrum("spectrum2", 430.8949,3, mz2, intensities2)
peptide = "SVGDLTSADLKGK"

spectrum1.annotate_proforma(peptide, 1, "Da", ion_types="aby")




fig, ax = plt.subplots(figsize=(10,6))
sup.mirror(spectrum1, spectrum2, ax=ax)
# Customize the plot
ax.legend(['Spectrum Identified', 'Spectrum Unidentified'])
ax.set_title('Mirror Plot of the Identified and Unidentified Spectra')
ax.set_xlabel('m/z')
ax.set_ylabel('Intensity')
plt.tight_layout()
plt.show()