import matplotlib.pyplot as plt
from pyteomics import mass
import pandas as pd

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
masses_B, intensities_B = read_spectrum('../data/Peptide_spec_demo/EK_Q_07_2.mzML_scan_4915.csv')
masses_A, intensities_A = read_spectrum('../data/Peptide_spec_demo/EK_Q_11.mzML_scan_5622.csv')

# Define peptide sequence for peptide A
peptide_A = 'VLALDVYNDKIK'
charge_state = 1

# Calculate b+ and y+ ions for peptide A
b_ions_A = []
y_ions_A = []
a_ions_A = []

for i in range(1, len(peptide_A)):
    a_ions_A.append(mass.calculate_mass(sequence=peptide_A[:i], ion_type='a', charge=charge_state))
    b_ions_A.append(mass.calculate_mass(sequence=peptide_A[:i], ion_type='b', charge=charge_state))
    y_ions_A.append(mass.calculate_mass(sequence=peptide_A[i:], ion_type='y', charge=charge_state))

# Plotting

# Plotting
plt.figure(figsize=(12, 8))

# Plot peptide A ions
# Plot peptide A ions in grey without markers
plt.stem(masses_A, intensities_A, linefmt='grey', markerfmt=' ', basefmt='grey', label='Spectrum A')

# Plot peptide B ions in grey without markers (mirror by negating the intensities)
plt.stem(masses_B, -intensities_B, linefmt='grey', markerfmt=' ', basefmt='grey', label='Spectrum B')

# Annotate the b+ and y+ ions for peptide A
for i, mass in enumerate(b_ions_A):
    closest_mass = find_closest_mass(masses_A, mass, tolerance=0.01)
    if closest_mass is not None:
        intensity_index = list(masses_A).index(closest_mass)
        intensity = intensities_A[intensity_index]
        plt.stem([closest_mass], [intensity], linefmt='green', markerfmt='go', basefmt=' ', label=f'b{i+1}+' if i == 0 else "")
        plt.annotate(f'b{i+1}+', (closest_mass, intensity), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8, color='green')

for i, mass in enumerate(y_ions_A):
    closest_mass = find_closest_mass(masses_A, mass, tolerance=0.01)
    if closest_mass is not None:
        intensity_index = list(masses_A).index(closest_mass)
        intensity = intensities_A[intensity_index]
        plt.stem([closest_mass], [intensity], linefmt='red', markerfmt='ro', basefmt=' ', label=f'y{i+1}+' if i == 0 else "")
        plt.annotate(f'y{i+1}+', (closest_mass, intensity), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8, color='red')

for i, mass in enumerate(a_ions_A):
    closest_mass = find_closest_mass(masses_A, mass, tolerance=0.01)
    if closest_mass is not None:
        intensity_index = list(masses_A).index(closest_mass)
        intensity = intensities_A[intensity_index]
        plt.stem([closest_mass], [intensity], linefmt='orange', markerfmt='yo', basefmt=' ', label=f'a{i+1}+' if i == 0 else "")
        plt.annotate(f'a{i+1}+', (closest_mass, intensity), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8, color='orange')

plt.axhline(0, color='black', linewidth=0.5)

plt.axhline(0, color='black', linewidth=0.5)
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title('Peptide Spectrum Mirror Plot')
plt.legend()
plt.grid(True)
plt.show()
