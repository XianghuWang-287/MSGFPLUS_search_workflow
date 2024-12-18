import matplotlib.pyplot as plt
from pyteomics import mass
import pandas as pd
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import math
import matplotlib.ticker as mticker
from typing import Dict, Optional
# Function to read the spectrum data from a CSV file
def read_spectrum(csv_file):
    data = pd.read_csv(csv_file)
    return data['mz'].values, data['intensity'].values

def find_closest_mass(mass_list, target_mass, tolerance):
    for mass in mass_list:
        if abs(mass - target_mass) <= tolerance:
            return mass
    return None

def mirror(
    spec_top: sus.MsmsSpectrum,
    spec_bottom: sus.MsmsSpectrum,
    spectrum_kws: Optional[Dict] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Mirror plot two MS/MS spectra.

    Parameters
    ----------
    spec_top : MsmsSpectrum
        The spectrum to be plotted on the top.
    spec_bottom : MsmsSpectrum
        The spectrum to be plotted on the bottom.
    spectrum_kws : Optional[Dict], optional
        Keyword arguments for `plot.spectrum`.
    ax : Optional[plt.Axes], optional
        Axes instance on which to plot the spectrum. If None the current Axes
        instance is used.

    Returns
    -------
    plt.Axes
        The matplotlib Axes instance on which the spectra are plotted.
    """
    if ax is None:
        ax = plt.gca()

    if spectrum_kws is None:
        spectrum_kws = {}
    # Top spectrum.
    sup.spectrum(spec_top, mirror_intensity=False,grid = False, ax=ax, **spectrum_kws)
    y_max = ax.get_ylim()[1]
    # Mirrored bottom spectrum.
    sup.spectrum(spec_bottom, mirror_intensity=True,grid = False, ax=ax, **spectrum_kws)
    y_min = ax.get_ylim()[0]
    ax.set_ylim(y_min, y_max)

    ax.axhline(0, color="#9E9E9E", zorder=10)

    max_mz_top = spec_top.mz[-1] if len(spec_top.mz) > 0 else 1
    max_mz_bottom = spec_bottom.mz[-1] if len(spec_bottom.mz) > 0 else 1
    # Update axes so that both spectra fit.
    round_mz = 50
    max_mz = max(
        [
            math.ceil(max_mz_top / round_mz + 1) * round_mz,
            math.ceil(max_mz_bottom / round_mz + 1) * round_mz,
        ]
    )
    ax.set_xlim(0, max_mz)
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, pos: f"{abs(x):.0%}")
    )

    # Disable the grid for the main axis
    ax.grid(False)

    # Disable the grid for all secondary axes if any
    for secondary_ax in ax.figure.get_axes():
        secondary_ax.grid(False)

    return ax

# Read the spectra for peptides A and B
mz1, intensities1 = read_spectrum('../data/Peptide_spec_demo/EK_Q_07.mzML_scan_3926.csv')
mz2, intensities2 = read_spectrum('../data/Peptide_spec_demo/EK_Q_07.mzML_scan_6535.csv')

# Define peptide sequence for peptide A

# Create Spectrum objects
spectrum1 = sus.MsmsSpectrum("spectrum1", 430.9007,3, mz1, intensities1)
spectrum2 = sus.MsmsSpectrum("spectrum2", 430.8949,3, mz2, intensities2)
peptide = "SVGDLTSADLKGK"

spectrum1.annotate_proforma(peptide, 1, "Da", ion_types="aby",)




fig, ax = plt.subplots(figsize=(8,6))
# Disable the grid for all secondary axes

ax = mirror(spectrum1, spectrum2, ax=ax)

ax.grid(False)
# Explicitly disable the grid for all axes in the figure
for axis in fig.get_axes():
    axis.grid(False)
# Customize the plot
ax.legend(['Spectrum Identified', 'Spectrum Unidentified'],fontsize=14)
ax.set_title('Mirror Plot of the Identified and Unidentified Spectra',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel('m/z',fontsize=14)
ax.set_ylabel('Intensity',fontsize=14)
current_ylim = plt.ylim()
plt.ylim(-1.2, 1.2)
plt.xlim(right = 1000)
plt.grid(False)

# plt.tight_layout()
plt.savefig("./pep_mirror.png", bbox_inches="tight", dpi=300, transparent=True)
plt.show()
