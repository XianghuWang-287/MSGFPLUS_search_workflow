
from pyteomics import mzml
import matplotlib.pyplot as plt


mzml_file = '/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_convert/mzML/EK_Q_11.mzML'


scan_number = 3695



def plot_spectrum(mzml_file, scan_number):
    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            scan_id = int(spectrum['id'].split('=')[-1])

            if scan_id == scan_number:
                mzs = spectrum['m/z array']
                intensities = spectrum['intensity array']

                plt.bar(mzs, intensities, width=1)  # Adjust the width as needed
                plt.xlabel('m/z')
                plt.ylabel('Intensity')
                plt.title(f'Spectrum for Scan {scan_number}')
                plt.show()
                return

    print(f"Scan {scan_number} not found in the mzML file.")

plot_spectrum(mzml_file, scan_number)

