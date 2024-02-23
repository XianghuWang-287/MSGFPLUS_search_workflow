from collections import defaultdict

def merge_scans(input_file):
    merged_scans = defaultdict(lambda: {"PEPMASS":float,"peaks": [], "intensities": []})

    current_series = None
    current_scan_num = None

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("PEPMASS="):
                PEPMASS = float(line.split('=')[1])
                merged_scans[(current_series, current_scan_num)]["PEPMASS"]=PEPMASS
            elif line.startswith("SCAN="):
                scan_info = line.split('=')[1].split('-')
                if len(scan_info) == 2:
                    current_series, current_scan_num = map(int, scan_info)
            elif line != 'BEGIN IONS' or line != 'END IONS':
                elements = line.split()
                if len(elements) >= 2:
                    peak = float(elements[0])
                    merged_scans[(current_series, current_scan_num)]["peaks"].append(peak)
                    intensity = float(elements[1])
                    merged_scans[(current_series, current_scan_num)]["intensities"].append(intensity)


    for (series, scan_num), data in merged_scans.items():
        PEPMASS = data["PEPMASS"]
        avg_intensity = sum(data["intensities"]) / len(data["intensities"]) if data["intensities"] else 0
        print("BEGIN IONS")
        print(f"PEPMASS={PEPMASS}")
        print("CHARGE=1+")
        print(f"SCAN={series}")
        print(f"TITLE=QCXMS_SCAN={series}")
        print("MSLEVEL=2")
        for peak, intensity in zip(data['peaks'], data["intensities"]):
            print(f"{peak} {avg_intensity if intensity == 0 else intensity}")
        print("END IONS")
        print()
if __name__ == "__main__":
    input_file = "../data/merged_predicted_spectra.mgf"  # Replace with your input file path
    merge_scans(input_file)

