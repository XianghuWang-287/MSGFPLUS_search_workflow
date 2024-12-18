import os
from pyopenms import *

from urllib.request import urlretrieve

def analyze_features_and_ms2(mzml_file, output_featurexml):
    """
    Perform untargeted metabolomics feature detection using FeatureFindingMetabo
    and analyze how many MS2 acquisitions are associated with each feature.

    Args:
        mzml_file (str): Path to the mzML file.
        output_featurexml (str): Path to save the featureXML output.

    Returns:
        None
    """
    # Load the mzML file
    print("Loading mzML file...")
    exp = MSExperiment()
    MzMLFile().load(mzml_file, exp)
    exp.sortSpectra(True)  # Ensure spectra are sorted by RT
    print("mzML file loaded and sorted successfully.")

    # Mass trace detection
    print("Performing mass trace detection...")
    mass_traces = []
    mtd = MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue("mass_error_ppm", 5.0)  # Set according to your instrument mass error
    mtd_params.setValue("noise_threshold_int", 3000.0)  # Adjust based on noise level in your data
    mtd.setParameters(mtd_params)
    mtd.run(exp, mass_traces, 0)
    print(f"Detected {len(mass_traces)} mass traces.")

    # Elution peak detection
    print("Detecting elution peaks...")
    mass_traces_split = []
    mass_traces_final = []
    epd = ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd_params.setValue("width_filtering", "fixed")
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)

    if epd.getParameters().getValue("width_filtering") == "auto":
        epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
    else:
        mass_traces_final = mass_traces_split
    print(f"Detected {len(mass_traces_final)} refined elution peaks.")

    # Feature finding
    print("Running FeatureFindingMetabo...")
    fm = FeatureMap()
    feat_chrom = []
    ffm = FeatureFindingMetabo()
    ffm_params = ffm.getDefaults()
    ffm_params.setValue("isotope_filtering_model", "none")
    ffm_params.setValue("remove_single_traces", "true")  # Keep or remove single mass traces
    ffm_params.setValue("mz_scoring_by_elements", "false")
    ffm_params.setValue("report_convex_hulls", "true")
    ffm.setParameters(ffm_params)
    ffm.run(mass_traces_final, fm, feat_chrom)
    print(f"Detected {fm.size()} features.")

    # Save features to file
    fm.setUniqueIds()
    fm.setPrimaryMSRunPath([mzml_file.encode()])  # Required for linking data
    FeatureXMLFile().store(output_featurexml, fm)
    print(f"Features saved to {output_featurexml}.")

    # Analyze MS2 mapping to features
    analyze_ms2_mapping(exp, fm)


# def analyze_ms2_mapping(exp, feature_map):
#
#     print("Mapping MS2 spectra to detected features...")
#
#     # Prepare a list of features with their m/z and RT range
#     features = []
#     for feature in feature_map:
#         mz = feature.getMZ()
#         rt_start = feature.getRT() - feature.getWidth() / 2.0
#         rt_end = feature.getRT() + feature.getWidth() / 2.0
#         features.append((mz, rt_start, rt_end))
#
#     # Count how many MS2 spectra map to each feature
#     ms2_feature_count = [0] * len(features)
#     for spectrum in exp:
#         if spectrum.getMSLevel() == 2:
#             precursor_mz = spectrum.getPrecursors()[0].getMZ()
#             precursor_rt = spectrum.getRT()
#             for i, (mz, rt_start, rt_end) in enumerate(features):
#                 if abs(precursor_mz - mz) <= 0.01 and rt_start <= precursor_rt <= rt_end:
#                     ms2_feature_count[i] += 1
#                     break
#
#     # Report statistics
#     features_with_multiple_ms2 = sum(1 for count in ms2_feature_count if count > 1)
#     total_features = len(features)
#     print(f"Total detected features: {total_features}")
#     print(f"Features with multiple MS2 acquisitions: {features_with_multiple_ms2}")
#     print(f"Percentage with multiple MS2: {features_with_multiple_ms2 / total_features * 100:.2f}%")

def analyze_ms2_mapping(exp, feature_map):
    """
    Map MS2 spectra to detected features and report statistics.

    Args:
        exp (MSExperiment): The MSExperiment object containing the spectra.
        feature_map (FeatureMap): The detected features.

    Returns:
        None
    """
    print("Mapping MS2 spectra to detected features...")

    # Prepare a list of features with their m/z and RT range
    features = []
    for feature in feature_map:
        mz = feature.getMZ()
        rt_start = feature.getRT() - feature.getWidth() / 2.0
        rt_end = feature.getRT() + feature.getWidth() / 2.0
        features.append((mz, rt_start, rt_end))

    # Create a mapping of MS2 spectra to features
    ms2_to_feature = []
    for spectrum in exp:
        if spectrum.getMSLevel() == 2:
            precursor_mz = spectrum.getPrecursors()[0].getMZ()
            precursor_rt = spectrum.getRT()
            for i, (mz, rt_start, rt_end) in enumerate(features):
                if abs(precursor_mz - mz) <= 0.01 and rt_start <= precursor_rt <= rt_end:
                    ms2_to_feature.append(i)
                    break

    # Count how many MS2 spectra are associated with each feature
    from collections import Counter
    feature_ms2_counts = Counter(ms2_to_feature)

    # Calculate the total MS2 spectra associated with features having multiple MS2
    total_ms2 = len(ms2_to_feature)
    ms2_from_multiple = sum(count for count in feature_ms2_counts.values() if count > 1)

    # Report statistics
    print(f"Total MS2 spectra: {total_ms2}")
    print(f"MS2 spectra from multiple acquisitions: {ms2_from_multiple}")
    print(f"Percentage of MS2 spectra from multiple acquisitions: {ms2_from_multiple / total_ms2 * 100:.2f}%")


if __name__ == "__main__":
    mzml_file = "/home/user/research/MSGFPLUS_search_workflow/data/Feature_find/isa_9.mzML"  # Replace with your mzML file path
    output_featurexml = "output.featureXML"  # Replace with desired output path
    analyze_features_and_ms2(mzml_file, output_featurexml)



