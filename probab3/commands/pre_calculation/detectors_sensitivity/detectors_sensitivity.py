import numpy as np
from matplotlib import pyplot as plt

detectors_sensitivity = {
    'LIGO A+': {'freq_range': [34.375506, 1814.7771],
                'min_strain': 1.5367463e-24,
                'strain_threshold': 4.6102389e-24},
    'Cosmic Explorer': {'freq_range': [10.624226762494413, 1899.1076530953678], 
                        'min_strain': 2.138284935904749e-25,
                        'strain_threshold': 6.414854807714247e-25},
    'Voyager': {'freq_range': [25.511, 1204.8],
                'min_strain': 8.3629e-25,
                'strain_threshold': 2.50887e-24},
    'LIGO (O3)': {'freq_range': [35.25, 1764.25],
                  'min_strain': 3.8613662e-24,
                  'strain_threshold': 1.1584098599999998e-23}
}


def get_detector_sensitivity(detector_name, file_path, delimiter=' ', skiprows=0):

    file_data = np.loadtxt(file_path, delimiter=delimiter, skiprows=skiprows)
    frequency = file_data[:,0]
    strain = file_data[:,-1]
    min_strain = np.min(strain)
    strain_threshold = 3*min_strain
    detectable_freq = frequency[strain <= strain_threshold]
    freq_range = [detectable_freq[0], detectable_freq[-1]]

    plt.scatter(frequency, strain, label=detector_name)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Strain noise [1/sqrt(Hz)]')
    plt.legend()

    return {detector_name: {'freq_range': freq_range, 'min_strain': min_strain, 'strain_threshold': strain_threshold}}