"""
1. in a folder, identify where the files with the pipeline RVs are, generate and plot them, save the plots.
2. wobble analysis
3. periodogram
"""

import wobble
from pathlib import Path
import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import fits

DATA_DIR = '/srv/scratch/astro/z5345592/data/tess_toi/'
RESULTS_DIR = '/srv/scratch/astro/z5345592/results/tess_toi'

TOI_LIST = ['5094_01'] #added this list so I could run analysis on just a subsection.

def main():

    for directory in Path(DATA_DIR).iterdir():
        if directory.is_dir:
            toi_name = directory.name
            if not toi_name in TOI_LIST: #this line is just for when I am working with a subset of TOIs
                continue
            else:
                if file_counter(directory, '*ccf*'):
                    make_toi_dir(toi_name)
                    data = create_wobble_data(toi_name)
                    plot_pipeline_rvs(data, toi_name)
                    results = wobble_analysis(data, toi_name)
                    plot_star_rvs(results, toi_name)

def file_check(directory): # a function that takes in a directory and returns the ccfs that have matching spectrum files that also have the correct wave calibs
    
    ccf_file_list = []
    wave_file_list = []
    valid_ccf_files = []
    obs_id_list = []

    for file in Path(directory).rglob('*_ccf_??_A.fits'):
        stellar_type = (file.name)[-9:-7]
        if not any(stellar_type == str for str in ['M2', 'K5', 'G2']): #checks to pick only the ccfs that relate to a stellar observation
            continue
        obs_id = file.name[:29] #slicing the HARPS observation ID
        ccf_file_list.append(str(file)) #making the list of ccfs
        obs_id_list.append(obs_id) #assigning the observation ID to the list

    for file in Path(directory).rglob('*_wave_A.fits'): #loop through the wave files in the directory
        wave_file_list.append(str(file.name)) #assign them to a list

    for file in Path(directory).rglob('*_e2ds_A.fits'):
        phdu = fits.open(file) #primary header data unit of the fits file
        header = phdu[0].header #the header
        wave_file = header['HIERARCH ESO DRS CAL TH FILE'] #wave file filename as found in the spectrum file header
        if wave_file in wave_file_list: #checks to see whether the wave file is in the files in the directory.
            if (file.name[:29]) in obs_id_list: #checks to make sure that the observation ID for the spectrum matches the observation IDs taken from the ccf.
                valid_ccf = [str for str in valid_ccf_files if file.name[:29] in str]
                valid_ccf = valid_ccf[0] #because valid ccf on the line above is a list, but I only need the one.
                valid_ccf_files.append(valid_ccf)

    return valid_ccf_files

def check_ccf_science(ccf_file):
    phdu = fits.open(ccf_file)
    cat = phdu[0].header['HIERARCH ESO DPR CATG']
    return cat == 'SCIENCE'

def make_toi_dir(toi_name):
    path = Path(RESULTS_DIR) / toi_name
    if not path.exists():
        path.mkdir()
        print('\n')
        print(f'Created results directory for TESS TOI {toi_name}.')
    else:
        print('\n')
        print(f'Results directory for {toi_name} already exists.')

    return

def create_wobble_data(toi_name):
    print('\n')
    print(f'Creating wobble data object for TOI {toi_name}.')
    data = wobble.Data()
    sp = wobble.Spectrum()
    path = Path(DATA_DIR) / toi_name
    for ccf_file in path.rglob('*ccf_M?_A.fits'):
        if check_ccf_science(ccf_file):
            sp.from_HARPS(str(ccf_file), process=False)
            data.append(sp)
        # data.drop_bad_epochs()
        # data.drop_bad_orders()

    return data

def plot_pipeline_rvs(data, toi_name):
    dates = data.dates
    rvs = data.pipeline_rvs
    rvs_err = data.pipeline_sigmas
    plt.figure()
    plt.errorbar(dates, rvs, yerr=rvs_err, fmt='k.')
    plt.title(f'HARPS Pipeline Radial Velocity for TOI {toi_name}')
    plt.xlabel('Date, [MJD]')
    plt.ylabel('Radial Velocity, [m/s]')
    plot_filename = f'{toi_name}_plrvs_harps.png'
    plot_path = str(Path(RESULTS_DIR) / toi_name / plot_filename)
    plt.savefig(f'{plot_path}', dpi=300)

def wobble_analysis(data, toi_name):
    results = wobble.Results(data = data)

    for r in range(len(data.orders)):
            
        print(f'starting order {r+1} of {len(data.orders)}')
        model = wobble.Model(data, results, r)
        model.add_star('star', learning_rate_rvs = 20)
        model.add_telluric('tellurics', variable_bases = 2)
        wobble.optimize_order(model, verbose=False)

    if len(data.orders) > 1:
        results.combine_orders('star')
    
    results.apply_bervs('star')
    results.apply_drifts('star')

    results.write_rvs('star', str(Path(RESULTS_DIR) / toi_name / f'{toi_name}_wobble_rvs.txt'), all_orders=True)
    results.write(str(Path(RESULTS_DIR) / toi_name / f'{toi_name}_wobble_results.hdf5'))

    for r in range(len(data.orders)):
        n = int(np.median(data.epochs)) #median epoch to plot
        results_path = Path(RESULTS_DIR) / toi_name
        filename = f'{toi_name}_spec_order_{r}.png'
        filepath = Path(results_path) / filename
        results.plot_spectrum(r, n, data, filepath)

    return results

def plot_star_rvs(results, toi_name):
    wobble_rvs = results.star_time_rvs
    wobble_rvs_err = results.star_time_sigmas
    pipeline_rvs = results.pipeline_rvs
    pipeline_rvs_err = results.pipeline_sigmas
    plt.figure()
    plt.subplot(2,1,1)
    plt.errorbar(results.dates, pipeline_rvs, yerr = pipeline_rvs_err, fmt='r,', label='HARPS Pipeline')
    plt.subplot(2,1,2)
    plt.errorbar(results.dates, wobble_rvs, yerr = wobble_rvs_err, fmt='k,', label='wobble')
    filename = f'{toi_name}_rvs.png'
    results_path = Path(RESULTS_DIR) / toi_name
    filepath = Path(results_path) / filename
    plt.savefig(filepath, dpi=300)    

def periodogram():
    return

def light_curve():
    return

def file_counter(folder_name, file_string):
    count = 0
    for file in Path(folder_name).rglob(file_string):
        if file.is_file:
            count += 1
    return count
    
if __name__ == '__main__':
    main()