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
            if not toi_name in TOI_LIST: #this line just checks to see if 
                continue
            else:
                if file_counter(directory, '*ccf*'):
                    make_toi_dir(toi_name)
                    data = create_wobble_data(toi_name)
                    plot_pipeline_rvs(data, toi_name)
                    results = wobble_analysis(data, toi_name)
                    plot_star_rvs(results, toi_name)

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

def calibration_file_check(directory):
    wave_file_string = '*_wave_A*'
    ccf_file_string = '*_ccf_??_A*'
    e2ds_file_string = '*_e2ds_*'
    wave_file_count = file_counter(directory, wave_file_string)
    ccf_file_count = file_counter(directory, ccf_file_string)
    e2ds_file_count = file_counter(directory, e2ds_file_string)
    return wave_file_count == ccf_file_count == e2ds_file_count



"""
Periodogram function lines2221_01
period = 1/freq
plt.plot(period, power)
plt.xlim([0,10])
plt.axvline(x=1.2089819, c='r')
plt.axvline(x=3.6480957, c='g')
plt.axvline(x=6.2014698, c='b')
plt.title('Periodogram for TOI 4017 all observations')
plt.xlabel('Period')
plt.ylabel('Periodogram Power')
"""
    
if __name__ == '__main__':
    main()