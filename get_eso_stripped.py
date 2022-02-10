"""
A stripped down version of get_eso to see if I can get away with just downloading the
files via 1 sort of loop. Don't use this script for actual work.
"""

import os
from astropy.extern.configobj.validate import _test
from astropy.io import fits
import numpy as np
import subprocess
import re
import time
import glob
from bs4 import BeautifulSoup
from astroquery.eso import Eso as ESO
import tarfile
from tqdm import tqdm
import sys
from astropy.visualization import astropy_mpl_style
from astropy import table
from astropy.coordinates import SkyCoord
from astropy.units import Quantity
from astropy import units as u
from pyvo.dal import tap
import pandas as pd
import requests
import cgi
import etta
from pyvo.utils.http import create_session
import pyvo
import eso_programmatic as esoprog
from pathlib import Path
import shutil

print("\n") #new line because of all the various errors you get from tensorflow importing

ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
ESO_USERNAME = "andrewjolly"
tapobs = tap.TAPService(ESO_TAP_OBS)

DATA_DIR = '/srv/scratch/astro/z5345592/data/tess_tois/'

arr_ind = int(sys.argv[1])
#print(arr_ind)
beg_ind = (arr_ind-1)*10
end_ind = beg_ind+10

# beg_ind = 0
# end_ind = 5

def main():

    tess_toi_df = get_tess_toi_df()

    for toi in tess_toi_df.index[beg_ind:end_ind]:

        toi_name = get_toi_name(toi)

        if not check_toi_dir(toi_name):
        
            toi_coords = get_coords(tess_toi_df, toi)

            if harps_data_available(toi_coords) > 0:
                
                data_dir = make_toi_dir(toi_name)
                download_ancillary_files(toi_coords, data_dir)
                extract_files(data_dir)
                move_files(data_dir)
                download_calibration_files(data_dir)
                extract_files(data_dir)
                move_files(data_dir)
                cleanup(data_dir)
                final_count(data_dir)

        else:
            print(f'Already have TOI {toi_name}.')

def get_tess_toi_df():
    tess_toi_df = etta.download_toi(sort = 'toi')
    print('data fetched!')
    tess_toi_df = tess_toi_df.sort_index(ascending=False) #puts the most recent TOIs first
    tess_toi_df_rows = tess_toi_df.shape[0]
    new_toi = tess_toi_df.first_valid_index()
    new_toi_date = tess_toi_df.iloc[0]['Date Modified']
    new_toi_comments = tess_toi_df.iloc[0]['Comments']
    print("\n")
    print(f'TESS TOI table has {tess_toi_df_rows} TOIs. The most recent is {new_toi} which was updated on {new_toi_date} with the notes: {new_toi_comments}.')

    return(tess_toi_df)

def get_toi_name(toi):
    toi_name = str(toi)
    toi_name = toi_name.replace('.','_')

    return(toi_name)

def check_toi_dir(toi_name):
    path = Path(DATA_DIR) / toi_name
    
    return path.exists()

def make_toi_dir(toi_name):
    path = os.path.join(os.path.join(DATA_DIR, toi_name + '/'))
    if not os.path.exists(path):
        os.mkdir(path)
        print('\n')
        print(f'Created data directory for TESS TOI {toi_name}.')
    else:
        print('\n')
        print(f'Data directory for {toi_name} already exists.')

    return path

def get_coords(tess_toi_df, toi):
    """gets the RA & Dec of the target to be able to do a cone search."""
    toi_ra = tess_toi_df.loc[toi].RA
    toi_dec = tess_toi_df.loc[toi].Dec
    toi_coords = SkyCoord(ra = toi_ra, dec = toi_dec, unit = (u.hour, u.degree), frame='icrs')

    return toi_coords

def harps_data_available(coords):
    sr = 1/60 

    query = """
    SELECT top 1000 access_url 
    FROM ivoa.ObsCore 
    WHERE instrument_name = 'HARPS'
    AND intersects(s_region, circle('', %f, %f, %f))=1
    """ % (coords.ra.degree , coords.dec.degree, sr)

    res = tapobs.search(query=query, maxrec=1000)

    return len(res)

def download_ancillary_files(coords, data_dir):

    sr = 1/60 

    query = """
    SELECT top 1000 access_url 
    FROM ivoa.ObsCore 
    WHERE instrument_name = 'HARPS'
    AND intersects(s_region, circle('', %f, %f, %f))=1
    """ % (coords.ra.degree , coords.dec.degree, sr)
    
    res = tapobs.search(query=query, maxrec=1000)
    if not len(res) == 0:
        print(f'{len(res)} ancillary files identified, commencing download.')
        print('\n')

        session = create_session() #this needed to be added - otherwise not defined in loop below.
        count = 0
        for rec in res:
            datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(rec['access_url'], session=session) #this line needed pyvo not vo
            ancillaries = datalink.bysemantics('#auxiliary')
            
            for anc in ancillaries:
                status_code, filepath = esoprog.downloadURL(anc['access_url'], data_dir, session=session)
                if not status_code == 200:
                    count += 1
        print(f'{count} files were not available.')

def count_files(directory):
    count = 0
    for path in Path(directory).rglob('*'):
        if path.is_file():
            count += 1
    return count

def extract_files(directory):
    """extracts all tar files in a directory into the same directory"""
    initial_count = count_files(directory)

    for file in Path(directory).rglob('*.tar'):
        tar = tarfile.open(file, 'r')
        tar.extractall(path = directory)
    
    final_count = count_files(directory)

    print(f'Extracted {final_count - initial_count} files to {directory}.')

    return

def move_files(directory):
    """moves all files in all subdirectories to be in the main directory, deletes the subdirectory"""
    print('Moving all files out of subdirectories, removing empty subdirectories')
    for path in Path(directory).rglob('*'):
        if path.is_file():
            destination = Path(directory) / path.name
            shutil.move(path, destination)

    #shutil.rmtree(Path(directory) / 'data') do I need to delete the dir? work out how to do this cleanly.

    return
    
def download_calibration_files(dir):

    count = 0
    
    for file in Path(dir).rglob('*e2ds_A.fits'):

        try:
            header = fits.open(file)[0].header
            wave_file = header['HIERARCH ESO DRS CAL TH FILE']
            wave_file = wave_file.split('_')[0]
            archive_id = wave_file + '.tar'

            data_url = "https://dataportal.eso.org/dataportal_new/file/%s" % archive_id

            status_code, filepath = esoprog.downloadURL(data_url, dirname = str(dir))

            if status_code == 200:
                count += 1
            else:
                print("ERROR: could not download %s" % (archive_id))

        except KeyError:
            obs_no_wave.append(file)
            # print(f'No valid wavelength calibration file for {file.name}')

    print(f"Downloaded {count} calibration files to directory")

def cleanup(dir):
    for file in Path(dir).rglob('*.tar'):
        os.remove(file)

def final_count(dir):
    spec_files = 0
    wave_files = 0
    ccf_files = 0
    for file in Path(dir).rglob('*e2ds_A.fits'):
        spec_files += 1
    for file in Path(dir).rglob('*wave_A.fits'):
        wave_files += 1
    for file in Path(dir).rglob('*ccf_??_A.fits'):
        ccf_files += 1

    print("\n")
    print(f'Downloaded {ccf_files} ccf files, {spec_files} spectrum files and {wave_files} wavelength calibration files for this TOI')

if __name__ == '__main__':
    main()
    