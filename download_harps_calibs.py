from pathlib import Path
from astropy.io import fits
import eso_programmatic as eso
import pandas as pd
from tqdm import tqdm

DATA_DIR = Path('/home/z5345592/data/tess_toi')

obs_no_wave = []

for dir in DATA_DIR.iterdir():
    
    toi_path = DATA_DIR / dir.name
    print(f'Getting calibration file details for TOI {dir.name}')

    for file in toi_path.rglob('*e2ds_A.fits'):

        try:
            header = fits.open(file)[0].header
            wave_file = header['HIERARCH ESO DRS CAL TH FILE']
            wave_file = wave_file.split('_')[0]
            archive_id = wave_file = wave_file + '.tar'

            toi_path = '/home/z5345592/katana_copy/scripts/test_folder'
            data_url = "https://dataportal.eso.org/dataportal_new/file/%s" % archive_id

            status_code, filepath = eso.downloadURL(data_url, dirname = str(dir))

            if status_code == 200:
                print("File {0} downloaded as {1}".format(archive_id, filepath))
            else:
                print("ERROR: could not download %s" % (archive_id))

        except KeyError:
            obs_no_wave.append(file)
            print(f'No valid wavelength calibration file for {file.name}')

def download_calibration_files(dir):

    for file in dir.rglob('e2ds_A.fits'):

        try:
            header = fits.open(file)[0].header
            wave_file = header['HIERARCH ESO DRS CAL TH FILE']
            wave_file = wave_file.split('_')[0]
            archive_id = wave_file + '.tar'

            data_url = "https://dataportal.eso.org/dataportal_new/file/%s" % archive_id

            status_code, filepath = eso.downloadURL(data_url, dirname = str(dir))

            if status_code == 200:
                print("File {0} downloaded as {1}".format(archive_id, filepath))
            else:
                print("ERROR: could not download %s" % (archive_id))

        except KeyError:
            obs_no_wave.append(file)
            # print(f'No valid wavelength calibration file for {file.name}')