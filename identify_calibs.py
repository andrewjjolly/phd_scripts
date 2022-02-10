from pathlib import Path
from astropy.io import fits
import shutil

TEST_TOI = '4304_01'
DATA_DIR = Path('/srv/scratch/astro/z5345592/data/data/tess_toi')
TOI_DIR = DATA_DIR / TEST_TOI
CALIB_DIR = 'srv/scratch/astro/z5345592/data/data/harps_calibration_files'

for file in TOI_DIR.rglob('*e2ds_A*.fits'):
    sp = fits.open(file)
    wave_file = sp[0].header['HIERARCH ESO DRS CAL TH FILE']
    try:
        origin = Path(CALIB_DIR) / wave_file
        destination = Path(TOI_DIR) / wave_file
        shutil.copy(origin, destination)
    except FileNotFoundError:
        print(f"{wave_file} not in the calibration files directory")
    except KeyError:
        print('Wave file not found in the header')