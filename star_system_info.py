from pathlib import Path
from posixpath import dirname
from astropy.io import fits
import etta
import pandas as pd
from tqdm import tqdm

DATA_DIR = Path('/srv/scratch/astro/z5345592/data/tess_tois')

tess_toi_df = etta.download_toi(sort = 'toi')
tess_toi_df = tess_toi_df.sort_index(ascending=False)

def get_num_obs(toi_path):
    num_observations = 0
    for path in toi_path.rglob('*e2ds_A.fits'):
        if path.is_file:
                num_observations += 1
    return num_observations

def get_baseline(toi_path):
    all_obs_dates = []
    for path in toi_path.rglob('*e2ds_A.fits'):
            sp = fits.open(path)
            obs_date = sp[0].header['MJD-OBS']
            all_obs_dates.append(obs_date)
    baseline = int(max(all_obs_dates) - min(all_obs_dates))
    return baseline

def sub_directory_counter(directory):
    count = 0
    for sub_dir in Path(directory).iterdir():
        if sub_dir.is_dir:
            count += 1
    return count

all_toi_names = []
all_num_observations = []
all_baselines = []

for dir in DATA_DIR.iterdir():
    if dir.is_dir:
        toi_name = dir.name
        toi_name = str(dir.name)
        all_toi_names.append(toi_name)
        toi_path = Path(DATA_DIR) / dir
        num_observations = get_num_obs(toi_path)
        all_num_observations.append(num_observations)
        if num_observations < 1:
            baseline = 0
        else:
            baseline = get_baseline(toi_path)
        all_baselines.append(baseline)
        if num_observations > 50:
            print(f'TOI {toi_name} has {num_observations} observations across {baseline} days.')

data_dictionary = {'Names': all_toi_names, 'No. Observations': all_num_observations, 'Baseline (days)': all_baselines}
toi_info_df = pd.DataFrame(data_dictionary)
toi_info_df.to_csv('toi_info.csv', sep=',', index=False)