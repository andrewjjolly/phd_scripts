import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import exoplanet as xo
from pathlib import Path
import pandas as pd
import etta
import pymc3 as pm
import pymc3_ext as pmx
import aesara_theano_fallback.tensor as tt
import arviz as az
import dataframe_image as dfi

#Wobble results for TOI 4517_01

RESULTS_DIR = '/srv/scratch/astro/z5345592/results/tess_toi/4517_01'
RESULTS_PATH = Path(RESULTS_DIR) / '4517_01_wobble_rvs.txt'
RESULTS_TABLE = pd.read_csv(RESULTS_PATH, sep=' ', header=3)

TOI_INFO_TABLE = etta.download_toi('4517.01')

#getting the useful information out of the tables.

dates = np.array(RESULTS_TABLE['dates'])
rvs_norm = np.array(RESULTS_TABLE['RV'] - np.median(RESULTS_TABLE['RV']))
rvs_err = np.array(RESULTS_TABLE['RV_err'])

x = dates
x_ref = 0.5 * (x.min() + x.max()) # this is needed for some pymc3 stuff later
t = np.linspace(x.min() - 5, x.max() + 5, 50000)

toi_period = TOI_INFO_TABLE['Period (days)']
toi_period_err = TOI_INFO_TABLE['Period (days) err']
toi_epoch = TOI_INFO_TABLE['Epoch (BJD)']
toi_epoch_err = TOI_INFO_TABLE['Epoch (BJD) err']

Ks = xo.estimate_semi_amplitude(toi_period, dates, rvs_norm, yerr=rvs_err)

t0s = toi_epoch
t0_errs = toi_epoch_err
periods = toi_period
period_errs = toi_period_err
num_planets = 1
yerr = rvs_err
x = dates
y = rvs_norm

dates_2 = [2458045, 2458051]

with pm.Model() as model:

    # Gaussian priors
    t0 = pm.Normal("t0", mu=np.array(t0s), sd=np.array(3*t0_errs), shape=num_planets)
    # t0 = pm.Uniform("t0", lower = float(toi_epoch - toi_epoch_err), upper = float(toi_epoch + toi_epoch_err))
    logP = pm.Normal("logP", mu=np.log(periods), sd=np.array(period_errs) / np.array(periods), shape=num_planets, testval=np.log(periods))
    P = pm.Deterministic("P", tt.exp(logP))

    # Wide log-normal prior for semi-amplitude
    logK = pm.Normal("logK", mu=np.log(Ks), sd=2.0, shape=num_planets, testval=np.log(Ks))

    # Eccentricity & argument of periasteron
    ecs = pmx.UnitDisk("ecs", shape=(2, num_planets), testval=0.01 * np.ones((2, num_planets)))
    ecc = pm.Deterministic("ecc", tt.sum(ecs ** 2, axis=0))
    omega = pm.Deterministic("omega", tt.arctan2(ecs[1], ecs[0]))
    xo.eccentricity.vaneylen19("ecc_prior", multi=True, shape=num_planets, fixed=True, observed=ecc)

    # Jitter & a quadratic RV trend
    logs = pm.Normal("logs", mu=np.log(np.median(yerr)), sd=5.0)
    trend = pm.Normal("trend", mu=0, sd=10.0 ** -np.arange(2)[:-1], shape=2)

    # Then we define the orbit
    orbit = xo.orbits.KeplerianOrbit(period=P, t0=t0, ecc=ecc, omega=omega)

    # And a function for computing the full RV model
    def get_rv_model(t, name=""):

        # First the RVs induced by the planets
        vrad = orbit.get_radial_velocity(t, K=tt.exp(logK))
        pm.Deterministic("vrad" + name, vrad)

        # Define the background model
        A = np.vander(t - x_ref, 2)
        bkg = pm.Deterministic("bkg" + name, tt.dot(A, trend))

        # Sum over planets and add the background to get the full model
        return pm.Deterministic("rv_model" + name, vrad + bkg)

    # Define the RVs at the observed times
    rv_model = get_rv_model(x)

    # Also define the model on a fine grid as computed above (for plotting)
    rv_model_pred = get_rv_model(t, name="_pred")

    # Finally add in the observation model. This next line adds a new contribution
    # to the log probability of the PyMC3 model
    err = tt.sqrt(yerr ** 2 + tt.exp(2 * logs))
    pm.Normal("obs", mu=rv_model, sd=err, observed=y)


plt.figure()
plt.errorbar(x, y, yerr=yerr, fmt=".k")

with model:
    plt.plot(t, pmx.eval_in_model(model.vrad_pred), "--k", alpha=0.5)
    plt.plot(t, pmx.eval_in_model(model.bkg_pred), ":k", alpha=0.5)
    plt.plot(t, pmx.eval_in_model(model.rv_model_pred), label="model")

plt.legend(fontsize=10)
plt.xlim(dates_2)
plt.xlabel("Date [MJD]")
plt.ylabel("Radial Velocity [m/s]")
plt.title("Initial model")
plt.savefig('initial_model.png')

with model:
    map_soln = pmx.optimize(start=model.test_point, vars=[trend])
    map_soln = pmx.optimize(start=map_soln, vars=[t0, trend, logK, logP, logs])
    map_soln = pmx.optimize(start=map_soln, vars=[ecs])
    map_soln = pmx.optimize(start=map_soln)

plt.figure()
plt.errorbar(x, y, yerr=yerr, fmt=".k")
plt.plot(t, map_soln["vrad_pred"], "--k", alpha=0.5)
plt.plot(t, map_soln["bkg_pred"], ":k", alpha=0.5)
plt.plot(t, map_soln["rv_model_pred"], label="model")

plt.legend(fontsize=10)
plt.xlim(dates_2)
plt.xlabel("Date [MJD]")
plt.ylabel("Radial Velocity [m/s]")
plt.title("MAP model")
plt.savefig('MAP_model.png')

np.random.seed(42)
with model:
    trace = pmx.sample(
        tune=1000,
        draws=1000,
        cores=4,
        chains=4,
        target_accept=0.9,
        return_inferencedata=True,)

sample_summary = az.summary(trace, var_names=["trend", "logs", "omega", "ecc", "t0", "logK", "P"])
sample_summary.to_csv('sample_summary.csv')