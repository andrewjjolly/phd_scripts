import eleanor

eleanor.Update(sector=42)

star = eleanor.multi_sectors(tic=301289516, sectors='all')

print('Found TIC {0} (Gaia {1}), with TESS magnitude {2}, RA {3}, and Dec {4}'
     .format(star.tic, star.gaia, star.tess_mag, star.coords[0], star.coords[1]))

#data = eleanor.TargetData(star, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True, regressors='corner')

#plt.figure(figsize=(15,5))

data = []
#plot_fmt = ['k.', 'r.','k.', 'r.']

for s in star:
    datum = eleanor.TargetData(s, height=15, width=15, bkg_size=31, do_psf=False, do_pca=False, aperture_mode='small')
    data.append(datum)

time = []
flux = []
flux_err = []
background = []
for sector, datum in enumerate(data):
    q = datum.quality == 0
    time.append(datum.time[q])
    flux.append(datum.corr_flux[q]/np.median(datum.corr_flux[q]))
    flux_err.append(datum.flux_err[q])
    background.append(datum.flux_bkg[q])


print("Data downloaded correctly")

    ###create lightcurve object

plt.figure()
lc = lk.LightCurve(time = time, flux = flux)
lc.plot()
plt.savefig('lightcurve')