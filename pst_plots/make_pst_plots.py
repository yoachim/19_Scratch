import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
import lsst.sims.maf.plots as plots
import healpy as hp
import subprocess


nside = 128
# Connect to an opsim database, bleeding-edge sims here:
# https://lsst-web.ncsa.illinois.edu/sim-data/sims_featureScheduler_runs/
opsdb = db.OpsimDatabase('baseline_1exp_pairsmix_10yrs.db')
outDir = 'maf_out'
resultsDb = db.ResultsDb(outDir=outDir)
plotFuncs = [plots.TwoDMap()]


filters = ['u', 'g', 'r', 'i', 'z', 'y']
day_max = np.round(365.25*2)
bins = np.arange(day_max)

nval = 3

for filtername in filters:
    metric = metrics.AccumulateCountMetric(bins=bins)
    slicer = slicers.HealpixSlicer(nside=nside)
    plotDict = {'xlabel': 'Night (days)', 'cbarTitle': 'N obs', 'colorMax': 75}
    # only use i-band
    sql = 'filter = "%s" and night < %i' % (filtername, day_max)
    bundle = metricBundles.MetricBundle(
        metric, slicer, sql, plotDict=plotDict, plotFuncs=plotFuncs)
    group = metricBundles.MetricBundleGroup(
        {0: bundle}, opsdb, outDir=outDir, resultsDb=resultsDb)
    group.runAll()

    n_val = 3

    # Ugh, can't figure out how to do this properly with slicing
    time_to_val = np.zeros(bundle.metricValues[:, 0].size)
    first_obs = np.zeros(time_to_val.size)
    for i in np.arange(time_to_val.size):
        good = np.where(bundle.metricValues[i, :] > 0)[0]
        if np.size(good) > 0:
            first_obs[i] = bins[np.min(good)]
        else:
            first_obs[i] = hp.UNSEEN
        good = np.where(bundle.metricValues[i, :] > n_val)[0]
        if good.size > 0:

            time_to_val[i] = bins[good.min()]
        else:
            time_to_val[i] = hp.UNSEEN

    hp.mollview(time_to_val, title='Time to %i obs' % n_val, unit='Days', max=365)
    plt.savefig('ttv_%s.pdf' % filtername)
    diff = time_to_val-first_obs

    mask = np.where((time_to_val == hp.UNSEEN) | (first_obs == hp.UNSEEN))
    diff[mask] = hp.UNSEEN
    hp.mollview(diff, title='Time to %i obs after first obs' % n_val, max=200)
    plt.savefig('ttn_afterfirst_%s.pdf' % filtername)

    outDir = 'movie_plots_%s' % filtername
    for i in np.arange(bundle.metricValues[0, :].size):
        hp.mollview(bundle.metricValues[:, i], max=80, unit='N obs', title='%s, Night %i' % (filtername, bins[i]))
        plt.savefig(outDir+'/%i.png' % i)
        plt.close()

    subprocess.call("ffmpeg -framerate 4 -pix_fmt yuv420p  -i movie_plots_%s/%%d.png %s_band.mp4" % (filtername, filtername), shell=True)
