#!/usr/bin/env python
import numpy as np
from ROOT import TFile, TChain
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit


def main():
    """
    # branch: Board                     6737
    # branch: Channel                 184333  0,1,2,3
    # branch: Timestamp              2389745
    # branch: Energy                  743518
    # branch: EnergyShort             725173
    # branch: Flags                    17481
    # compassF_run_1.root
    """
    # energy_1d()
    # energy_2d()
    # plot_waveforms()
    # psa_cut()
    calibrate_energy()

    # TODO:
    # - figure out parameters of Energy and EnergyShort
    # - exact calculation of livetime for each file


def energy_1d():
    # TODO: make plots to compare weekend 1 and weekend 2 data
    # (be sure to edit filenames)

    # runList = [1,3,4,5] # weekend 1
    # runList = [1,3]
    runList = [13] # weekend 2
    weekend = 2

    fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox"
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))

    plt.figure(figsize=(10,6),facecolor='w')

    for chan in range(4):
    # for chan in range(1):
        print("Plotting channel",chan,"...")

        n = ch.Draw("Energy","Channel==%d" % (chan),"goff")
        ene = ch.GetV1()
        ene = np.asarray([ene[i] for i in range(n)])

        # full energy range
        hEne, xEne = np.histogram(ene, bins=1000, range=(0,20000))

        plt.cla() # clear plot from last time
        plt.semilogy(xEne[1:], hEne, ls='steps', c='r', label="Channel %d" % chan)
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        plt.tight_layout()
        # plt.show()
        plt.savefig("./plots/energy_1d_ch%d_week%d.pdf" % (chan, weekend))

        # zoom in on low e region
        hEne2, xEne2 = np.histogram(ene, bins=100, range=(0,2000))

        plt.cla() # clear plot from last time
        plt.semilogy(xEne2[1:], hEne2, ls='steps', c='r', label="Channel %d" % chan)
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        plt.legend(loc='best')
        plt.tight_layout()
        # plt.show()
        plt.savefig("./plots/energy_1d_ch%d_zoom.pdf" % (chan))


def energy_2d():

    # runList = [1,3,4,5] # weekend 1
    runList = [13] # weekend 2
    # runList = [1,3]

    fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox"
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))

    for chan in range(4):

        n = ch.Draw("Energy:EnergyShort","Channel==%d" % (chan),"goff")
        ene, eshort = ch.GetV1(), ch.GetV2()
        ene = np.asarray([ene[i] for i in range(n)])
        eshort = np.asarray([eshort[i] for i in range(n)])

        fig = plt.figure(figsize=(10,6),facecolor='w')
        plt.hist2d(ene, eshort, bins=(1000,1000), range=((0,8000),(0,8000)), norm=LogNorm())
        plt.xlabel("Energy")
        plt.ylabel("EnergyShort")
        plt.tight_layout()
        # plt.show()
        plt.savefig("./plots/energy_2d_ch%d.pdf" % (chan))

        # fig = plt.figure(figsize=(10,6),facecolor='w')
        # plt.hist2d(ene, eshort, bins=(1000,1000), range=((0,1500),(0,1500)), norm=LogNorm())
        # plt.xlabel("Energy")
        # plt.ylabel("EnergyShort")
        # plt.tight_layout()
        # # plt.show()
        # plt.savefig("./plots/energy_2d_ch%d_zoom.pdf" % (chan))


def plot_waveforms():

    # run 15 unfiltered is the waveform data in csv format
    # TIMETAG;ENERGY;ENERGYSHORT;FLAGS;SAMPLES

    chan = 0
    with open("./run_15/UNFILTERED/%d@DT5730 #2-11-1463_Data_run_15.csv" % chan) as f:
        raw = f.readlines()

    for i, line in enumerate(raw[2:3]):

        data = line.rstrip().split(";")

        timestamp = data[0]
        energy = data[1]
        energy_short = data[2]
        flags = data[3]
        wf = data[4:]

        print("Event number",i)
        print(timestamp)
        print(energy)
        print(energy_short)
        print(flags)
        print("wf is %d samples long" % len(wf))
        print("")

        ts_wf = np.arange(len(wf))

        fig = plt.figure(figsize=(10,6),facecolor='w')
        plt.plot(ts_wf, wf, "-r", label="Channel %d" % chan)
        plt.xlabel("Time [arb]")
        plt.ylabel("ADC value")
        plt.legend(loc='best')
        plt.show()

        exit()


def psa_cut(ene = None, eshort = None, chan = 1):

    def linear(x, m, b):
        return m*x + b

    def poly(x, a, b, c):
        return a*x**2 + b*x + c

    if ene == None:
        # runList = [1,3,4,5] # weekend 1
        runList = [13] # weekend 2
        # runList = [1,3]

        fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox/data"
        ch = TChain("Data_F")
        for run in runList:
        #run = 1
            fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
            ch.Add(fName)
        print("Found %d entries" % (ch.GetEntries()))

        n = ch.Draw("Energy:EnergyShort","Channel==%d" % (chan),"goff")
        ene, eshort = ch.GetV1(), ch.GetV2()
        ene = np.asarray([ene[i] for i in range(n)])
        eshort = np.asarray([eshort[i] for i in range(n)])

    fit_lo, fit_hi = 1000, 2500 # fit line where we don't have alphas in the spectrum

    init_m = (1421-735)/1000.
    print("init_m",init_m)

    idx = np.where((ene > fit_lo) & (ene < fit_hi) & (eshort > 10)) # (to cut out the noise)
    par, pcov = curve_fit(linear, ene[idx], eshort[idx]) #, p0=[init_m,0,1])
    print("Fit done.")
    print("Pars",par)
    print("Covariance",pcov)

    fig = plt.figure(figsize=(10,6),facecolor='w')
    plt.hist2d(ene, eshort, bins=(1000,1000), range=((0,8000),(0,8000)), norm=LogNorm())

    fit_hi = 7000
    buf = 120 # offset of linear cut

    xf = np.arange(fit_lo, fit_hi, 0.1)
    plt.plot(xf, linear(xf, par[0] - 0.01, par[1] + buf), '-r')
    # plt.plot(xf, poly(xf, *par), '-r')

    # Cut everything below our fit, and then some
    alphas = np.where(eshort > par[0] * ene + par[1] + buf)
    gammas = np.where(eshort < par[0] * ene + par[1] + buf)

    fig2 = plt.figure(figsize=(10,6),facecolor='w')
    plt.hist2d(ene[alphas], eshort[alphas], bins=(1000,1000), range=((0,8000),(0,8000)), norm=LogNorm())

    hEne, xEne = np.histogram(ene, bins=1000, range=(0, 20000))
    hEneG, xEneG = np.histogram(ene[gammas], bins=1000, range=(0, 20000))
    hEneA, xEneA = np.histogram(ene[alphas], bins=1000, range=(0, 20000))

    # Plot the histograms for the alpha and gamma events
    fig3 = plt.figure(figsize=(10,6),facecolor='w')

    plt.semilogy(xEne[1:], hEne, ls='steps', c='g', label="Totals for Channel %d" % chan)
    plt.semilogy(xEneG[1:], hEneG, ls='steps', c='r', label="Gamma events")
    plt.semilogy(xEneA[1:], hEneA, ls='steps', c='b', label = "Alpha events")
    plt.xlabel("Energy")
    plt.ylabel("Counts")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()


def calibrate_energy():
    e_peak = 1460.0  # Could be real, it feels right
    runList = [13]
    chan = 3

    # fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox/data"
    fileDir = "./data"

    det_mass = 7.698  # kg

    runtime = 0
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
        runtime += get_runtime(fName)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))
    print("Total runtime (hrs): {:.2f}".format(runtime))

    n = ch.Draw("Energy:EnergyShort", "Channel==%d" % (chan), "goff")
    ene, eshort = ch.GetV1(), ch.GetV2()
    ene = np.asarray([ene[i] for i in range(n)])
    eshort = np.asarray([eshort[i] for i in range(n)])

    idx = np.where((ene > 10) & (eshort > 10))

    fig = plt.figure(figsize=(10, 6), facecolor='w')

    hSene, xSene = np.histogram(ene[idx], bins=2000, range=(0, 20000))

    # ctTotal = sum(hSene[idx1:idx2])

    mx = xSene[np.argmax(hSene)]
    scale = e_peak / mx

    xCal = xSene * scale
    hCal = hSene / (det_mass * runtime)

    print("Found calibration constant for bin {}: {:.2f}".format(mx, scale))

    plt.semilogy(xCal[1:], hCal, ls='steps', c='r', label="Totals for Channel %d" % chan)
    plt.axvline(e_peak, c='b', label="40K: {:.2f}".format(e_peak))
    plt.legend(loc='best')
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (arb)")
    # plt.tight_layout()
    plt.show()

    # Perform a PSA cut on our data
    psa_cut(ene, eshort, chan)


if __name__=="__main__":
    main()
