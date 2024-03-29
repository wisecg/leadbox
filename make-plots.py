#!/usr/bin/env python
import numpy as np
from scipy.integrate import quad
from ROOT import TFile, TChain, TTree
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from array import array
import glob


def main():
    """
    branch: Board                     6737
	branch: Channel                 184333  0,1,2,3
	branch: Timestamp              2389745
	branch: Energy                  743518
	branch: EnergyShort             725173
	branch: Flags                    17481
	compassF_run_1.root

    calibrate_energy will calibrate, perform the PSA cut, and fit the alphas. It is likely
	all you will need. Method headers note graphs generated.
    """
    fig = plt.figure(figsize=(10, 6), facecolor='w')

    # skim_data(0)
    # energy_1d()
    # energy_2d()
    # plot_waveforms()
    # psa_cut()
    # check_rate()
    # skim_1d()

    calibrate_energy(2, [33])


def get_runtime(file_name, verbose=False):

    tf = TFile(file_name)
    tt = tf.Get("Data_F")
    n_ent = tt.GetEntries()
    ts = array('l', [0])
    tt.SetBranchAddress("Timestamp", ts)
    tt.GetEntry(0)
    ts_start = ts[0] / 1e12  # caen's timestamp is evidently in picoseconds
    tt.GetEntry(n_ent - 1)
    ts_stop = ts[0] / 1e12
    runtime = (ts_stop - ts_start) / 3600  # sec to hours
    if verbose:
        print("Start:{}  Stop:{}  Runtime (hrs):{}".format(
            ts_start, ts_stop, runtime))
    return runtime


def skim_1d():
    # checking the data rate for the new runs

    # runList = [1,3,4,5] # weekend 1
    # runList = [1,3]
    # runList = [9,13] # weekend 2
    # runList = [27] #weekend 3
    runList = [33]  #weekend 4
    weekend = 4

    runtime = 0
    fileDir = "./skim"
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d.root" % (fileDir, run)
        runtime += get_runtime(fName)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))
    print("Total runtime (hrs): {:.2f}".format(runtime))

    plt.figure(figsize=(10, 6), facecolor='w')

    for chan in range(4):
        # for chan in range(1,4):
        print("Plotting channel", chan, "...")

        ch.SetEstimate(ch.GetEntries() + 1)
        n = ch.Draw("Energy", "Channel==%d" % (chan), "goff")
        ene = ch.GetV1()
        ene = np.asarray([ene[i] for i in range(n)])

        # full energy range
        hEne, xEne = np.histogram(ene, bins=1000, range=(0, 20000))

        plt.cla()  # clear plot from last time
        plt.semilogy(
            xEne[1:], hEne, ls='steps', c='r', label="Channel %d" % chan)
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        plt.tight_layout()
        # plt.show()
        plt.savefig(
            "./plots/skim_plots/energy_1d_ch%d_week%d.pdf" % (chan, weekend))
        plt.clf()


def skim_data(weekend):
    rdata = [[1, 3, 4, 5], [9, 13], [27], [33], [34]]
    runList = rdata[weekend]

    thresh = {}
    thresh[1] = {0: 80, 1: 60, 2: 70, 3: 80}
    thresh[1] = {0: 80, 1: 60, 2: 70, 3: 80}
    thresh[3] = {0: 80, 1: 60, 2: 70, 3: 80}
    thresh[4] = {0: 80, 1: 60, 2: 70, 3: 80}
    thresh[5] = {0: 80, 1: 60, 2: 70, 3: 80}
    thresh[9] = {0: 80, 1: 70, 2: 100, 3: 90}
    thresh[13] = {0: 80, 1: 70, 2: 100, 3: 90}
    thresh[27] = {0: 80, 1: 75, 2: 70, 3: 90}
    thresh[33] = {0: 80, 1: 230, 2: 300, 3: 430}
    thresh[34] = {0: 80, 1: 230, 2: 300, 3: 430}

    # print(thresh)
    # for key in thresh:
    #     print(key, thresh[key])

    fileDir, skimDir = "./data", "./skim"

    for run in runList:
        filelist = glob.glob("{}/run_{}/FILTERED/compassF_run_{}".format(
            fileDir, run, run) + "*.root")
        print(filelist)
        for f in sorted(filelist):

            fcut = "(Channel==0 && Energy > {}) || ".format(thresh[run][0])
            fcut += "(Channel==1 && Energy > {}) || ".format(thresh[run][1])
            fcut += "(Channel==2 && Energy > {}) || ".format(thresh[run][2])
            fcut += "(Channel==3 && Energy > {})".format(thresh[run][3])

            ch = TChain("Data_F")
            ch.Add(f)
            print("Found %d entries" % (ch.GetEntries()))

            outName = "{}/run_{}.root".format(skimDir, run)
            outFile = TFile(outName, "RECREATE")
            outTree = TTree()
            outTree = ch.CopyTree(fcut)
            outTree.Write()
            outFile.Close()

            f2 = TFile(outName)
            tt2 = f2.Get("Data_F")
            print(tt2.GetEntries())


def check_rate():

    runList = [33]  # weekend 4
    weekend = 4

    runtime = 0
    fileDir = "./data"
    ch = TChain("Data_F")
    for run in runList:
        filelist = glob.glob("%s/run_%d/FILTERED/compassF_run_33" %
                             (fileDir, run) + "*.root")
        for f in filelist[:2]:
            runtime += get_runtime(f)
            ch.Add(f)
    print("Found %d entries" % (ch.GetEntries()))
    print("Total runtime (hrs): {:.2f}".format(runtime))

    # ch.SetEstimate(ch.GetEntries() + 1)
    # n = ch.Draw("Energy:Channel","","goff")
    # ene = ch.GetV1()
    # chan = ch.GetV2()
    #
    # ene = np.asarray([ene[i] for i in range(n)])
    # chan = np.asarray([chan[i] for i in range(n)])

    # print(len(ene))
    # idx = np.where(chan == 1)
    # print(type(idx),idx)
    # print(len(ene[idx]))

    # cts1 = len(ene[np.where(chan == 0)])
    # cts2 = len(ene[np.where(chan == 1)])
    # cts3 = len(ene[np.where(chan == 2)])
    # cts4 = len(ene[np.where(chan == 3)])
    #
    # for i in range(0,5):
    #     print("cts", i, len(ene[np.where(chan == i)]))

    # text_file = open("./counts.txt", "a")
    # text_file.write("cts1: %s\n" % cts1)
    # text_file.write("cts2: %s\n" % cts2)
    # text_file.write("cts3: %s\n" % cts3)
    # text_file.write("cts4: %s\n" % cts4)
    # text_file.close()

    for chan in range(1, 4):
        # for chan in range(1,2):
        print("Plotting channel", chan, "...")

        ch.SetEstimate(ch.GetEntries() + 1)
        n = ch.Draw("Energy", "Channel==%d" % (chan), "goff")
        ene = ch.GetV1()
        ene = np.asarray([ene[i] for i in range(n)])

        #energy cut here
        # ene = ene[ene > 100]

        # full energy range
        hEne, xEne = np.histogram(ene, bins=1000, range=(0, 20000))

        plt.cla()  # clear plot from last time
        plt.semilogy(
            xEne[1:], hEne, ls='steps', c='r', label="Channel %d" % chan)
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        plt.tight_layout()
        plt.show()
        # plt.savefig("./plots/energy_1d_ch%d_week%d.pdf" % (chan, weekend))
        plt.clf()


def energy_1d():
    # TODO: make plots to compare weekend 1 and weekend 2 data
    # (be sure to edit filenames)

    # runList = [1,3,4,5] # weekend 1
    # runList = [1,3]
    runList = [34]  # weekend 2
    # runList = [27] #weekend 3
    weekend = 4

    runtime = 0
    fileDir = "./data"
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
        runtime += get_runtime(fName)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))
    print("Total runtime (hrs): {:.2f}".format(runtime))

    plt.figure(figsize=(10, 6), facecolor='w')

    for chan in range(4):
        # for chan in range(1,4):
        print("Plotting channel", chan, "...")

        ch.SetEstimate(ch.GetEntries() + 1)
        n = ch.Draw("Energy", "Channel==%d" % (chan), "goff")
        ene = ch.GetV1()
        ene = np.asarray([ene[i] for i in range(n)])

        # full energy range
        hEne, xEne = np.histogram(ene, bins=1000, range=(0, 20000))

        plt.cla()  # clear plot from last time
        plt.semilogy(
            xEne[1:], hEne, ls='steps', c='r', label="Channel %d" % chan)
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        plt.tight_layout()
        plt.show()
        # plt.savefig("./plots/energy_1d_ch%d_week%d.pdf" % (chan, weekend))
        plt.clf()

        # zoom in on low e region
        hEne2, xEne2 = np.histogram(ene, bins=100, range=(0, 2000))

        plt.cla()  # clear plot from last time
        plt.semilogy(
            xEne2[1:], hEne2, ls='steps', c='r', label="Channel %d" % chan)
        plt.xlabel("Energy")
        plt.ylabel("Counts")
        plt.legend(loc='best')
        plt.tight_layout()
        plt.show()
        # plt.savefig("./plots/energy_1d_ch%d_zoom.pdf" % (chan))
        plt.clf()


def energy_2d():

    # runList = [1,3,4,5] # weekend 1
    runList = [13]  # weekend 2
    # runList = [1,3]

    fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox"
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))

    for chan in range(4):

        n = ch.Draw("Energy:EnergyShort", "Channel==%d" % (chan), "goff")
        ene, eshort = ch.GetV1(), ch.GetV2()
        ene = np.asarray([ene[i] for i in range(n)])
        eshort = np.asarray([eshort[i] for i in range(n)])

        fig = plt.figure(figsize=(10, 6), facecolor='w')
        # plt.hist2d(ene, eshort, bins=(1000,1000), range=((0,8000),(0,8000)), norm=LogNorm())
        plt.hist2d(
            ene,
            eshort,
            bins=(1000, 1000),
            range=((0, 1500), (0, 1500)),
            norm=LogNorm())
        plt.xlabel("Energy")
        plt.ylabel("EnergyShort")
        plt.tight_layout()
        # plt.show()
        plt.savefig("./plots/energy_2d_ch%d.pdf" % (chan))
        plt.clf()


def plot_waveforms():

    # run 15 unfiltered is the waveform data in csv format
    # TIMETAG;ENERGY;ENERGYSHORT;FLAGS;SAMPLES

    chan = 0
    with open("./run_15/UNFILTERED/%d@DT5730 #2-11-1463_Data_run_15.csv" %
              chan) as f:
        raw = f.readlines()

    for i, line in enumerate(raw[2:3]):

        data = line.rstrip().split(";")

        timestamp = data[0]
        energy = data[1]
        energy_short = data[2]
        flags = data[3]
        wf = data[4:]

        print("Event number", i)
        print(timestamp)
        print(energy)
        print(energy_short)
        print(flags)
        print("wf is %d samples long" % len(wf))
        print("")

        ts_wf = np.arange(len(wf))

        fig = plt.figure(figsize=(10, 6), facecolor='w')
        plt.plot(ts_wf, wf, "-r", label="Channel %d" % chan)
        plt.xlabel("Time [arb]")
        plt.ylabel("ADC value")
        plt.legend(loc='best')
        # plt.show()
        plt.savefig("./plots/wf_test.pdf")
        plt.clf()

        exit()


def calibrate_energy(chan, runList):
    # Calibrates the spectum of channel "chan" across all runs in runList. Generates a 1D energy
    # spectrum with a line at the K40 peak. Then calls psa_cut.

    e_peak = 1460.0  # K40 energy

    # fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox/data"
    fileDir = "./skim"

    det_mass = 7.698  # kg

    runtime = 0
    ch = TChain("Data_F")
    for run in runList:
        fName = "%s/run_%d.root" % (fileDir, run)
        runtime += get_runtime(fName)
        ch.Add(fName)
    print("Found %d entries" % (ch.GetEntries()))
    print("Total runtime (hrs): {:.2f}".format(runtime))

    n = ch.Draw("Energy:EnergyShort", "Channel==%d" % (chan), "goff")
    ene, eshort = ch.GetV1(), ch.GetV2()
    ene = np.asarray([ene[i] for i in range(n)])
    eshort = np.asarray([eshort[i] for i in range(n)])

    idx = np.where((ene > 10) & (eshort > 10))
    print("Total events: {}  Num. kept after eshort=0: {}".format(
        len(ene), len(ene[idx])))

    # diagnostic plot
    # plt.hist(ene, bins=(1000), range=(0,10000), histtype='step', color='b', label="no cut")
    # plt.hist(ene[idx], bins=1000, range=(0, 10000), histtype='step', color='r', label="with cut")
    # plt.legend()
    # plt.show()

    hSene, xSene = np.histogram(ene[idx], bins=2000, range=(0, 20000))

    # ctTotal = sum(hSene[idx1:idx2])

    # # placeholder
    # params = {
    #     "run or channel":[guess values, ],
    #     "run 2":[guess 2]
    # }

    mx = xSene[np.argmax(hSene)]
    scale = e_peak / mx
    hscale = 1 / (det_mass * runtime)

    xCal = xSene * scale
    hCal = hSene / (det_mass * runtime)

    print("Found calibration constant for bin {}: {:.2f}".format(mx, scale))

    plt.semilogy(
        xCal[1:], hCal, ls='steps', c='r', label="Totals for Channel %d" % chan)
    plt.axvline(e_peak, c='b', label="40K: {:.2f}".format(e_peak))
    plt.legend(loc='best')
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts (arb)")
    # plt.tight_layout()
    plt.savefig("./plots/ecal_spectrum.png")
    plt.show()
    plt.clf()

    # Perform a PSA cut on our data
    psa_cut(ene, eshort, chan, runList, scale, hscale)


def psa_cut(ene=None, eshort=None, chan=1, runList=[13], scale=1, hscale=1):
    # Fit our alphas

    # Generates a 2D graph of energy vs. eshort with our cut line drawn in, another graph of
    # the same after performing our cut, and then ploits the separate histograms for the
    # alpha and gamma events on top of one another. If called from calibrate_energy, all
    # energies will be scaled from the K40 peak.

    def linear(x, m, b):
        return m * x + b

    # def poly(x, a, b, c):
    # return a*x**2 + b*x + c

    # if ene == None:
    #     print("i'm here!")
    #     # runList = [1,3,4,5] # weekend 1
    #     runList = [13] # weekend 2
    #     # runList = [1,3]
    #
    #     fileDir = "/Users/ccenpa/Desktop/coherent/Analysis/leadbox/data"
    #     ch = TChain("Data_F")
    #     for run in runList:
    #     #run = 1
    #         fName = "%s/run_%d/FILTERED/compassF_run_%d.root" % (fileDir, run, run)
    #         ch.Add(fName)
    #     print("Found %d entries" % (ch.GetEntries()))
    #
    #     n = ch.Draw("Energy:EnergyShort","Channel==%d" % (chan),"goff")
    #     ene, eshort = ch.GetV1(), ch.GetV2()
    #     ene = np.asarray([ene[i] for i in range(n)])
    #     eshort = np.asarray([eshort[i] for i in range(n)])

    # Dictionary of fit bounds and offsets for the various 'good' datasets. Runs not in here
    # are bad, and don't behave nicely; if a run is in here, we use this dict for cleaner cuts.
    fit_dict = {
        "1, [13]": (1000, 2500, 160),
        "3, [13]": (1000, 2500, 100),
        "2, [27]": (2500, 3000, 60),
        "(1, [33])": (1000, 2500, 180)
    }

    idx_noisecut = np.where((ene > 10) & (eshort > 10))
    ene = ene[idx_noisecut]
    eshort = eshort[idx_noisecut]

    fits = []
    ind = str(chan) + ", " + str(runList)
    if ind in fit_dict:
        fits = fit_dict[str(chan) + ", " + str(runList)]
    else:
        fits = (
            1000, 2500, 100
        )  # Just use any fit params, the data not in the dict ==> bad data
    fit_lo, fit_hi = fits[0], fits[
        1]  # fit line where we don't have alphas in the spectrum

    init_m = (1421 - 735) / 1000.
    print("init_m", init_m)
    idx = np.where((ene > fit_lo) & (ene < fit_hi) &
                   (eshort > 10))  # (to cut out the noise)
    par, pcov = curve_fit(linear, ene[idx], eshort[idx])
    print("Fit done.")
    print("Pars", par)
    print("Covariance", pcov)

    plt.cla()
    # fig = plt.figure(figsize=(10,6),facecolor='w')
    plt.hist2d(
        ene,
        eshort,
        bins=(1000, 1000),
        range=((0, 8000), (0, 8000)),
        norm=LogNorm())

    fit_hi = 7000
    fit_lo = 1000
    buf = fits[2]  # offset of linear cut

    xf = np.arange(fit_lo, fit_hi, 0.1)
    plt.plot(xf, linear(xf, par[0], par[1] + buf), '-r')
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("E-short (keV)", ha='right', y=1)
    plt.ylim(0, 4000)  # Note this is the zoomed version!!!
    plt.xlim(2000, 6000)  # Also zoomed
    plt.savefig("./plots/psa2d_id{},{}.pdf".format(
        chan, runList))  # 2D hist of e vs eshort

    ene = ene[~np.isnan(ene)]
    eshort = eshort[~np.isnan(eshort)]

    # Cut everything below our fit, and then some
    alphas = np.where((eshort > par[0] * ene + par[1] + buf))
    gammas = np.where((eshort < par[0] * ene + par[1] + buf))

    plt.clf()
    plt.hist2d(
        ene[alphas],
        eshort[alphas],
        bins=(1000, 1000),
        range=((0, 8000), (0, 8000)),
        norm=LogNorm())
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.savefig("./plots/psa_2d_id{},{}ac.pdf".format(
        chan, runList))  # 2D hist after cutting
    plt.show()

    hEne, xEne = np.histogram(ene, bins=1000, range=(0, 20000))
    hEneG, xEneG = np.histogram(ene[gammas], bins=1000, range=(0, 20000))
    hEneA, xEneA = np.histogram(ene[alphas], bins=1000, range=(0, 20000))

    # Plot the histograms for the alpha and gamma events
    plt.clf()
    plt.semilogy(
        xEne[1:] * scale,
        hEne * hscale,
        ls='steps',
        c='g',
        label="Totals for Channel %d" % chan)
    plt.semilogy(
        xEneG[1:] * scale,
        hEneG * hscale,
        ls='steps',
        c='r',
        label="Gamma events")
    plt.semilogy(
        xEneA[1:] * scale,
        hEneA * hscale,
        ls='steps',
        c='b',
        label="Alpha events")
    plt.xlabel("Energy")
    plt.ylabel("Counts/kg-hr")
    plt.legend(loc='best')
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/psa_cut.png")

    plt.savefig("./plots/agspec {},{}.pdf".format(
        chan, runList))  # alpha vs gamma hists
    plt.show()

    fit_alphas(hEneA * hscale, xEneA * scale)


def fit_alphas(ha, xa):
    # Actually performs the alpha fit. Creates a graph of the fitted alphas.

    def gauss(x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params) - 1, 3):
            x0 = params[i]
            a = params[i + 1]
            sigma = params[i + 2]
            y += a * np.exp(-(x - x0)**2 / (2 * sigma**2))
        y = y + params[-1]
        return y

    # # pk_list = [2798, 3212]
    # pk_list = [3212]
    # p0_list = []
    # bnd = [[],[]]
    # for pk in pk_list:
    #     p0_list.extend([pk, .6, 100])
    #     # bnd[0].extend([pk*.9, ])
    #     # bnd[1].extend([pk*1.1])
    # p0_list.append(.005)
    # # bnd = tuple(bnd)

    p0_list = [
        2798, .5, 50, 3212, .68, 100, 3688, .46, 75, 4256, .058, 50, 4900, 1,
        100, .005
    ]

    # 0, np.inf, -np.inf, 100000

    bnd = ((p0_list[0] * .9, .4, 0, p0_list[3] * .9, .65, 0, p0_list[6] * .9,
            .4, 0, p0_list[9] * .9, .04, 0, p0_list[12] * .9, .3, 50, 0),
           (p0_list[0] * 1.1, 2, 500, p0_list[3] * 1.1, 2, 500,
            p0_list[6] * 1.1, 2, 500, p0_list[9] * 1.1, 2, 500,
            p0_list[12] * 1.1, 2, 500, .02))

    # bnd = ((p0_list[0]*.9, .4, 0, p0_list[3]*.9, .65, 0, p0_list[6]*.9, .4, 0, p0_list[9]*.9, .04, 0, p0_list[12]*.9, .8, 100, 0),
    #         (p0_list[0]*1.1, 2, 500, p0_list[3]*1.1, 2, 500, p0_list[6]*1.1, 2, 500, p0_list[9]*1.1, 2, 500, p0_list[12]*1.005, 1.25, 500, .02))

    # bnd = ((p0_list[0]*.5, .6, 0, 0),
    #         (p0_list[0]*1.4, .7, 400, .02))

    # 100,320

    par, pcov = curve_fit(
        gauss, xa[100:320], ha[100:320], p0=p0_list, bounds=bnd)
    print(par)

    # text_file = open("./params.txt", "a")
    # text_file.write("{}".format(par))
    # text_file.close()

    counts = 0
    n = 1

    for i in range(0, len(par) - 1, 3):
        mu, amp, sig, bkg = par[i], par[i + 1], par[i + 2], par[-1]
        print("Scanning peak ", n, " at energy", mu)
        ans = quad(gauss, mu - 5 * sig, mu + 5 * sig, args=(mu, amp, sig, bkg))
        # print(ans[0], " counts/(kg*hr)")
        answer = ans[0]
        print("R = {:.4f} +/- {:.2e} cts/kg-hr-keV".format(answer, ans[1]))
        n += 1

    plt.cla()
    xf = np.arange(0, 10000, 0.1)

    plt.plot(xf, gauss(xf, *par), '-r')
    plt.plot(xa[100:320], ha[100:320], ls='steps', c='b', label="Alpha events")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Cts/kg-hr-keV")

    plt.show()


if __name__ == "__main__":
    main()
