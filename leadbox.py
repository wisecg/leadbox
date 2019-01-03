#!/usr/bin/env python3
import os
import json
import glob
import numpy as np
from scipy.optimize import curve_fit
from pprint import pprint
from array import array
from ROOT import TFile, TTree, TChain
import matplotlib.pyplot as plt
plt.style.use("clint.mplstyle")
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle


def main():
    """
    NOTE:
    here are the CAEN DT1502 digitizer's:
    --> filename convention: compassF_run_1.root
    --> "filtered" ROOT data format:
    branch: Board                     6737
	branch: Channel                 184333  [0,1,2,3]
	branch: Timestamp              2389745
	branch: Energy                  743518
	branch: EnergyShort             725173
	branch: Flags                    17481
    """
    # load our metadata
    global run_info
    with open('run_info.json') as f:
        run_info = json.load(f)

    # -- analysis routines --
    # plot_waveforms()
    # get_thresholds()
    # skim_data("1")
    # calibrate_channel("4", 1)
    # psd_channel("2",3)
    # alpha_spec("2", 1)
    plot_thresholds()

    # -- batch mode --
    # process()


def process(debug=False):
    """
    wrapper function to process multiple runs,
    comment lines in/out as necessary
    """

    # -- create skim files --
    # for wk in run_info["goodres_chans"]:
        # skim_data(wk)

    # -- analyze spectra --
    if debug:
        data_list = [("1",1),]

    else:
        data_list = []
        for wk in run_info["goodres_chans"]:
            for ch in run_info["goodres_chans"][wk]:
                data_list.append((wk, ch))

    for wk, ch in data_list:
        # calibrate_channel(wk, ch)
        # psd_channel(wk, ch)
        alpha_spec(wk, ch)


def plot_waveforms():
    """
    in run 15, we saved the waveform data in csv format:
    TIMETAG;ENERGY;ENERGYSHORT;FLAGS;SAMPLES[...]

    note: typical PMT rise time is 20--50 ns
    """
    chan = 1
    wf_file = "./data/run_15/UNFILTERED/"\
              "{}@DT5730 #2-11-1463_Data_run_15.csv".format(chan)

    with open(wf_file) as f:
        raw = f.readlines()

    # for i, line in enumerate(raw[50:51]):
    for i, line in enumerate(raw[1:100]):

        data = line.rstrip().split(";")
        time, ene, eshort, flags = data[0], data[1], data[2], data[3]

        if int(ene) < 300:
            continue

        wf = np.asarray([int(d) for d in data[4:]]) * -1
        ts = np.arange(len(wf)) * 4 # sample rate: 250 MHz --> 4 ns period
        bl = np.mean(wf[:20])
        wf = wf - bl # baseline subtract

        trap_wf = trap_filter(wf, rt=40, ft=20) # long trap for energy
        trap_short = trap_filter(wf, rt=20, ft=1) # short trap for psd

        # -- plot 1: lotsa waveforms --
        # plt.plot(ts, wf, "-b", lw=2, alpha=0.5)
        # plt.plot(ts, trap_wf, "-r", lw=2, alpha=0.5)
        # plt.plot(ts, trap_short, "-m", alpha=0.5)

    # -- plot 2: illustrate long and short trap --
    plt.plot(ts, wf , "-b", lw=2, label="week 2, channel {}   ".format(chan))
    plt.plot(ts, trap_wf, "-r", lw=2, label="40--20--40 energy trap")
    plt.plot(ts, trap_short, "-g", label="20--1--20 short trap")

    # print(len(wf))

    plt.xlabel("Time [ns]", ha='right', x=1)
    plt.ylabel("ADC", ha='right', y=1)
    plt.ylim(-500, np.max(wf) * 1.1)
    plt.legend()
    plt.tight_layout()
    # plt.show()

    # plt.savefig("./plots/trap_wfs.pdf") # plot 1
    plt.savefig("./plots/trap_demo.pdf") # plot 2


def trap_filter(wfs, rt=400, ft=200, dt=0, test=False):
    """
    assumes `wfs` is an ndarray (vectorizable),
    where each row is a waveform: with shape (nwfs, nsamp)
    in the case of a single waveform w/ shape (nsamp,),
    return a trapezoid filter w/ the same shape.
    """

    flip = False
    if len(wfs.shape) == 1:
        flip = True
        wfs = wfs.reshape(1, len(wfs))

    nwfs, nsamp = wfs.shape[0], wfs.shape[1]

    wfs_minus_ramp = np.zeros_like(wfs)
    wfs_minus_ramp[:, :rt] = 0
    wfs_minus_ramp[:, rt:] = wfs[:, :nsamp - rt]

    wfs_minus_ft_and_ramp = np.zeros_like(wfs)
    wfs_minus_ft_and_ramp[:, :(ft + rt)] = 0
    wfs_minus_ft_and_ramp[:, (ft + rt):] = wfs[:, :nsamp - ft - rt]

    wfs_minus_ft_and_2ramp = np.zeros_like(wfs)
    wfs_minus_ft_and_2ramp[:, :(ft + 2 * rt)] = 0
    wfs_minus_ft_and_2ramp[:, (ft + 2 * rt):] = wfs[:, :nsamp - ft - 2 * rt]

    scratch = wfs - (wfs_minus_ramp + wfs_minus_ft_and_ramp + wfs_minus_ft_and_2ramp)

    trap_wfs = np.zeros_like(wfs)
    trap_wfs = np.cumsum(trap_wfs + scratch, axis=1) / rt

    if test:
        # diagnostic plot
        iwf = 0
        plt.plot(np.arange(nsamp), wfs[iwf], '-b', lw=2, label='raw wf')
        # plt.plot(np.arange(nsamp), wfs_minus_ramp[iwf], '-b', label='wf-ramp')
        # plt.plot(np.arange(nsamp), wfs_minus_ft_and_ramp[iwf], '-g', label='wf-ft-ramp')
        # plt.plot(np.arange(nsamp), wfs_minus_ft_and_2ramp[iwf], '-m', label='wf-ft-2ramp')
        # plt.plot(np.arange(nsamp), scratch[iwf], '-b', label='scratch')
        plt.plot(np.arange(nsamp), trap_wfs[iwf], '-r', lw=2, label='trap wf')

        plt.ylim(-200, 1.2*np.amax(wfs[iwf]))
        plt.legend()
        plt.show()
        exit()

    if flip:
        trap_wfs = trap_wfs.reshape(nsamp)

    return trap_wfs


def get_runtime(file_name, verbose=False):
    """ get the runtime by finding the first and last timestamps in
    the TTree we got from the CAEN digitizer
    """
    from ROOT import TFile, TTree

    tf = TFile(file_name)
    tt = tf.Get("Data_F")
    n_ent = tt.GetEntries()
    ts = array('l', [0])
    tt.SetBranchAddress("Timestamp", ts)
    tt.GetEntry(0)
    ts_start = ts[0] / 1e12 # caen's timestamp is evidently in picoseconds
    tt.GetEntry(n_ent - 1)
    ts_stop = ts[0] / 1e12
    runtime = (ts_stop - ts_start) / 3600  # sec to hours
    if verbose:
        print("Start:{}  Stop:{}  Runtime (hrs):{}".format(
            ts_start, ts_stop, runtime))
    return runtime


def get_thresholds():
    """
    manually set a reasonable software threshold for each run/channel spectrum
    then write the results in run_info.json
    """
    wk = "4"

    # load metadata
    runs = run_info["runs"][wk]

    # load data
    files = []
    for r in runs:
        # handle runs where we have multiple subfiles (_x.root)
        tmp = glob.glob("./data/run_{}/FILTERED/compassF_run_{}*.root".format(r,r))
        files.extend(tmp)

    print("Loading week", wk," ...")
    runtime = 0
    ch = TChain("Data_F")
    for f in files:
        runtime += get_runtime(f)
        ch.Add(f)

    # loop over all channels, not just good-res channels.  Assumes chans 0--3.
    for chan in range(4):
        print("Week {}, channel {}".format(wk, chan))

        ch.SetEstimate(ch.GetEntries() + 1)
        n = ch.Draw("Energy", "Channel==%d" % (chan), "goff")
        ene = ch.GetV1()
        ene = np.asarray([ene[i] for i in range(n)])

        hEne, xEne = np.histogram(ene, bins=1000, range=(0, 20000))

        plt.cla()

        plt.semilogy(
            xEne[1:], hEne, ls='steps', c='b', lw=2,
            label="Week {}, chan {}".format(wk, chan)
            )
        plt.xlabel("Energy", ha='right')
        plt.ylabel("Counts")
        plt.tight_layout()
        plt.show()

        # exit()
        # plt.savefig(
            # "./plots/skim_plots/energy_1d_ch%d_week%d.pdf" % (chan, weekend))

        print("Threshold? (enter an ADC value, or q to quit).")
        thresh = input()
        if "q" in thresh:
            exit()
        print("Good resolution?")
        goodres = input()
        print("thresh:",thresh, "goodres:",goodres)

        # could automatically write to run_info.json here if I wanted to


def skim_data(wk):
    """ run this automatically w/ process()
    apply energy thresholds to create reduced size "skim" ROOT files
    for faster plotting
    --> output files: "./data/skim/week_{}.root"
    """
    wk = "4"
    out_name = "./data/skim/week_{}.root".format(wk)

    # load metadata
    runs = run_info["runs"][wk]
    chans = run_info["goodres_chans"][wk]
    thresh = run_info["thresholds"][wk]

    # build a ROOT TCut string
    print("Skimming data from week", wk, "run list:", runs)
    cut = "("
    for chan in chans:
        cut += "(Channel=={} && Energy>{}) || ".format(chan,thresh[int(chan)])
    cut = cut[:-4]
    cut += ")"
    print("  Using cut:",cut)

    # load data
    files = []
    for r in runs:
        # handle runs where we have multiple subfiles (_x.root)
        tmp = glob.glob("./data/run_{}/FILTERED/compassF_run_{}*.root".format(r,r))
        files.extend(tmp)

    runtime = 0
    ch = TChain("Data_F")
    for f in sorted(files):
        runtime += get_runtime(f)
        ch.Add(f)
    print("  Data loaded.  Found",ch.GetEntries(),"entries.")
    print("  Runtime (hrs): {:.2f}".format(runtime))

    # create a TTree from the TCut and write to the output file
    print("Writing to output file:\n  ",out_name)
    out_file = TFile(out_name, "RECREATE")
    out_tree = TTree()
    out_tree = ch.CopyTree(cut)
    out_tree.Write()
    out_file.Close()

    # check it
    f2 = TFile(out_name)
    tt2 = f2.Get("Data_F")
    print("Done.  Skim file has",tt2.GetEntries(),"entries.")


def get_hist(np_arr, x_lo, x_hi, xpb, nb=None, shift=False, wts=None):
    """ wrapper function for numpy's histogram """
    if nb is None:
        nb = int((x_hi - x_lo) / xpb)
    y, x = np.histogram(np_arr, bins=nb, range=(x_lo, x_hi), weights=wts)
    y = np.insert(y, 0, 0, axis=0)
    if shift:
        x = x - xpb / 2.
    return x, y


def load_data(wk, chan, use_skim=True):
    """ get the appropriate TChain for a wk/chan """

    runs = run_info["runs"][wk]

    if use_skim:
        files = ["./data/skim/week_{}.root".format(wk)]
    else:
        files = []
        for r in runs:
            # handle runs where we have multiple subfiles (_x.root)
            tmp = glob.glob("./data/run_{}/FILTERED/compassF_run_{}*.root".format(r,r))
            files.extend(tmp)

    ch = TChain("Data_F")
    for f in files:
        ch.Add(f)

    print("\nWeek {}, channel {}, runs: {}, Entries: {}"
          .format(wk, chan, runs, ch.GetEntries()))

    return ch


def clean_data(ch, ene, eshort):
    """ sometimes eshort can have bad values, we want to flag that """
    idx = np.where(eshort < 1)
    if len(idx[0] > 0):
        nzero = len(idx[0])
        ntot = ch.GetEntries()
        print("  Warning, {} eshort==0 entries found ({:.2f}%)"
              .format(nzero, 100*nzero/ntot))
        idx = np.where(eshort > 1)
        ene, eshort = ene[idx], eshort[idx]
    return ene, eshort


def calibrate_channel(wk, chan, use_skim = True):
    """
    run this automatically w/ process()
    calibrate a run+channel, and determine its PSA parameters
    """

    # load metadata
    pk_ene = float(run_info["e_peak"])
    pk_guess = run_info["k40_peak"][wk][chan]
    runs = run_info["runs"][wk]
    runtime = float(run_info["runtime"][wk])

    # load data
    ch = load_data(wk, chan)

    # -- calibrate the channel to K40 peak --
    ch.SetEstimate(ch.GetEntries())
    n = ch.Draw("Energy","Channel=={}".format(chan),"goff")
    ene = ch.GetV1()
    ene = np.asarray([ene[i] for i in range(n)])

    x_lo, x_hi, xpb = 0, 20000, 10
    xE, hE = get_hist(ene, x_lo, x_hi, xpb)

    idx = np.where(xE > pk_guess-10)
    h_max = np.argmax(hE[idx])
    x_max = xE[idx][h_max]

    gain = pk_ene / x_max
    xE = xE * gain
    print("   gain: {:.4f}".format(gain))

    idx = np.where(xE < 30000)

    plt.semilogy(xE[idx], hE[idx], ls='steps', lw=2, c='b',
                 label="week {}, chan {}".format(wk, chan))

    plt.axvline(pk_ene, c='r', lw=2, alpha=0.6,
                label = "40K peak, gain {:.2f}".format(gain))

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # i could automatically update run_info.json here,
    # but let's just write it in by hand for now.


def psd_channel(wk, chan, use_skim = True):
    """
    run this automatically w/ process()
    determine a channel's PSA parameters
    """
    # load metadata
    runs = run_info["runs"][wk]
    gain = run_info["gain"][wk][chan]
    runtime = float(run_info["runtime"][wk])

    # load data
    ch = load_data(wk, chan)

    # load data into numpy arrays
    ch.SetEstimate(ch.GetEntries())
    n = ch.Draw("Energy:EnergyShort","Channel=={}".format(chan),"goff")
    ene, eshort = ch.GetV1(), ch.GetV2()
    ene = np.asarray([ene[i] for i in range(n)])
    eshort = np.asarray([eshort[i] for i in range(n)])

    # data cleaning step
    ene, eshort = clean_data(ch, ene, eshort)

    # apply gain and calculate epsd
    ene = ene * gain
    epsd = eshort / ene # kinda like A/E ...

    # -- plot 1: check 1d hist
    # x_lo, x_hi, xpb = 0, 25000, 10
    # nbx = int((x_hi-x_lo)/xpb)
    # # plt.hist(ene, nbx, range=(x_lo, x_hi), color='b', histtype='step',
    #          # label="hist")
    # xE, hE = get_hist(ene, x_lo, x_hi, xpb, shift=False)
    # plt.plot(xE, hE, c='r', lw=1, ls='steps',
    #          label='wk {}, chan {}'.format(wk, chan))
    # plt.xlabel("Energy (keV)", ha='right', x=1)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    # exit()
    # plt.cla()

    # -- plot 2: Energy vs E_psd

    # testing
    # x_lo, x_hi, xpb = 0, 8000, 10
    # y_lo, y_hi, ypb = 0, 1.5, 0.001

    x_lo, x_hi, xpb, y_lo, y_hi, ypb = run_info["psd_view"][wk][chan]

    nbx = int((x_hi-x_lo)/xpb)
    nby = int((y_hi-y_lo)/ypb)

    plt.clf()
    fig = plt.figure()
    plt.hist2d(ene, epsd,
               bins=[nbx, nby],
               range=((x_lo, x_hi), (y_lo, y_hi)),
               norm=LogNorm(),
               cmap='jet')

    # if we have PSD pars, draw the alpha box
    try:
        bx1, by1, bx2, by2 = run_info["alpha_box"][wk][chan]
        plt.gca().add_patch(
            Rectangle((bx1, by1), (bx2 - bx1), (by2 - by1),
                      lw=1, edgecolor='r', facecolor='none'))
    except:
        pass

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel(r"$\mathregular{E_{PSD}}$ (arb)", ha='right', y=1)
    plt.colorbar()
    plt.tight_layout()

    # plt.show()
    plt.savefig("./plots/alphas_wk{}_ch{}.png".format(wk, chan), dpi=300)


def multi_gauss(x, *pars):
    """
    multipeak gaussian function.
    pars should be = [mu1, a1, sig1, mu2, a2, sig2, ... , bkg]
    """
    y = np.zeros_like(x)

    for i in range(0, len(pars) - 1, 3):

        x0, a, sig = pars[i], pars[i+1], pars[i+2]

        y = y + a * np.exp(-(x - x0)**2 / (2 * sig**2))

    y = y + pars[-1]
    return y


def alpha_spec(wk, chan):
    """
    run this automatically w/ process()
    plot the alpha spectrum for each good channel (depends on successful PSD)
    """
    # load metadata
    runs = run_info["runs"][wk]
    gain = run_info["gain"][wk][chan]
    runtime = float(run_info["runtime"][wk])
    mass = float(run_info["det_mass_kg"])
    expo = (runtime/24) * mass # kg-d

    # load data
    ch = load_data(wk, chan)
    ch.SetEstimate(ch.GetEntries())
    n = ch.Draw("Energy:EnergyShort","Channel=={}".format(chan),"goff")
    ene, eshort = ch.GetV1(), ch.GetV2()
    ene = np.asarray([ene[i] for i in range(n)])
    eshort = np.asarray([eshort[i] for i in range(n)])

    # data cleaning step
    ene, eshort = clean_data(ch, ene, eshort)

    # apply gain and calculate epsd
    ene = ene * gain
    epsd = eshort / ene # kinda like T/E ...

    # -- select alpha events --
    bx1, by1, bx2, by2 = run_info["alpha_box"][wk][chan]
    ixa = np.where((ene > bx1) & (epsd > by1) & (ene < bx2) & (epsd < by2) )
    ene, eshort = ene[ixa], eshort[ixa]

    x_lo, x_hi, xpb = bx1, bx2, 10
    xE, hE = get_hist(ene, x_lo, x_hi, xpb, shift=False)

    # -- fit peaks --
    pars = [
        # alpha peaks (mu, amp, sig)
        2798, 100, 50,
        3212, 100, 100,
        3688, 100, 75,
        4256, 10, 50,
        5000, 100, 100,

        # flat background
        0.005
        ]

    # bounds should be = ((mu1lo, a1lo, ...), (mu1hi, a1hi, ...))
    bnds = [[],[]]
    for i in [0,3,6,9,12]:
        bnds[0].extend([pars[i] * .8, 0, 0])
        bnds[1].extend([pars[i] * 1.1, 10000, 500])
    bnds[0].append(0)
    bnds[1].append(1)

    # fit function x values
    xF = np.arange(x_lo, x_hi, xpb)

    # initial guess
    fGuess = multi_gauss(xF, *pars)

    # run fit and apply result
    par, pcov = curve_fit(multi_gauss, xE, hE, p0=pars, bounds=bnds)
    perr = np.sqrt(np.diag(pcov))
    fFit = multi_gauss(xF, *par)

    # -- update run_info.json with the results and errors --

    # paste these into run_info.json manually
    # print("""run_info["fit_results"][{}][{}] = """.format(wk, chan))
    # for p in par:
    #     print("{:.4e}".format(p) + ", ")

    # print("""run_info["fit_errors"][{}][{}] = """.format(wk, chan))
    # for p in par:
    #     print("{:.4e}".format(p) + ", ")

    # don't use this, it messes up the formatting of run_info.json
    # run_info["fit_results"] = {} <- also this wouldn't work in batch
    # run_info["fit_results"][wk][ch] = par
    # run_info["fit_errors"] = {}
    # run_info["fit_errors"][wk][ch] = perr
    # update_info()

    # -- calculate alpha rates --
    print("Week {}, chan {}, exposure {:.2f}".format(wk, chan, expo))

    # normalize counts to cts/kev/kg-d
    dE = xE[1] - xE[0]

    hE = hE / dE / expo
    fGuess = fGuess / dE / expo
    fFit = fFit / dE / expo

    # -- make plot --
    plt.clf()
    fig = plt.figure()

    plt.plot(xE, hE, c='b', lw=1, ls='steps',
             label='alphas, wk {}, chan {}'.format(wk, chan))

    plt.plot(xF, fGuess, '-g', alpha=0.6, label='guess')
    plt.plot(xF, fFit, '-r', label='fit')

    plt.ylabel(r"$\alpha$ Rate (cts / keV / kg-d)", ha='right', y=1)

    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/alpharate_wk{}_ch{}.pdf".format(wk,chan))


def update_info():
    """
    this *could* update run_info.json,
    but it messes up my nice formatting,
    so let's not use this for now.
    """
    global run_info
    with open("run_info.json", "wt") as f:

        # this just prints a big block of text
        # json.dump(run_info, f)

        # this indents the embedded lists, which looks weird
        json.dump(run_info, f, sort_keys=True, indent=2, separators=(',', ': '))


def plot_thresholds():
    """
    take our good channels, and plot the low energy region of
    their calibrated spectra.
    """
    # get list of wk/chan pairs
    pairs = []
    for wk in run_info["goodres_chans"]:
        for chan in run_info["goodres_chans"][wk]:
            pairs.append((wk, chan))

    cmap = plt.cm.get_cmap('nipy_spectral', len(pairs)+5)

    for i, (wk, chan) in enumerate(pairs):

        # load metadata
        runs = run_info["runs"][wk]
        gain = run_info["gain"][wk][chan]
        runtime = float(run_info["runtime"][wk])
        mass = float(run_info["det_mass_kg"])
        expo = (runtime/24) * mass # kg-d

        # load and clean data
        ch = load_data(wk, chan)
        ch.SetEstimate(ch.GetEntries())
        n = ch.Draw("Energy","Channel=={}".format(chan),"goff")
        ene = ch.GetV1()
        ene = np.asarray([ene[i] for i in range(n)])

        x_lo, x_hi, xpb = 0, 5000, 1
        xE, hE = get_hist(ene, x_lo, x_hi, xpb)

        xE = xE * gain
        epb = xE[1] - xE[0]
        hE = hE / expo / epb

        print("gain {:.2f}  rt {:.2f}  expo {:.2f}".format(gain, runtime, expo))

        plt.plot(xE, hE, c=cmap(i), ls='steps', lw=2,
                     label="wk {}, chan {}".format(wk, chan))

    plt.axvline(1460, c='r', label=r"$\mathregular{{}^{40}K}$")
    plt.xlim(10, 2000)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.ylabel("Counts / keV / kg-d", ha='right', y=1)
    plt.legend(loc=2, fontsize=12)
    plt.tight_layout()
    # plt.show()
    plt.savefig("./plots/thresholds.pdf")


def time_coins(wk, chan):
    """ identify alpha energies using sam hedges' delta-t method """
    print("tbd")


if __name__=="__main__":
    main()