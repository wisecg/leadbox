#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("~/dev/lat/pltReports.mplstyle")
from matplotlib.colors import LogNorm
from ROOT import TFile, TTree
from array import array

def main():
    """
    branch: Board                    79511
    branch: Channel                2877766
    branch: Timestamp             39708350
    branch: Energy                12015582
    branch: EnergyShort           11803163
    branch: Flags                   204907
    """
    run = 31
    chan = 1
    caen_file = "/Users/wisecg/project/leadbox/run_{}/FILTERED/compassF_run_{}.root".format(run,run)
    tf = TFile(caen_file)
    tt = tf.Get("Data_F")
    n = tt.Draw("Energy:EnergyShort","Channel==%d"%chan,"goff")
    ene, eshort = tt.GetV1(), tt.GetV2()
    ene = np.asarray([ene[i] for i in range(n)])
    eshort = np.asarray([eshort[i] for i in range(n)])

    rt = get_runtime(caen_file)
    izero = np.where(eshort == 0)
    nzero = len(izero[0])
    zero_rate = nzero/rt
    nent = len(ene)
    data_rate = nent/rt

    print("run {} runtime {:.2f}  nent {}  nzero {}  zero_rate {:.2f}  data_rate {:.2f}  ({:.2f}%)"
          .format(run, rt, nent, nzero, zero_rate, data_rate, 100*(zero_rate/data_rate)))

    icut = np.where((~np.isnan(eshort)) & (eshort > 10) & (ene > 10))

    plt.close()
    x_lo, x_hi, xpb = 0, 15000, 10
    nbx = int((x_hi-x_lo)/xpb)
    plt.hist(ene[icut], bins=nbx, range=(x_lo, x_hi), histtype='step', color='b',
             label="channel {}".format(chan))
    plt.xlabel("Energy [uncal]", ha='right', x=1)
    plt.ylabel("Counts", ha='right', y=1)
    plt.legend()
    plt.tight_layout()
    plt.savefig("./ene_1d_run{}.png".format(run), dpi=300)

    plt.close()
    x_lo, x_hi, xpb, y_lo, y_hi, ypb = 0, 6000, 10, 0, 6000, 10
    nbx, nby = int((y_hi-y_lo)/ypb), int((x_hi-x_lo)/xpb)
    plt.hist2d(ene, eshort, bins=[nby,nbx], range=[[x_lo,x_hi],[y_lo,y_hi]], norm=LogNorm())
    plt.colorbar()
    plt.xlabel("Energy", ha='right', x=1)
    plt.ylabel("EnergyShort", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("./eshort_2d_run{}.png".format(run), dpi=300)


def get_runtime(file_name, verbose=False):

    tf = TFile(file_name)
    tt = tf.Get("Data_F")
    n_ent = tt.GetEntries()
    ts = array('l',[0])
    tt.SetBranchAddress("Timestamp",ts)
    tt.GetEntry(0)
    ts_start = ts[0]/1e12 # caen's timestamp is evidently in picoseconds ... ugh
    tt.GetEntry(n_ent - 1)
    ts_stop = ts[0]/1e12
    runtime = (ts_stop - ts_start)/3600 # sec to hours
    if verbose:
        print("Start:{}  Stop:{}  Runtime (hrs):{}".format(ts_start, ts_stop, runtime))
    return runtime


if __name__=="__main__":
    main()