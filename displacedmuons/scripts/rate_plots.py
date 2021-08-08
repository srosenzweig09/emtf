# from configparser import ConfigParser
# import awkward as ak
from icecream import ic
import numpy as np
np.seterr(all='ignore')
# import matplotlib as mpl
# from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
# from myuproot import open_up
import uproot
import sys

N_bx = 2736 # IP5, CMS (Run 3)
# https://indico.cern.ch/event/751857/contributions/3259413/attachments/1781638/3257666/ExperimentsInRun3_proceedings.pdf
f_LHC = 11.246; # khZ

prefix = 'root://eoscms.cern.ch/'

# filename = prefix + '/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/L1TMenuStudies/Run3full_rates_bdt/NuGun_11.05.21/nuGun_11May21.root'
filename = 'NuGun_11_3_0_pre5_NNv6_5M.root'

# tree, ak_table, np_table = open_up(filename, 'EMTFNtuple/tree')

t = uproot.open(filename + ":tree")

bdt_pt  = t['emtfTrack_pt'].array()
nn_pt   = t['emtfTrack_pt_dxy'].array()
nn_dxy  = t['emtfTrack_dxy'].array()
total = len(t['emtfTrack_pt'].array())

# ic(np.sum(nn_dxy >= 0))
# sys.exit()

# print out rate vs. pt_threshold
pt_thresholds = np.arange(300)

# full_rates_bdt = []
full_rates_nn = []
dxy1_rates_nn = []
dxy2_rates_nn = []
dxy3_rates_nn = []

for pt_cut in pt_thresholds:
    # bdt_rate = N_bx * f_LHC * np.sum(bdt_pt > pt_cut) / len(bdt_pt)
    pt_mask  = nn_pt  >  pt_cut
    dxy1_mask = abs(nn_dxy) >= 25
    dxy2_mask = abs(nn_dxy) >= 50
    dxy3_mask = abs(nn_dxy) >= 75

    nn_rate_all  = N_bx * f_LHC * np.sum(pt_mask) / len(nn_pt)
    nn_rate_dxy1 = N_bx * f_LHC * np.sum(pt_mask & dxy1_mask) / len(nn_pt)
    nn_rate_dxy2 = N_bx * f_LHC * np.sum(pt_mask & dxy2_mask) / len(nn_pt)
    nn_rate_dxy3 = N_bx * f_LHC * np.sum(pt_mask & dxy3_mask) / len(nn_pt)

    # full_rates_bdt.append(bdt_rate)
    full_rates_nn.append(nn_rate_all)
    dxy1_rates_nn.append(nn_rate_dxy1)
    dxy2_rates_nn.append(nn_rate_dxy2)
    dxy3_rates_nn.append(nn_rate_dxy3)

# full_rates_bdt = np.asarray(full_rates_bdt)
full_rates_nn  = np.asarray(full_rates_nn)
dxy1_rates_nn  = np.asarray(dxy1_rates_nn)
dxy2_rates_nn  = np.asarray(dxy2_rates_nn)
dxy3_rates_nn  = np.asarray(dxy3_rates_nn)

fig, ax = plt.subplots()
# ax.hist(pt_thresholds, bins=pt_thresholds, weights=full_rates_bdt, histtype='step', label='BDT')
ax.hist(pt_thresholds, bins=pt_thresholds, weights=full_rates_nn, histtype='step', label=r'L1 $d_{xy}$ > 0 cm')
ax.hist(pt_thresholds, bins=pt_thresholds, weights=dxy1_rates_nn, histtype='step', label=r'L1 $d_{xy}$ > 25 cm')
ax.hist(pt_thresholds, bins=pt_thresholds, weights=dxy2_rates_nn, histtype='step', label=r'L1 $d_{xy}$ > 50 cm')
ax.hist(pt_thresholds, bins=pt_thresholds, weights=dxy3_rates_nn, histtype='step', label=r'L1 $d_{xy}$ > 75 cm')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'L1 $p_T$ Threshold [GeV]')
ax.set_ylabel('Rate [kHz]')
ax.set_ylim(10**4)
ax.set_xlim(300)
ax.legend()

fig.savefig('rates/NuGun_11_3_0_pre5_NNv6_5M.pdf')
