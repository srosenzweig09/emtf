from configparser import ConfigParser
import numpy as np
import awkward0 as awk
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

from uproot_open import get_uproot_Table
from kinematics import calcDeltaR, calcStar, convert_emtf, convert_emtf_pt, calc_d0
from logger import info
from consistent_plots import hist, hist2d

#############################################################
# Import file info from config file
config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

HtoLL    = config['HtoLL']['filename']
treename = config['HtoLL']['treename']

tree = get_uproot_Table(HtoLL, treename)

nevents = len(tree['genPart_pt'])
muon = (abs(tree['genPart_ID']) == 13) & (tree['genPart_parentID'] == 6000113)
nmuons = len(tree['genPart_pt'][muon].flatten())

gen_emtf_matching  = np.load('matching/gen_gmt/gen_gmt_matching_masks_0pt6.npz', allow_pickle=True)
gen_index = gen_emtf_matching['gen_indices']
emtf_index = gen_emtf_matching['emtf_indices']
gmt_index = gen_emtf_matching['gmt_indices']

#############################################################
## Matched gen kinematics

matched_gen_pt  = tree['genPart_pt'] [muon].flatten()[gen_index]
matched_gen_phi = tree['genPart_phi'][muon].flatten()[gen_index]
matched_gen_eta = tree['genPart_eta'][muon].flatten()[gen_index]
matched_gen_vx  = tree['genPart_vx'] [muon].flatten()[gen_index]
matched_gen_vy  = tree['genPart_vy'] [muon].flatten()[gen_index]
matched_gen_vz  = tree['genPart_vz'] [muon].flatten()[gen_index]
matched_gen_q   = tree['genPart_q']  [muon].flatten()[gen_index]

matched_gen_d0 = calc_d0(matched_gen_pt, matched_gen_phi, matched_gen_vx, matched_gen_vy, matched_gen_q, B=3.811)

#############################################################
# Raw gen kinematics

gen_pt  = tree['genPart_pt'] [muon].flatten()
gen_phi = tree['genPart_phi'][muon].flatten()
gen_eta = tree['genPart_eta'][muon].flatten()
gen_vx  = tree['genPart_vx'] [muon].flatten()
gen_vy  = tree['genPart_vy'] [muon].flatten()
gen_vz  = tree['genPart_vz'] [muon].flatten()
gen_q   = tree['genPart_q']  [muon].flatten()

gen_d0 = calc_d0(gen_pt, gen_phi, gen_vx, gen_vy, gen_q, B=3.811)

#############################################################
# pT from EMTF

BDT_pt, L1_eta, L1_phi = convert_emtf(tree['emtfTrack_pt'],  tree['emtfTrack_eta'], tree['emtfTrack_phi'], flatten=True)

BDT_pt = BDT_pt[emtf_index]
NN_pt = tree['emtfTrack_pt_dxy'].flatten()[emtf_index]
NN_eta = tree['gmtMuon_eta'].flatten()[gmt_index]
# NN_pt = convert_emtf_pt(tree['emtfTrack_pt_dxy'], flatten=True)

BDT_mask = BDT_pt > 20 # GeV
NN_mask  = NN_pt > 20 # GeV

#############################################################
# 1D Distributions

ptbins = np.linspace(0, 80, 40)

fig, axs = plt.subplots(nrows=1, ncols=2)

hist(axs[0], NN_pt, bins=ptbins)
hist(axs[1], matched_gen_pt, bins=ptbins)

axs[0].set_xlabel("NN pT [GeV]")
axs[1].set_xlabel("gen pT [GeV]")

fig.savefig('distributions/matched_pt')

fig, ax = plt.subplots()
bins = np.linspace(-3, 3, 40)
# hist(ax, BDT_pt, bins=bins)
hist(ax, NN_eta, bins=bins)
fig.savefig('distributions/matched_nn_eta')

fig, ax = plt.subplots()
bins = np.linspace(0, 80, 40)
# hist(ax, BDT_pt, bins=bins)
hist2d(ax, NN_pt, matched_gen_pt, bins, bins)
fig.savefig('distributions/matched_pt_scatter')

fig, ax = plt.subplots()
bins = np.linspace(-3, 3, 40)
# hist(ax, BDT_pt, bins=bins)
hist2d(ax, NN_eta, matched_gen_eta, bins, bins)
fig.savefig('distributions/matched_eta_scatter')

#############################################################
# Efficiency plots

matched_gen_eta_mask = (abs(matched_gen_eta) > 1.2) & (abs(matched_gen_eta) < 2.4)
matched_gen_z0_d0_mask =  (abs(matched_gen_d0) < 100) & (abs(matched_gen_vz) < 100) # cm, cm
matched_gen_pt_mask = matched_gen_pt > 15 # GeV

gen_eta_mask = (abs(gen_eta) > 1.24) & (abs(gen_eta) < 2.5)
gen_z0_d0_mask =  (abs(gen_d0) < 100) & (abs(gen_vz) < 100) # cm, cm
gen_pt_mask = gen_pt > 15 # GeV

# fig, ax = plt.subplots()
num_mask = (BDT_mask | NN_mask) &  matched_gen_eta_mask & matched_gen_z0_d0_mask
denom_mask = gen_eta_mask & gen_z0_d0_mask

ptbins = np.linspace(0, 60, 25)
n_num, n_edges = np.histogram(matched_gen_pt[num_mask].flatten(), bins=ptbins)
n_den, d_edges = np.histogram(gen_pt[denom_mask].flatten(), bins=ptbins)

fig, ax = plt.subplots()

hist(ax, matched_gen_pt[num_mask].flatten(), bins=ptbins)
hist(ax, gen_pt[denom_mask].flatten(), bins=ptbins)

fig.savefig('distributions/masked_gen_pt')

print(n_num)
print(n_den)

eff = np.where(n_den == 0, 0, n_num/n_den)
print(eff)
x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots()
err = np.sqrt(eff*(1-eff) / n_den)
ax.set_ylim(0,1.2)
    
for i,w in enumerate(eff):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='C0')
    ax.plot([x[i], x[i]], [w-err[i], w+err[i]], color='C0')
fig.savefig('efficiencies/efficiency_pt', bbox_inches=None)


# # n_num, n_edges = np.histogram()