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

matched_gen_eta, matched_gen_phi = calcStar(matched_gen_eta, matched_gen_phi, matched_gen_vx, matched_gen_vy, matched_gen_vz, is2darray=False)

matched_gen_d0 = calc_d0(matched_gen_pt, matched_gen_phi, matched_gen_vx, matched_gen_vy, matched_gen_q, B=3.811)

#############################################################
# Raw gen kinematics

gen_pt  = tree['genPart_pt'] [muon].flatten()
# print(len(gen_pt))
gen_phi = tree['genPart_phi'][muon].flatten()
gen_eta = tree['genPart_eta'][muon].flatten()
gen_vx  = tree['genPart_vx'] [muon].flatten()
gen_vy  = tree['genPart_vy'] [muon].flatten()
gen_vz  = tree['genPart_vz'] [muon].flatten()
gen_q   = tree['genPart_q']  [muon].flatten()

gen_eta, gen_phi = calcStar(gen_eta, gen_phi, gen_vx, gen_vy, gen_vz, is2darray=False)

gen_d0 = calc_d0(gen_pt, gen_phi, gen_vx, gen_vy, gen_q, B=3.811)

#############################################################
# pT from EMTF

BDT_pt, L1_eta, L1_phi = convert_emtf(tree['emtfTrack_pt'],  tree['emtfTrack_eta'], tree['emtfTrack_phi'], flatten=True)

BDT_pt = BDT_pt[emtf_index]
NN_pt = tree['emtfTrack_pt_dxy'].flatten()[emtf_index]
NN_eta = tree['gmtMuon_eta'].flatten()[gmt_index]

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
bdt_mask = BDT_mask &  matched_gen_eta_mask & matched_gen_z0_d0_mask
nn_mask = (NN_mask) &  matched_gen_eta_mask & matched_gen_z0_d0_mask
comb_mask = (BDT_mask | NN_mask) &  matched_gen_eta_mask & matched_gen_z0_d0_mask
denom_mask = gen_eta_mask & gen_z0_d0_mask
print(np.sum((comb_mask)*1))

ptbins = np.linspace(0, 60, 30)
n_comb, n_edges = np.histogram(matched_gen_pt[comb_mask].flatten(), bins=ptbins)
n_bdt, n_edges = np.histogram(matched_gen_pt[bdt_mask].flatten(), bins=ptbins)
n_nn, n_edges = np.histogram(matched_gen_pt[nn_mask].flatten(), bins=ptbins)
n_den, d_edges = np.histogram(gen_pt[denom_mask].flatten(), bins=ptbins)

eff_comb = np.where(n_den == 0, 0, n_comb/n_den)
eff_bdt = np.where(n_den == 0, 0, n_bdt/n_den)
eff_nn = np.where(n_den == 0, 0, n_nn/n_den)

x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots()
ax.set_ylim(0,1.2)
ax.set_xlabel(r'gen muon $p_T$ [GeV]')
ax.set_ylabel('L1T Efficiency')

err_comb = np.sqrt(eff_comb*(1-eff_comb) / n_den)
err_bdt = np.sqrt(eff_bdt*(1-eff_bdt) / n_den)
err_nn = np.sqrt(eff_nn*(1-eff_nn) / n_den)

for i,w in enumerate(eff_comb):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='C0', label='Combined')
    ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='C0')
for i,w in enumerate(eff_bdt):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='C1', label='BDT')
    ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='C1')
for i,w in enumerate(eff_nn):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='C2',  label='NN')
    ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='C2')

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='C0', lw=2),
                Line2D([0], [0], color='C1', lw=2),
                Line2D([0], [0], color='C2', lw=2)]

ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])

fig.savefig('efficiencies/efficiency_pt', bbox_inches=None)

