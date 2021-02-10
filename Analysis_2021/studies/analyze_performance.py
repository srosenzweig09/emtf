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

config_file = 'input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

HtoLL    = config['HtoLL']['filename']
treename = config['HtoLL']['treename']

tree = get_uproot_Table(HtoLL, treename)

nevents = len(tree['genPart_pt'])
muon = (abs(tree['genPart_ID']) == 13) & (tree['genPart_parentID'] == 6000113)
nmuons = len(tree['genPart_pt'][muon].flatten())

matching = np.load('gen_emtf_matching_masks_0pt6.npz', allow_pickle=True)
muon_mask = matching['gen_muons']
gen_mask = matching['gen_mask']
emtf_mask = matching['emtf_index']

gen_pt  = tree['genPart_pt'][muon].flatten()[gen_mask]
gen_phi = tree['genPart_phi'][muon].flatten()[gen_mask]
gen_eta = tree['genPart_eta'][muon].flatten()[gen_mask]
gen_vx  = tree['genPart_vx'][muon].flatten()[gen_mask]
gen_vy  = tree['genPart_vy'][muon].flatten()[gen_mask]
gen_vz  = tree['genPart_vz'][muon].flatten()[gen_mask]
gen_q   = tree['genPart_q'][muon].flatten()[gen_mask]

gen_d0 = calc_d0(gen_pt, gen_phi, gen_vx, gen_vy, gen_q, B=3.811)

BDT_pt, L1_eta, L1_phi = convert_emtf(tree['emtfTrack_pt'],  tree['emtfTrack_eta'], tree['emtfTrack_phi'], flatten=True)

BDT_pt = BDT_pt[emtf_mask]
NN_pt = tree['emtfTrack_pt_dxy'].flatten()[emtf_mask]
# NN_pt = convert_emtf_pt(tree['emtfTrack_pt_dxy'], flatten=True)


BDT_mask = BDT_pt > 20 # GeV
NN_mask  = NN_pt > 20 # GeV


fig, ax = plt.subplots()

bins = np.linspace(0, 80, 40)
# hist(ax, BDT_pt, bins=bins)
hist(ax, NN_pt, bins=bins)
fig.savefig('test')

gen_eta_mask = (abs(gen_pt) > 1.2) & (abs(gen_eta) < 2.4)
gen_z0_d0_mask =  (abs(gen_d0) < 100) & (abs(gen_vz) < 100) # cm, cm
gen_pt_mask = gen_pt > 15 # GeV

# fig, ax = plt.subplots()
num_mask = BDT_mask | NN_mask
denom_mask = gen_eta_mask & gen_z0_d0_mask & gen_pt_mask

bins = np.linspace(0, 60, 25)
n_num, n_edges = np.histogram(NN_pt[num_mask], bins=bins)
n_den, d_edges = np.histogram(gen_pt[denom_mask], bins=bins)

print(n_num)
print(n_den)

eff = np.where(n_den == 0, 0, n_num/n_den)
print(eff)
x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots()
err = np.sqrt(eff*(1-eff) / n_den)
    
for i,w in enumerate(eff):
    ax.plot([bins[i], bins[i+1]], [w,w], color='C0')
    ax.plot([x[i], x[i]], [w-err[i], w+err[i]], color='C0')
fig.savefig('test2')


# n_num, n_edges = np.histogram()