from configparser import ConfigParser
import awkward1 as ak
import numpy as np
np.seterr(all='ignore')
import matplotlib.pyplot as plt
import matplotlib as mpl

from myuproot import open_up
from kinematics import calcStar, calc_d0
from logger import info
from consistent_plots import hist, hist2d

#############################################################
# Import file info from config file
config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

HtoLL    = config['MuonGun']['filename']
treename = config['MuonGun']['treename']

tree, table, nptab = open_up(HtoLL, treename)

nevents = len(tree)

#############################################################
# Raw gen kinematics

gen_pt  = ak.flatten(table['gen_pt'] )
gen_phi = ak.flatten(table['gen_phi'])
gen_eta = ak.flatten(table['gen_eta'])
gen_vx  = ak.flatten(table['gen_vx'] )
gen_vy  = ak.flatten(table['gen_vy'] )
gen_vz  = ak.flatten(table['gen_vz'] )
gen_q   = ak.flatten(table['gen_q']  )

gen_eta, gen_phi = calcStar(gen_eta, gen_phi, gen_vx, gen_vy, gen_vz)

gen_d0 = calc_d0(gen_pt, gen_phi, gen_vx, gen_vy, gen_q, B=3.811)

print("matched gen/emtf",len(gen_pt))

#############################################################
# pT from EMTF

BDT_pt = ak.flatten(table['emtfTrack_pt'])
NN_pt  = ak.flatten(table['emtfTrack_pt_dxy'])

#############################################################
# Efficiency plots

gen_eta_mask = (abs(gen_eta) > 1.24) & (abs(gen_eta) < 2.5)
gen_z0_d0_mask =  (abs(gen_d0) < 100) & (abs(gen_vz) < 100) # cm, cm
gen_pt_mask = gen_pt > 15 # GeV
gen_d0_mask = abs(gen_d0) < 20

BDT_mask = BDT_pt > 10 # GeV
NN_mask  = NN_pt > 10 # GeV

def make_masks(pt_mask=True, eta_mask=False, d0_mask=False):
    bdt_mask = BDT_mask &  gen_z0_d0_mask
    nn_mask = NN_mask & gen_z0_d0_mask
    comb_mask = (BDT_mask | NN_mask) & gen_z0_d0_mask
    denom_mask = gen_eta_mask & gen_z0_d0_mask

    if eta_mask: 
        bdt_mask = gen_eta_mask & bdt_mask
        nn_mask = gen_eta_mask & nn_mask
        comb_mask = gen_eta_mask & comb_mask
        denom_mask = gen_eta_mask & denom_mask

    if pt_mask: 
        bdt_mask = gen_pt_mask & bdt_mask
        nn_mask = gen_pt_mask & nn_mask
        comb_mask = gen_pt_mask & comb_mask
        denom_mask = gen_pt_mask & denom_mask

    if d0_mask: 
        bdt_mask = gen_d0_mask & bdt_mask

    bdt_denom_mask = denom_mask & gen_d0_mask

    return bdt_mask, nn_mask, comb_mask, denom_mask, bdt_denom_mask

print("gen w/ cuts",np.sum((gen_eta_mask & gen_z0_d0_mask)*1))


pt_bdt_mask, pt_nn_mask, pt_comb_mask, pt_denom_mask, pt_bdt_denom_mask = make_masks(pt_mask=False, eta_mask=True, d0_mask=True)
eta_bdt_mask, eta_nn_mask, eta_comb_mask, eta_denom_mask, eta_bdt_denom_mask = make_masks(pt_mask=True, eta_mask=False, d0_mask=True)
d0_bdt_mask, d0_nn_mask, d0_comb_mask, d0_denom_mask, d0_bdt_denom_mask = make_masks(pt_mask=True, eta_mask=True, d0_mask=False)


###                     PT PLOTS


ptbins = np.linspace(0, 60, 30)
pt_n_comb, n_edges = np.histogram(ak.to_list(gen_pt[pt_comb_mask]), bins=ptbins)
pt_n_bdt, n_edges = np.histogram(ak.to_list(gen_pt[pt_bdt_mask]), bins=ptbins)
pt_n_nn, n_edges = np.histogram(ak.to_list(gen_pt[pt_nn_mask]), bins=ptbins)
pt_n_den, d_edges = np.histogram(ak.to_list(gen_pt[pt_denom_mask]), bins=ptbins)
pt_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_pt[pt_bdt_denom_mask]), bins=ptbins)

pt_eff_comb = np.where(pt_n_den == 0, 0, pt_n_comb/pt_n_den)
pt_eff_bdt = np.where(pt_n_den == 0, 0, pt_n_bdt/pt_n_den_bdt)
pt_eff_nn = np.where(pt_n_den == 0, 0, pt_n_nn/pt_n_den)

x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots(figsize=(10,10))
ax.set_ylim(0,1.05)
ax.set_xlabel(r'gen muon $p_T$ [GeV]')
ax.set_ylabel('L1T Efficiency')

pt_err_comb = np.sqrt(pt_eff_comb*(1-pt_eff_comb) / pt_n_den)
pt_err_bdt = np.sqrt(pt_eff_bdt*(1-pt_eff_bdt) / pt_n_den_bdt)
pt_err_nn = np.sqrt(pt_eff_nn*(1-pt_eff_nn) / pt_n_den)

for i,w in enumerate(pt_eff_comb):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='blue', label='Combined')
    ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='blue')
# for i,w in enumerate(eff_bdt):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
for i,w in enumerate(pt_eff_bdt):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='green', label='BDT')
    ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='green')
for i,w in enumerate(pt_eff_nn):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='red',  label='NN')
    ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='red')

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='black', lw=2),
                Line2D([0], [0], color='red', lw=2)]

ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
fig.savefig('efficiencies/pt_efficiency.pdf', bbox_inches=None)


fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$p_T^\mathrm{NN} - p_T^\mathrm{GEN}$ [GeV]')
ax.set_ylabel('Count')
hist(ax, ak.to_list(NN_pt[pt_nn_mask]-gen_pt[pt_nn_mask]), bins=np.linspace(-60, 60, 60))
fig.savefig('efficiencies/pt_nn_difference.pdf', bbox_inches=None)


fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$\Delta p_T^\mathrm{NN}/p_T^\mathrm{GEN}$')
ax.set_ylabel('Count')
hist(ax, (NN_pt[pt_nn_mask]-gen_pt[pt_nn_mask])/gen_pt[pt_nn_mask], bins=np.linspace(-5, 15, 60))
fig.savefig('efficiencies/pt_nn_diff_norm.pdf', bbox_inches=None)



fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$p_T^\mathrm{BDT} - p_T^\mathrm{GEN}$ [GeV]')
ax.set_ylabel('Count')
hist(ax, BDT_pt[pt_bdt_mask]-gen_pt[pt_bdt_mask], bins=np.linspace(-60, 60, 60))
fig.savefig('efficiencies/pt_bdt_difference.pdf', bbox_inches=None)


fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT}/p_T^\mathrm{GEN}$')
ax.set_ylabel('Count')
hist(ax, (BDT_pt[pt_bdt_mask]-gen_pt[pt_bdt_mask])/gen_pt[pt_bdt_mask], bins=np.linspace(-5, 15, 60))
fig.savefig('efficiencies/pt_bdt_diff_norm.pdf', bbox_inches=None)



###                     ETA PLOTS


etabins = [-2.5, -2.1, -1.6, -1.24]
for eta in etabins[::-1]:
    etabins.append(-eta)
eta_n_comb, n_edges = np.histogram(ak.to_list(gen_eta[eta_comb_mask]), bins=etabins)
eta_n_bdt, n_edges = np.histogram(ak.to_list(gen_eta[eta_bdt_mask]), bins=etabins)
eta_n_nn, n_edges = np.histogram(ak.to_list(gen_eta[eta_nn_mask]), bins=etabins)
eta_n_den, d_edges = np.histogram(ak.to_list(gen_eta[eta_denom_mask]), bins=etabins)
eta_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_eta[eta_bdt_denom_mask]), bins=etabins)

eta_eff_comb = np.where(eta_n_den == 0, 0, eta_n_comb/eta_n_den)
eta_eff_bdt = np.where(eta_n_den == 0, 0, eta_n_bdt/eta_n_den_bdt)
eta_eff_nn = np.where(eta_n_den == 0, 0, eta_n_nn/eta_n_den)

x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots(figsize=(10,10))
ax.set_ylim(0,1.05)
ax.set_xlabel(r'gen muon $p_T$ [GeV]')
ax.set_ylabel('L1T Efficiency')

eta_err_comb = np.sqrt(eta_eff_comb*(1-eta_eff_comb) / eta_n_den)
eta_err_bdt = np.sqrt(eta_eff_bdt*(1-eta_eff_bdt) / eta_n_den_bdt)
eta_err_nn = np.sqrt(eta_eff_nn*(1-eta_eff_nn) / eta_n_den)

for i,w in enumerate(eta_eff_comb):
    if i == 3: continue
    ax.plot([etabins[i], etabins[i+1]], [w,w], color='blue', label='Combined')
    ax.plot([x[i], x[i]], [w-eta_err_comb[i], w+eta_err_comb[i]], color='blue')
# for i,w in enumerate(eff_bdt):
#     ax.plot([etabins[i], etabins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
for i,w in enumerate(eta_eff_bdt):
    if i == 3: continue
    ax.plot([etabins[i], etabins[i+1]], [w,w], color='green', label='BDT')
    ax.plot([x[i], x[i]], [w-eta_err_comb[i], w+eta_err_comb[i]], color='green')
for i,w in enumerate(eta_eff_nn):
    if i == 3: continue
    ax.plot([etabins[i], etabins[i+1]], [w,w], color='red',  label='NN')
    ax.plot([x[i], x[i]], [w-eta_err_comb[i], w+eta_err_comb[i]], color='red')

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='black', lw=2),
                Line2D([0], [0], color='red', lw=2)]

ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
fig.savefig('efficiencies/eta_efficiency.pdf', bbox_inches=None)



###                     D0 PLOTS


d0bins = np.linspace(0,100,25)
d0_n_comb, n_edges = np.histogram(ak.to_list(gen_d0[d0_comb_mask]), bins=d0bins)
d0_n_bdt, n_edges = np.histogram(ak.to_list(gen_d0[d0_bdt_mask]), bins=d0bins)
d0_n_nn, n_edges = np.histogram(ak.to_list(gen_d0[d0_nn_mask]), bins=d0bins)
d0_n_den, d_edges = np.histogram(ak.to_list(gen_d0[d0_denom_mask]), bins=d0bins)
d0_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_d0[d0_bdt_denom_mask]), bins=d0bins)

d0_eff_comb = np.where(d0_n_den == 0, 0, d0_n_comb/d0_n_den)
d0_eff_bdt = np.where(d0_n_den == 0, 0, d0_n_bdt/d0_n_den_bdt)
d0_eff_nn = np.where(d0_n_den == 0, 0, d0_n_nn/d0_n_den)

x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots(figsize=(10,10))
ax.set_ylim(0,1.05)
ax.set_xlabel(r'gen muon $p_T$ [GeV]')
ax.set_ylabel('L1T Efficiency')

d0_err_comb = np.sqrt(d0_eff_comb*(1-d0_eff_comb) / d0_n_den)
d0_err_bdt = np.sqrt(d0_eff_bdt*(1-d0_eff_bdt) / d0_n_den_bdt)
d0_err_nn = np.sqrt(d0_eff_nn*(1-d0_eff_nn) / d0_n_den)

for i,w in enumerate(d0_eff_comb):
    ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='blue', label='Combined')
    ax.plot([x[i], x[i]], [w-d0_err_comb[i], w+d0_err_comb[i]], color='blue')
# for i,w in enumerate(eff_bdt):
#     ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
for i,w in enumerate(d0_eff_bdt):
    ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='green', label='BDT')
    ax.plot([x[i], x[i]], [w-d0_err_comb[i], w+d0_err_comb[i]], color='green')
for i,w in enumerate(d0_eff_nn):
    ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='red',  label='NN')
    ax.plot([x[i], x[i]], [w-d0_err_comb[i], w+d0_err_comb[i]], color='red')

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='black', lw=2),
                Line2D([0], [0], color='red', lw=2)]

ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
fig.savefig('efficiencies/d0_efficiency.pdf', bbox_inches=None)