from configparser import ConfigParser
import awkward1 as ak
import numpy as np
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
print(HtoLL)

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

print(gen_eta)

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

print("gen w/ cuts",np.sum((gen_eta_mask & gen_z0_d0_mask)*1))

BDT_mask = BDT_pt > 10 # GeV
NN_mask  = NN_pt > 10 # GeV

# fig, ax = plt.subplots()
bdt_mask = BDT_mask &  gen_eta_mask & gen_z0_d0_mask
print("bdt_mask",type(bdt_mask))
nn_mask = NN_mask &  gen_eta_mask & gen_z0_d0_mask
comb_mask = (BDT_mask | NN_mask) &  gen_eta_mask & gen_z0_d0_mask
denom_mask = gen_eta_mask & gen_z0_d0_mask
print("match emtf/gen w/ cuts",(np.sum((gen_eta_mask & gen_z0_d0_mask)*1)))
print("nn cuts",np.sum(nn_mask*1))
print("bdt cuts",np.sum(bdt_mask*1))
print("comb cuts",np.sum(comb_mask*1))

ptbins = np.linspace(0, 60, 30)
tester = gen_pt[comb_mask]
print(type(tester))
n_comb, n_edges = np.histogram(tester, bins=ptbins)
n_bdt, n_edges = np.histogram(gen_pt[bdt_mask], bins=ptbins)
n_nn, n_edges = np.histogram(gen_pt[nn_mask], bins=ptbins)
n_den, d_edges = np.histogram(gen_pt[denom_mask], bins=ptbins)

eff_comb = np.where(n_den == 0, 0, n_comb/n_den)
eff_bdt = np.where(n_den == 0, 0, n_bdt/n_den)
eff_nn = np.where(n_den == 0, 0, n_nn/n_den)

x = (n_edges[1:] + n_edges[:-1]) / 2
fig, ax = plt.subplots(figsize=(10,10))
ax.set_ylim(0,1.05)
ax.set_xlabel(r'gen muon $p_T$ [GeV]')
ax.set_ylabel('L1T Efficiency')

err_comb = np.sqrt(eff_comb*(1-eff_comb) / n_den)
err_bdt = np.sqrt(eff_bdt*(1-eff_bdt) / n_den)
err_nn = np.sqrt(eff_nn*(1-eff_nn) / n_den)

for i,w in enumerate(eff_comb):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='blue', label='Combined')
    ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='blue')
for i,w in enumerate(eff_bdt):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
    ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
for i,w in enumerate(eff_nn):
    ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='red',  label='NN')
    ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='red')

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='black', lw=2),
                Line2D([0], [0], color='red', lw=2)]

ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])

fig.savefig('efficiencies/efficiency_pt.pdf', bbox_inches=None)



fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$p_T^\mathrm{NN} - p_T^\mathrm{GEN}$ [GeV]')
ax.set_ylabel('Count')

hist(ax, ak.flatten(NN_pt[nn_mask])-ak.flatten(gen_pt[nn_mask]), bins=np.linspace(-60, 60, 60))

fig.savefig('efficiencies/nn_pt_difference.pdf', bbox_inches=None)


fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$\Delta p_T^\mathrm{NN}/p_T^\mathrm{GEN}$')
ax.set_ylabel('Count')

hist(ax, (ak.flatten(NN_pt[nn_mask])-ak.flatten(gen_pt[nn_mask]))/ak.flatten(gen_pt[nn_mask]), bins=np.linspace(-5, 15, 60))


fig.savefig('efficiencies/nn_pt_diff_norm.pdf', bbox_inches=None)



fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$p_T^\mathrm{BDT} - p_T^\mathrm{GEN}$ [GeV]')
ax.set_ylabel('Count')

hist(ax, ak.flatten(BDT_pt[bdt_mask])-ak.flatten(gen_pt[bdt_mask]), bins=np.linspace(-60, 60, 60))

fig.savefig('efficiencies/bdt_pt_difference.pdf', bbox_inches=None)


fig, ax = plt.subplots(figsize=(10,10))
ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT}/p_T^\mathrm{GEN}$')
ax.set_ylabel('Count')

hist(ax, (ak.flatten(BDT_pt[bdt_mask])-ak.flatten(gen_pt[bdt_mask]))/ak.flatten(gen_pt[bdt_mask]), bins=np.linspace(-5, 15, 60))


fig.savefig('efficiencies/bdt_pt_diff_norm.pdf', bbox_inches=None)


#############################################################
# 1D Distributions

# ptbins = np.linspace(0, 80, 40)

# fig, axs = plt.subplots(nrows=1, ncols=2)

# hist(axs[0], NN_pt, bins=ptbins)
# hist(axs[1], gen_pt, bins=ptbins)

# axs[0].set_xlabel("NN pT [GeV]")
# axs[1].set_xlabel("gen pT [GeV]")

# fig.savefig('distributions/pt')

# fig, ax = plt.subplots()
# bins = np.linspace(-3, 3, 40)
# # hist(ax, BDT_pt, bins=bins)
# hist(ax, NN_eta, bins=bins)
# fig.savefig('distributions/nn_eta')

# fig, ax = plt.subplots()
# bins = np.linspace(0, 80, 40)
# # hist(ax, BDT_pt, bins=bins)
# hist2d(ax, NN_pt, gen_pt, bins, bins)
# fig.savefig('distributions/pt_scatter')

# fig, ax = plt.subplots()
# bins = np.linspace(-3, 3, 40)
# # hist(ax, BDT_pt, bins=bins)
# hist2d(ax, NN_eta, gen_eta, bins, bins)
# fig.savefig('distributions/eta_scatter')
