from consistent_plots import hist, norm_hist
from icecream import ic
import matplotlib as mpl
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
import matplotlib.pyplot as plt
import numpy as np
np.seterr(all='ignore')
import sys


##########################################################################
## Load v6 

f6 = np.load(f"dm_pt_v6.npz")

v6_gen_pt    = f6['gen_pt']
v6_gen_d0    = f6['gen_d0']
v6_bdt_pt    = f6['bdt_pt']
v6_nn_pt     = f6['nn_pt']
v6_nn_d0     = -f6['nn_d0']
v6_gen_invpt = f6['gen_invpt']
v6_bdt_invpt = f6['bdt_invpt']
v6_nn_invpt  = f6['nn_invpt']
v6_gen_eta  = f6['gen_eta']
v6_lo_eta = np.logical_and(abs(v6_gen_eta) > 1.2, abs(v6_gen_eta) < 1.6)
v6_md_eta = np.logical_and(abs(v6_gen_eta) > 1.6, abs(v6_gen_eta) < 2.1)
v6_hi_eta = np.logical_and(abs(v6_gen_eta) > 2.1, abs(v6_gen_eta) < 2.5)

v6_delta_pt_bdt = -v6_gen_pt + v6_bdt_pt
v6_delta_pt_nn  = -v6_gen_pt + v6_nn_pt
v6_delta_d0_nn  = -v6_gen_d0 + v6_nn_d0
v6_delta_invpt_bdt = -v6_gen_invpt + v6_bdt_invpt
v6_delta_invpt_nn  = -v6_gen_invpt + v6_nn_invpt

##########################################################################
## Load v10

f10 = np.load(f"dm_pt_v10.npz")

v10_gen_pt    = f10['gen_pt']
v10_gen_d0    = f10['gen_d0']
v10_bdt_pt    = f10['bdt_pt']
v10_nn_pt     = f10['nn_pt']
v10_nn_d0     = -f10['nn_d0']
v10_gen_invpt = f10['gen_invpt']
v10_bdt_invpt = f10['bdt_invpt']
v10_nn_invpt  = f10['nn_invpt']
v10_gen_eta  = f10['gen_eta']
v10_lo_eta = np.logical_and(abs(v10_gen_eta) > 1.2, abs(v10_gen_eta) < 1.6)
v10_md_eta = np.logical_and(abs(v10_gen_eta) > 1.6, abs(v10_gen_eta) < 2.1)
v10_hi_eta = np.logical_and(abs(v10_gen_eta) > 2.1, abs(v10_gen_eta) < 2.5)

v10_delta_pt_bdt = -v10_gen_pt + v10_bdt_pt
v10_delta_pt_nn  = -v10_gen_pt + v10_nn_pt
v10_delta_d0_nn  = -v10_gen_d0 + v10_nn_d0
v10_delta_invpt_bdt = -v10_gen_invpt + v10_bdt_invpt
v10_delta_invpt_nn  = -v10_gen_invpt + v10_nn_invpt

savedir = 'resolutions/'


##########################################################################
## Plotting

d0_bins = np.linspace(-250,250,100)
pt_bins = np.linspace(-25,25,100)
invpt_bins = np.linspace(-0.5,0.5,100)
frac_d0_bins = np.linspace(-10,10,30)
frac_pt_bins = np.linspace(-5,10,30)
frac_invpt_bins = np.linspace(-1,1,30)

bdt_pt_n_v6, b, x_pt = norm_hist(v6_delta_pt_bdt, pt_bins)
bdt_pt_n_v10, b, x = norm_hist(v10_delta_pt_bdt, pt_bins)
nn_pt_n_v6, b, x = norm_hist(v6_delta_pt_nn, pt_bins)
nn_pt_n_v10, b, x = norm_hist(v10_delta_pt_nn, pt_bins)
nn_d0_n_v6, b, x_d0 = norm_hist(v6_delta_d0_nn, d0_bins)
nn_d0_n_v10, b, x = norm_hist(v10_delta_d0_nn, d0_bins)
nn_d0_lo_eta_n_v6, b, x = norm_hist(v6_delta_d0_nn[v6_lo_eta], d0_bins)
nn_d0_lo_eta_n_v10, b, x = norm_hist(v10_delta_d0_nn[v10_lo_eta], d0_bins)
nn_d0_md_eta_n_v6, b, x = norm_hist(v6_delta_d0_nn[v6_md_eta], d0_bins)
nn_d0_md_eta_n_v10, b, x = norm_hist(v10_delta_d0_nn[v10_md_eta], d0_bins)
nn_d0_hi_eta_n_v6, b, x = norm_hist(v6_delta_d0_nn[v6_hi_eta], d0_bins)
nn_d0_hi_eta_n_v10, b, x = norm_hist(v10_delta_d0_nn[v10_hi_eta], d0_bins)
bdt_invpt_n_v6, b, x_invpt = norm_hist(v6_delta_invpt_bdt, invpt_bins)
bdt_invpt_n_v10, b, x = norm_hist(v10_delta_invpt_bdt, invpt_bins)
nn_invpt_n_v6, b, x = norm_hist(v6_delta_invpt_nn, invpt_bins)
nn_invpt_n_v10, b, x = norm_hist(v10_delta_invpt_nn, invpt_bins)
nn_invpt_lo_eta_n_v6, b, x = norm_hist(v6_delta_invpt_nn[v6_lo_eta], invpt_bins)
nn_invpt_lo_eta_n_v10, b, x = norm_hist(v10_delta_invpt_nn[v10_lo_eta], invpt_bins)
nn_invpt_md_eta_n_v6, b, x = norm_hist(v6_delta_invpt_nn[v6_md_eta], invpt_bins)
nn_invpt_md_eta_n_v10, b, x = norm_hist(v10_delta_invpt_nn[v10_md_eta], invpt_bins)
nn_invpt_hi_eta_n_v6, b, x = norm_hist(v6_delta_invpt_nn[v6_hi_eta], invpt_bins)
nn_invpt_hi_eta_n_v10, b, x = norm_hist(v10_delta_invpt_nn[v10_hi_eta], invpt_bins)

bdt_pt_n_v6_frac, b, x_frac_pt = norm_hist(v6_delta_pt_bdt/v6_gen_pt, frac_pt_bins)
bdt_pt_n_v10_frac, b, x = norm_hist(v10_delta_pt_bdt/v10_gen_pt, frac_pt_bins)
nn_pt_n_v6_frac, b, x = norm_hist(v6_delta_pt_nn/v6_gen_pt, frac_pt_bins)
nn_pt_n_v10_frac, b, x = norm_hist(v10_delta_pt_nn/v10_gen_pt, frac_pt_bins)
nn_d0_n_v6_frac, b, x_frac_d0 = norm_hist(v6_delta_d0_nn/v6_gen_d0, frac_d0_bins)
nn_d0_n_v10_frac, b, x = norm_hist(v10_delta_d0_nn/v10_gen_d0, frac_d0_bins)
nn_d0_lo_eta_n_v6_frac, b, x = norm_hist((v6_delta_d0_nn/v6_gen_eta)[v6_lo_eta], frac_d0_bins)
nn_d0_lo_eta_n_v10_frac, b, x = norm_hist((v10_delta_d0_nn/v10_gen_eta)[v10_lo_eta], frac_d0_bins)
nn_d0_md_eta_n_v6_frac, b, x = norm_hist((v6_delta_d0_nn/v6_gen_eta)[v6_md_eta], frac_d0_bins)
nn_d0_md_eta_n_v10_frac, b, x = norm_hist((v10_delta_d0_nn/v10_gen_eta)[v10_md_eta], frac_d0_bins)
nn_d0_hi_eta_n_v6_frac, b, x = norm_hist((v6_delta_d0_nn/v6_gen_eta)[v6_hi_eta], frac_d0_bins)
nn_d0_hi_eta_n_v10_frac, b, x = norm_hist((v10_delta_d0_nn/v10_gen_eta)[v10_hi_eta], frac_d0_bins)
bdt_invpt_n_v6_frac, b, x_frac_invpt = norm_hist(v6_delta_invpt_bdt/v6_gen_invpt, frac_invpt_bins)
bdt_invpt_n_v10_frac, b, x = norm_hist(v10_delta_invpt_bdt/v10_gen_invpt, frac_invpt_bins)
nn_invpt_n_v6_frac, b, x = norm_hist(v6_delta_invpt_nn/v6_gen_invpt, frac_invpt_bins)
nn_invpt_n_v10_frac, b, x = norm_hist(v10_delta_invpt_nn/v10_gen_invpt, frac_invpt_bins)
nn_invpt_lo_eta_n_v6_frac, b, x = norm_hist((v6_delta_invpt_nn/v6_gen_invpt)[v6_lo_eta], frac_invpt_bins)
nn_invpt_lo_eta_n_v10_frac, b, x = norm_hist((v10_delta_invpt_nn/v10_gen_invpt)[v10_lo_eta], frac_invpt_bins)
nn_invpt_md_eta_n_v6_frac, b, x = norm_hist((v6_delta_invpt_nn/v6_gen_invpt)[v6_md_eta], frac_invpt_bins)
nn_invpt_md_eta_n_v10_frac, b, x = norm_hist((v10_delta_invpt_nn/v10_gen_invpt)[v10_md_eta], frac_invpt_bins)
nn_invpt_hi_eta_n_v6_frac, b, x = norm_hist((v6_delta_invpt_nn/v6_gen_invpt)[v6_hi_eta], frac_invpt_bins)
nn_invpt_hi_eta_n_v10_frac, b, x = norm_hist((v10_delta_invpt_nn/v10_gen_invpt)[v10_hi_eta], frac_invpt_bins)

fig, ax = plt.subplots()
hist(ax, x_pt, weights=bdt_pt_n_v6 , label='v6' , bins=pt_bins)
hist(ax, x_pt, weights=bdt_pt_n_v10, label='v10', bins=pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT-gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta pT bdt resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_frac_pt, weights=bdt_pt_n_v6_frac , label='v6' , bins=frac_pt_bins)
hist(ax, x_frac_pt, weights=bdt_pt_n_v10_frac, label='v10', bins=frac_pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT-gen}/p_T^\mathrm{gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta pT bdt frac resolution.pdf')


fig, ax = plt.subplots()
hist(ax, x_pt, weights=nn_pt_n_v6, label='v6'  , bins=pt_bins)
hist(ax, x_pt, weights=nn_pt_n_v10, label='v10', bins=pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{NN-gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta pT nn resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_frac_pt, weights=nn_pt_n_v6_frac, label='v6'  , bins=frac_pt_bins)
hist(ax, x_frac_pt, weights=nn_pt_n_v10_frac, label='v10', bins=frac_pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{NN-gen}/p_T^\mathrm{gen}$')
ax.legend()
fig.savefig(savedir + 'normalized delta pT nn frac resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_d0, weights=nn_d0_n_v6, label='v6'  , bins=d0_bins)
hist(ax, x_d0, weights=nn_d0_n_v10, label='v10', bins=d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn resolution.pdf')
fig, ax = plt.subplots()
hist(ax, x_d0, weights=nn_d0_lo_eta_n_v6, label='v6'  , bins=d0_bins)
hist(ax, x_d0, weights=nn_d0_lo_eta_n_v10, label='v10', bins=d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn resolution lo eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_d0, weights=nn_d0_md_eta_n_v6, label='v6'  , bins=d0_bins)
hist(ax, x_d0, weights=nn_d0_md_eta_n_v10, label='v10', bins=d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn resolution md eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_d0, weights=nn_d0_hi_eta_n_v6, label='v6'  , bins=d0_bins)
hist(ax, x_d0, weights=nn_d0_hi_eta_n_v10, label='v10', bins=d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn resolution hi eta.pdf')

fig, ax = plt.subplots()
hist(ax, x_frac_d0, weights=nn_d0_n_v6_frac, label='v6'  , bins=frac_d0_bins)
hist(ax, x_frac_d0, weights=nn_d0_n_v10_frac, label='v10', bins=frac_d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}/d_0^\mathrm{gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn frac resolution.pdf')
fig, ax = plt.subplots()
hist(ax, x_frac_d0, weights=nn_d0_lo_eta_n_v6_frac, label='v6'  , bins=frac_d0_bins)
hist(ax, x_frac_d0, weights=nn_d0_lo_eta_n_v10_frac, label='v10', bins=frac_d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn frac resolution lo eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_frac_d0, weights=nn_d0_md_eta_n_v6_frac, label='v6'  , bins=frac_d0_bins)
hist(ax, x_frac_d0, weights=nn_d0_md_eta_n_v10_frac, label='v10', bins=frac_d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn frac resolution md eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_frac_d0, weights=nn_d0_hi_eta_n_v6_frac, label='v6'  , bins=frac_d0_bins)
hist(ax, x_frac_d0, weights=nn_d0_hi_eta_n_v10_frac, label='v10', bins=frac_d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{NN-gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn frac resolution hi eta.pdf')


fig, ax = plt.subplots()
hist(ax, x_invpt, weights=bdt_invpt_n_v6 , label='v6' , bins=invpt_bins)
hist(ax, x_invpt, weights=bdt_invpt_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT-gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT bdt resolution.pdf')


fig, ax = plt.subplots()
hist(ax, x_frac_invpt, weights=bdt_invpt_n_v6_frac , label='v6' , bins=frac_invpt_bins)
hist(ax, x_frac_invpt, weights=bdt_invpt_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT-gen}/p_T^\mathrm{gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT bdt frac resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_invpt, weights=nn_invpt_n_v6, label='v6'  , bins=invpt_bins)
hist(ax, x_invpt, weights=nn_invpt_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn resolution.pdf')
fig, ax = plt.subplots()
hist(ax, x_invpt, weights=nn_invpt_lo_eta_n_v6, label='v6'  , bins=invpt_bins)
hist(ax, x_invpt, weights=nn_invpt_lo_eta_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn resolution lo eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_invpt, weights=nn_invpt_md_eta_n_v6, label='v6'  , bins=invpt_bins)
hist(ax, x_invpt, weights=nn_invpt_md_eta_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn resolution md eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_invpt, weights=nn_invpt_hi_eta_n_v6, label='v6'  , bins=invpt_bins)
hist(ax, x_invpt, weights=nn_invpt_hi_eta_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn resolution hi eta.pdf')

fig, ax = plt.subplots()
hist(ax, x_frac_invpt, weights=nn_invpt_n_v6_frac, label='v6'  , bins=frac_invpt_bins)
hist(ax, x_frac_invpt, weights=nn_invpt_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})/(q^\mathrm{gen}/p_T^\mathrm{gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn frac resolution.pdf')
fig, ax = plt.subplots()
hist(ax, x_frac_invpt, weights=nn_invpt_lo_eta_n_v6_frac, label='v6'  , bins=frac_invpt_bins)
hist(ax, x_frac_invpt, weights=nn_invpt_lo_eta_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn frac resolution lo eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_frac_invpt, weights=nn_invpt_md_eta_n_v6_frac, label='v6'  , bins=frac_invpt_bins)
hist(ax, x_frac_invpt, weights=nn_invpt_md_eta_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn frac resolution md eta.pdf')
fig, ax = plt.subplots()
hist(ax, x_frac_invpt, weights=nn_invpt_hi_eta_n_v6_frac, label='v6'  , bins=frac_invpt_bins)
hist(ax, x_frac_invpt, weights=nn_invpt_hi_eta_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{NN-gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn frac resolution hi eta.pdf')