from consistent_plots import hist, norm_hist
from icecream import ic
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
v6_nn_d0     = f6['nn_d0']
v6_gen_invpt = f6['gen_invpt']
v6_bdt_invpt = f6['bdt_invpt']
v6_nn_invpt  = f6['nn_invpt']


v6_delta_pt_bdt = v6_gen_pt - v6_bdt_pt
v6_delta_pt_nn  = v6_gen_pt - v6_nn_pt
v6_delta_d0_nn  = v6_gen_d0 - v6_nn_d0
v6_delta_invpt_bdt = v6_gen_invpt - v6_bdt_invpt
v6_delta_invpt_nn  = v6_gen_invpt - v6_nn_invpt

##########################################################################
## Load v10

f10 = np.load(f"dm_pt_v10.npz")

v10_gen_pt    = f10['gen_pt']
v10_gen_d0    = f10['gen_d0']
v10_bdt_pt    = f10['bdt_pt']
v10_nn_pt     = f10['nn_pt']
v10_nn_d0     = f10['nn_d0']
v10_gen_invpt = f10['gen_invpt']
v10_bdt_invpt = f10['bdt_invpt']
v10_nn_invpt  = f10['nn_invpt']

v10_delta_pt_bdt = v10_gen_pt - v10_bdt_pt
v10_delta_pt_nn  = v10_gen_pt - v10_nn_pt
v10_delta_d0_nn  = v10_gen_d0 - v10_nn_d0
v10_delta_invpt_bdt = v10_gen_invpt - v10_bdt_invpt
v10_delta_invpt_nn  = v10_gen_invpt - v10_nn_invpt

ic(v6_gen_d0, v6_nn_d0)
ic(v10_gen_d0, v10_nn_d0)

savedir = 'resolutions/'


##########################################################################
## Plotting

d0_bins = np.linspace(-100,100,100)
pt_bins = np.linspace(-25,25,100)
invpt_bins = np.linspace(-0.5,0.5,100)
frac_d0_bins = np.linspace(-1,1,100)
frac_pt_bins = np.linspace(-1,1,100)
frac_invpt_bins = np.linspace(-1,1,100)

bdt_pt_n_v6, b, x_pt = norm_hist(v6_delta_pt_bdt, pt_bins)
bdt_pt_n_v10, b, x = norm_hist(v10_delta_pt_bdt, pt_bins)
nn_pt_n_v6, b, x = norm_hist(v6_delta_pt_nn, pt_bins)
nn_pt_n_v10, b, x = norm_hist(v10_delta_pt_nn, pt_bins)
nn_d0_n_v6, b, x = norm_hist(v6_delta_d0_nn, d0_bins)
nn_d0_n_v10, b, x = norm_hist(v10_delta_d0_nn, d0_bins)
bdt_invpt_n_v6, b, x_invpt = norm_hist(v6_delta_invpt_bdt, invpt_bins)
bdt_invpt_n_v10, b, x = norm_hist(v10_delta_invpt_bdt, invpt_bins)
nn_invpt_n_v6, b, x = norm_hist(v6_delta_invpt_nn, invpt_bins)
nn_invpt_n_v10, b, x = norm_hist(v10_delta_invpt_nn, invpt_bins)

bdt_pt_n_v6_frac, b, x_pt = norm_hist(v6_delta_pt_bdt/v6_gen_pt, pt_bins)
bdt_pt_n_v10_frac, b, x = norm_hist(v10_delta_pt_bdt/v10_gen_pt, pt_bins)
nn_pt_n_v6_frac, b, x_pt = norm_hist(v6_delta_pt_nn/v6_gen_pt, pt_bins)
nn_pt_n_v10_frac, b, x = norm_hist(v10_delta_pt_nn/v10_gen_pt, pt_bins)
nn_d0_n_v6_frac, b, x_pt = norm_hist(v6_delta_d0_nn/v6_gen_d0, d0_bins)
nn_d0_n_v10_frac, b, x = norm_hist(v10_delta_d0_nn/v10_gen_d0, d0_bins)
bdt_invpt_n_v6_frac, b, x_pt = norm_hist(v6_delta_invpt_bdt/v6_gen_invpt, invpt_bins)
bdt_invpt_n_v10_frac, b, x = norm_hist(v10_delta_invpt_bdt/v10_gen_invpt, invpt_bins)
nn_invpt_n_v6_frac, b, x_pt = norm_hist(v6_delta_invpt_nn/v6_gen_invpt, invpt_bins)
nn_invpt_n_v10_frac, b, x = norm_hist(v10_delta_invpt_nn/v10_gen_invpt, invpt_bins)

fig, ax = plt.subplots()
hist(ax, x_pt, weights=bdt_pt_n_v6 , label='v6' , bins=pt_bins)
hist(ax, x_pt, weights=bdt_pt_n_v10, label='v10', bins=pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{gen-BDT}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta pT bdt resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_pt, weights=bdt_pt_n_v6_frac , label='v6' , bins=frac_pt_bins)
hist(ax, x_pt, weights=bdt_pt_n_v10_frac, label='v10', bins=frac_pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{gen-BDT}/p_T^\mathrm{gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta pT bdt frac resolution.pdf')


fig, ax = plt.subplots()
hist(ax, x_pt, weights=nn_pt_n_v6, label='v6'  , bins=pt_bins)
hist(ax, x_pt, weights=nn_pt_n_v10, label='v10', bins=pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{gen-NN}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta pT nn resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_pt, weights=nn_pt_n_v6_frac, label='v6'  , bins=frac_pt_bins)
hist(ax, x_pt, weights=nn_pt_n_v10_frac, label='v10', bins=frac_pt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{gen-NN}/p_T^\mathrm{gen}$')
ax.legend()
fig.savefig(savedir + 'normalized delta pT nn frac resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_pt, weights=nn_d0_n_v6, label='v6'  , bins=d0_bins)
hist(ax, x_pt, weights=nn_d0_n_v10, label='v10', bins=d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{gen-NN}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_pt, weights=nn_d0_n_v6_frac, label='v6'  , bins=frac_d0_bins)
hist(ax, x_pt, weights=nn_d0_n_v10_frac, label='v10', bins=frac_d0_bins)
ax.set_xlabel(r'$\Delta d_0^\mathrm{gen-NN}/d_0^\mathrm{gen}$ [cm]')
ax.legend()
fig.savefig(savedir + 'normalized delta d0 nn frac resolution.pdf')


fig, ax = plt.subplots()
hist(ax, x_invpt, weights=bdt_invpt_n_v6 , label='v6' , bins=invpt_bins)
hist(ax, x_invpt, weights=bdt_invpt_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{gen-BDT}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT bdt resolution.pdf')


fig, ax = plt.subplots()
hist(ax, x_invpt, weights=bdt_invpt_n_v6_frac , label='v6' , bins=frac_invpt_bins)
hist(ax, x_invpt, weights=bdt_invpt_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta p_T^\mathrm{gen-BDT}/p_T^\mathrm{gen}$ [GeV]')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT bdt frac resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_invpt, weights=nn_invpt_n_v6, label='v6'  , bins=invpt_bins)
hist(ax, x_invpt, weights=nn_invpt_n_v10, label='v10', bins=invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{gen-NN})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn resolution.pdf')

fig, ax = plt.subplots()
hist(ax, x_invpt, weights=nn_invpt_n_v6_frac, label='v6'  , bins=frac_invpt_bins)
hist(ax, x_invpt, weights=nn_invpt_n_v10_frac, label='v10', bins=frac_invpt_bins)
ax.set_xlabel(r'$\Delta (q^\mathrm{gen}/p_T^\mathrm{gen-NN})/(q^\mathrm{gen}/p_T^\mathrm{gen})$')
ax.legend()
fig.savefig(savedir + 'normalized delta invpT nn frac resolution.pdf')