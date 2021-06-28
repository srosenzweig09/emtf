from configparser import ConfigParser
import awkward as ak
from icecream import ic
import numpy as np
np.seterr(all='ignore')
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import uproot

# from myuproot import open_up
from kinematics import calcStar, calc_d0
from logger import info

def get_n(arr, bins, mask):
    n, edges = np.histogram(ak.to_list(arr[mask]), bins=bins)
    return n

def plot_eff(n_bdt, n_nn, n_both, n_den, bins, title, xlabel):
    global folder
    eff_both = np.where(n_den == 0, 0, n_both/n_den)
    eff_bdt = np.where(n_den == 0, 0, n_bdt/n_den)
    eff_nn = np.where(n_den == 0, 0, n_nn/n_den)

    err_both = np.sqrt(eff_both*(1-eff_both) / n_den)
    err_bdt = np.sqrt(eff_bdt*(1-eff_bdt) / n_den)
    err_nn = np.sqrt(eff_nn*(1-eff_nn) / n_den)

    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_ylim(0,1.05)

    x = (bins[1:] + bins[:-1]) / 2

    for i,w in enumerate(eff_both):
        ax.plot([bins[i], bins[i+1]], [w,w], color='blue') # BDT or NN
        ax.plot([x[i], x[i]], [w-err_both[i], w+err_both[i]], color='blue')
    for i,w in enumerate(eff_bdt):
        ax.plot([bins[i], bins[i+1]], [w,w], color='black') # BDT
        ax.plot([x[i], x[i]], [w-err_bdt[i], w+err_bdt[i]], color='black')
    for i,w in enumerate(eff_nn):
        ax.plot([bins[i], bins[i+1]], [w,w], color='red') # NN
        ax.plot([x[i], x[i]], [w-err_nn[i], w+err_nn[i]], color='red')

    custom_lines = [Line2D([0], [0], color='blue', lw=2),
                    Line2D([0], [0], color='black', lw=2),
                    Line2D([0], [0], color='red', lw=2)]

    ax.legend(custom_lines, [f'{sum(n_both)} : BBDT or NN', f'{sum(n_bdt)} : BDT', f'{sum(n_nn)} : NN'])
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Efficiency')

    plt.grid(color='lightgrey', linestyle='--', linewidth=1, which='both')

    props = dict(boxstyle='round')
    ax.text(0.85, 0.8, f'Denom: {sum(n_den)}')

    fig.savefig(f'{folder}{title}.pdf', bbox_inches=None)

#############################################################
# Import file info from config file
config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

filename1 = config['MuonGun']['filename'] + "_1.root"
filename2 = config['MuonGun']['filename'] + "_2.root"
filename3 = config['MuonGun']['filename'] + "_3.root"
filename4 = config['MuonGun']['filename'] + "_4.root"
treename = config['MuonGun']['treename']
version = config['MuonGun']['version']

folder = f'efficiencies/v{version}/'

# filelist = [filename1, filename2, filename3, filename4]
filelist = [filename1]

n_pt_bdt  = 0
n_pt_nn   = 0
n_pt_both = 0

ptbins = np.linspace(0, 60, 30)
etabins = [-2.5, -2.1, -1.6, -1.24]
for eta in etabins[::-1]:
    etabins.append(-eta)
etabins = np.array((etabins))
ic(etabins)

for filename in filelist:
    info("Loading file")
    # info(filename)

    ic(filename)

    with uproot.open(filename) as file:
        info("Loaded file.")

        tree = file['tree']
        table = tree.arrays()

        info("Extracting branches.")

        gen_pt      = ak.flatten(table['gen_pt'])
        gen_phi     = ak.flatten(table['gen_phi'])
        gen_eta     = ak.flatten(table['gen_etaStar'])
        gen_phiStar = ak.flatten(table['gen_phiStar'])
        # gen_etaStar = ak.flatten(table['gen_etaStar'])
        gen_q       = ak.flatten(table['gen_q'])
        gen_d0      = ak.flatten(table['gen_d0'])
        gen_vz      = ak.flatten(table['gen_vz'])
        gen_dR      = ak.flatten(table['gen_dR'])

        matched_gen_pt      = ak.flatten(table['matched_gen_pt'])
        matched_gen_phi     = ak.flatten(table['matched_gen_phi'])
        matched_gen_eta     = ak.flatten(table['matched_gen_etaStar'])
        matched_gen_phiStar = ak.flatten(table['matched_gen_phiStar'])
        # matched_gen_etaStar = ak.flatten(table['matched_gen_etaStar'])
        matched_gen_q       = ak.flatten(table['matched_gen_q'])
        matched_gen_d0      = ak.flatten(table['matched_gen_d0'])
        matched_gen_dR      = ak.flatten(table['matched_gen_dR'])
        matched_gen_vz      = ak.flatten(table['matched_gen_vz'])

        BDT_pt = ak.flatten(table['emtfTrack_pt'])
        NN_pt  = ak.flatten(table['emtfTrack_pt_dxy'])
        NN_d0  = ak.flatten(table['emtfTrack_dxy'])
        matched_dR  = ak.flatten(table['emtfTrack_dR'])

        matched_gen_dR = gen_dR[gen_dR > -1]
        gen_dR_mask = (matched_gen_dR < 0.6)
        emtf_dR_mask = matched_dR < 5

        denom_mask = (abs(gen_vz) < 100) & (abs(gen_d0) < 100)
        num_mask = (abs(matched_gen_vz) < 100) & (abs(matched_gen_d0) < 100) 
        match_mask = (emtf_dR_mask & gen_dR_mask)

        denom_eta_mask = (abs(gen_eta) > 1.2) & (abs(gen_eta) < 2.5)
        denom_lo_eta_mask = (abs(gen_eta) > 1.2) & (abs(gen_eta) < 1.6)
        denom_md_eta_mask = (abs(gen_eta) > 1.6) & (abs(gen_eta) < 2.1)
        denom_hi_eta_mask = (abs(gen_eta) > 2.1) & (abs(gen_eta) < 2.5)

        num_eta_mask = (abs(matched_gen_eta) > 1.2) & (abs(matched_gen_eta) < 2.5)
        num_lo_eta_mask = (abs(matched_gen_eta) > 1.2) & (abs(matched_gen_eta) < 1.6)
        num_md_eta_mask = (abs(matched_gen_eta) > 1.6) & (abs(matched_gen_eta) < 2.1)
        num_hi_eta_mask = (abs(matched_gen_eta) > 2.1) & (abs(matched_gen_eta) < 2.5)

        denom_pt_mask = gen_pt > 15
        
        pt_mask_bdt = BDT_pt > 10
        pt_mask_nn = NN_pt > 10
        pt_mask_both = pt_mask_bdt | pt_mask_nn

        denom_d0_mask = gen_d0 < 20

        ##############################################################################
        ## pT efficiency
        den_pt_mask  = denom_mask & denom_eta_mask

        # no matching eff (eta cut but no matching cut)
        num_pt_bdt_mask  = num_mask & num_eta_mask & pt_mask_bdt
        num_pt_nn_mask   = num_mask & num_eta_mask & pt_mask_nn
        num_pt_both_mask = num_mask & num_eta_mask & pt_mask_both

        n_pt_num_bdt  = get_n(matched_gen_pt, ptbins, num_pt_bdt_mask)
        n_pt_num_nn   = get_n(matched_gen_pt, ptbins, num_pt_nn_mask)
        n_pt_num_both = get_n(matched_gen_pt, ptbins, num_pt_both_mask)
        n_pt_denom    = get_n(gen_pt, ptbins, den_pt_mask)

        # matching eff, 0.6, 5
        matched_num_pt_bdt_mask  = num_mask & num_eta_mask & match_mask & pt_mask_bdt
        matched_num_pt_nn_mask   = num_mask & num_eta_mask & match_mask & pt_mask_nn
        matched_num_pt_both_mask = num_mask & num_eta_mask & match_mask & pt_mask_both

        matched_n_pt_num_bdt  = get_n(matched_gen_pt, ptbins, matched_num_pt_bdt_mask)
        matched_n_pt_num_nn   = get_n(matched_gen_pt, ptbins, matched_num_pt_nn_mask)
        matched_n_pt_num_both = get_n(matched_gen_pt, ptbins, matched_num_pt_both_mask)
        matched_n_pt_denom    = get_n(gen_pt, ptbins, den_pt_mask)

        # no cuts (no eta)
        nocuts_num_pt_bdt_mask  = num_mask & pt_mask_bdt
        nocuts_num_pt_nn_mask   = num_mask & pt_mask_nn
        nocuts_num_pt_both_mask = num_mask & pt_mask_both

        nocuts_n_pt_num_bdt  = get_n(matched_gen_pt, ptbins, nocuts_num_pt_bdt_mask)
        nocuts_n_pt_num_nn   = get_n(matched_gen_pt, ptbins, nocuts_num_pt_nn_mask)
        nocuts_n_pt_num_both = get_n(matched_gen_pt, ptbins, nocuts_num_pt_both_mask)
        nocuts_n_pt_denom    = get_n(gen_pt, ptbins, denom_mask)


        ##############################################################################
        ## eta efficiency
        den_eta_mask  = denom_mask & denom_eta_mask & denom_pt_mask

        # no matching eff (eta cut but no matching cut)
        num_eta_bdt_mask  = num_mask & num_eta_mask & pt_mask_bdt
        num_eta_nn_mask   = num_mask & num_eta_mask & pt_mask_nn
        num_eta_both_mask = num_mask & num_eta_mask & pt_mask_both

        n_eta_num_bdt  = get_n(matched_gen_eta, etabins, num_eta_bdt_mask)
        n_eta_num_nn   = get_n(matched_gen_eta, etabins, num_eta_nn_mask)
        n_eta_num_both = get_n(matched_gen_eta, etabins, num_eta_both_mask)
        n_eta_denom    = get_n(gen_eta, etabins, den_eta_mask)

        # matching eff, 0.6, 5
        matched_num_eta_bdt_mask  = num_mask & num_eta_mask & match_mask & pt_mask_bdt
        matched_num_eta_nn_mask   = num_mask & num_eta_mask & match_mask & pt_mask_nn
        matched_num_eta_both_mask = num_mask & num_eta_mask & match_mask & pt_mask_both

        matched_n_eta_num_bdt  = get_n(matched_gen_eta, etabins, matched_num_eta_bdt_mask)
        matched_n_eta_num_nn   = get_n(matched_gen_eta, etabins, matched_num_eta_nn_mask)
        matched_n_eta_num_both = get_n(matched_gen_eta, etabins, matched_num_eta_both_mask)
        matched_n_eta_denom    = get_n(gen_eta, etabins, den_eta_mask)

        # no cuts (no eta)
        nocuts_num_eta_bdt_mask  = num_mask & pt_mask_bdt
        nocuts_num_eta_nn_mask   = num_mask & pt_mask_nn
        nocuts_num_eta_both_mask = num_mask & pt_mask_both

        nocuts_n_eta_num_bdt  = get_n(matched_gen_eta, ptbins, nocuts_num_eta_bdt_mask)
        nocuts_n_eta_num_nn   = get_n(matched_gen_eta, ptbins, nocuts_num_eta_nn_mask)
        nocuts_n_eta_num_both = get_n(matched_gen_eta, ptbins, nocuts_num_eta_both_mask)
        nocuts_n_eta_denom    = get_n(gen_eta, ptbins, denom_mask)

## pt plotting
plot_eff(n_pt_num_bdt, n_pt_num_nn, n_pt_num_both, n_pt_denom, ptbins, r'pt efficiency, all eta, no matching', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')

plot_eff(matched_n_pt_num_bdt, matched_n_pt_num_nn, matched_n_pt_num_both, matched_n_pt_denom, ptbins, r'pt efficiency, all eta', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')

plot_eff(nocuts_n_pt_num_bdt, nocuts_n_pt_num_nn, nocuts_n_pt_num_both, nocuts_n_pt_denom, ptbins, r'pt efficiency, no cuts', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')

## eta plotting
plot_eff(n_eta_num_bdt, n_eta_num_nn, n_eta_num_both, n_eta_denom, etabins, r'eta efficiency, all eta, no matching', r'$\eta^{*,\mathrm{gen}}$')

plot_eff(matched_n_eta_num_bdt, matched_n_eta_num_nn, matched_n_eta_num_both, matched_n_eta_denom, etabins, r'eta efficiency, all eta', r'$\eta^{*\mathrm{gen}}$')

plot_eff(nocuts_n_eta_num_bdt, nocuts_n_eta_num_nn, nocuts_n_eta_num_both, nocuts_n_eta_denom, etabins, r'eta efficiency, no cuts', r'$\eta{*^\mathrm{gen}}$')