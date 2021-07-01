"""
Project: EMTF / Displaced Muons
Author: Suzanne Rosenzweig

Denominator mask begins with |z0| < 100 cm and |d0| < 100 cm. 
"""

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

def start_plots():
    global pt_bdt_mask, pt_nn_mask, pt_comb_mask, pt_denom_mask, pt_bdt_denom_mask, eta_bdt_mask, eta_nn_mask, eta_comb_mask, eta_denom_mask, eta_bdt_denom_mask, d0_bdt_mask, d0_nn_mask, d0_comb_mask, d0_denom_mask, d0_bdt_denom_mask

    pt_bdt_mask, pt_nn_mask, pt_comb_mask, pt_denom_mask, pt_bdt_denom_mask = make_masks(pt_mask=False, eta_mask=True, d0_mask=False)
    eta_bdt_mask, eta_nn_mask, eta_comb_mask, eta_denom_mask, eta_bdt_denom_mask = make_masks(pt_mask=True, d0_mask=False)
    d0_bdt_mask, d0_nn_mask, d0_comb_mask, d0_denom_mask, d0_bdt_denom_mask = make_masks(pt_mask=True, eta_mask=True)

    global pt_low_eta_bdt_mask, pt_low_eta_nn_mask, pt_low_eta_comb_mask, pt_low_eta_denom_mask, pt_low_eta_bdt_denom_mask, pt_med_eta_bdt_mask, pt_med_eta_nn_mask, pt_med_eta_comb_mask, pt_med_eta_denom_mask, pt_med_eta_bdt_denom_mask, pt_high_eta_bdt_mask, pt_high_eta_nn_mask, pt_high_eta_comb_mask, pt_high_eta_denom_mask, pt_high_eta_bdt_denom_mask

    pt_low_eta_bdt_mask, pt_low_eta_nn_mask, pt_low_eta_comb_mask, pt_low_eta_denom_mask, pt_low_eta_bdt_denom_mask = make_masks(pt_mask=False, d0_mask=True, low_eta=True)
    pt_med_eta_bdt_mask, pt_med_eta_nn_mask, pt_med_eta_comb_mask, pt_med_eta_denom_mask, pt_med_eta_bdt_denom_mask = make_masks(pt_mask=False, d0_mask=True, med_eta=True)
    pt_high_eta_bdt_mask, pt_high_eta_nn_mask, pt_high_eta_comb_mask, pt_high_eta_denom_mask, pt_high_eta_bdt_denom_mask = make_masks(pt_mask=False, d0_mask=True, high_eta=True)

    global d0_low_eta_bdt_mask, d0_low_eta_nn_mask, d0_low_eta_comb_mask, d0_low_eta_denom_mask, d0_low_eta_bdt_denom_mask, d0_med_eta_bdt_mask, d0_med_eta_nn_mask, d0_med_eta_comb_mask, d0_med_eta_denom_mask, d0_med_eta_bdt_denom_mask, d0_high_eta_bdt_mask, d0_high_eta_nn_mask, d0_high_eta_comb_mask, d0_high_eta_denom_mask, d0_high_eta_bdt_denom_mask

    d0_low_eta_bdt_mask, d0_low_eta_nn_mask, d0_low_eta_comb_mask, d0_low_eta_denom_mask, d0_low_eta_bdt_denom_mask = make_masks(pt_mask=True, d0_mask=True, low_eta=True)
    d0_med_eta_bdt_mask, d0_med_eta_nn_mask, d0_med_eta_comb_mask, d0_med_eta_denom_mask, d0_med_eta_bdt_denom_mask = make_masks(pt_mask=True, d0_mask=True, med_eta=True)
    d0_high_eta_bdt_mask, d0_high_eta_nn_mask, d0_high_eta_comb_mask, d0_high_eta_denom_mask, d0_high_eta_bdt_denom_mask = make_masks(pt_mask=True, d0_mask=True, high_eta=True)

def make_masks(pt_mask=True, eta_mask=False, d0_mask=False, low_eta=False, med_eta=False, high_eta=False):
    global gen_eta, gen_d0, gen_vz, gen_pt, gen_d0, NN_pt, BDT_pt, matched_gen_eta, matched_gen_pt, matched_gen_d0, matched_gen_vz
    
    gen_eta_mask = (abs(gen_eta) > 1.24) & (abs(gen_eta) < 2.5)
    gen_z0_d0_mask =  (abs(gen_d0) < 100) & (abs(gen_vz) < 100) # cm, cm
    gen_pt_mask = gen_pt > 15 # GeV
    gen_d0_mask = abs(gen_d0) < 15
    low_eta_mask = (abs(gen_eta) > 1.24) & (abs(gen_eta) < 1.6)
    med_eta_mask = (abs(gen_eta) > 1.6) & (abs(gen_eta) < 2.1)
    high_eta_mask = (abs(gen_eta) > 2.1) & (abs(gen_eta) < 2.5)

    matched_gen_z0_d0_mask =  (abs(matched_gen_d0) < 100) & (abs(matched_gen_vz) < 100) # cm, cm
    matched_gen_pt_mask = matched_gen_pt > 15 # GeV
    matched_gen_d0_mask = abs(matched_gen_d0) < 15
    matched_gen_eta_mask = (abs(matched_gen_eta) > 1.24) & (abs(matched_gen_eta) < 2.5)
    matched_low_eta_mask = (abs(matched_gen_eta) > 1.24) & (abs(matched_gen_eta) < 1.6)
    matched_med_eta_mask = (abs(matched_gen_eta) > 1.6) & (abs(matched_gen_eta) < 2.1)
    matched_high_eta_mask = (abs(matched_gen_eta) > 2.1) & (abs(matched_gen_eta) < 2.5)

    BDT_mask = BDT_pt > 10 # GeV
    NN_mask  = NN_pt > 10 # GeV

    bdt_mask = BDT_mask &  matched_gen_z0_d0_mask 
    nn_mask = NN_mask & matched_gen_z0_d0_mask
    comb_mask = (BDT_mask | NN_mask) & matched_gen_z0_d0_mask
    denom_mask = matched_gen_z0_d0_mask ##
    bdt_denom_mask = denom_mask

    if low_eta:
        bdt_mask = matched_low_eta_mask & bdt_mask
        nn_mask = matched_low_eta_mask & nn_mask
        comb_mask = matched_low_eta_mask & comb_mask
        denom_mask = matched_low_eta_mask & denom_mask ##
    if med_eta:
        bdt_mask = matched_med_eta_mask & bdt_mask
        nn_mask = matched_med_eta_mask & nn_mask
        comb_mask = matched_med_eta_mask & comb_mask 
        denom_mask = matched_med_eta_mask & denom_mask ##
    if high_eta:
        bdt_mask = matched_high_eta_mask & bdt_mask
        nn_mask = matched_high_eta_mask & nn_mask
        comb_mask = matched_high_eta_mask & comb_mask
        denom_mask = matched_gen_eta_mask & denom_mask ##

    if eta_mask: 
        bdt_mask = matched_gen_eta_mask & bdt_mask
        nn_mask = matched_gen_eta_mask & nn_mask
        comb_mask = matched_gen_eta_mask & comb_mask
        denom_mask = matched_gen_eta_mask & denom_mask ##

    if pt_mask: 
        bdt_mask = matched_gen_pt_mask & bdt_mask
        nn_mask = matched_gen_pt_mask & nn_mask
        comb_mask = matched_gen_pt_mask & comb_mask
        denom_mask = matched_gen_pt_mask & denom_mask ##

    if d0_mask: 
        bdt_mask = matched_gen_d0_mask & bdt_mask

    return bdt_mask, nn_mask, comb_mask, denom_mask, bdt_denom_mask

def get_n_eff_err(bins, arr, comb_mask, bdt_mask, nn_mask, denom_mask, comb_n, bdt_n, nn_n, den_n):
    n_comb, n_edges = np.histogram(ak.to_list(arr[comb_mask]), bins=bins)
    n_bdt, n_edges = np.histogram(ak.to_list(arr[bdt_mask]), bins=bins)
    n_nn, n_edges = np.histogram(ak.to_list(arr[nn_mask]), bins=bins)
    n_den, d_edges = np.histogram(ak.to_list(arr[denom_mask]), bins=bins)

    n_comb += comb_n
    n_bdt += bdt_n
    n_nn += nn_n
    n_den += den_n

    eff_comb = np.where(n_den == 0, 0, n_comb/n_den)
    eff_bdt = np.where(n_den == 0, 0, n_bdt/n_den)
    eff_nn = np.where(n_den == 0, 0, n_nn/n_den)

    err_comb = np.sqrt(eff_comb*(1-eff_comb) / n_den)
    err_bdt = np.sqrt(eff_bdt*(1-eff_bdt) / n_den)
    err_nn = np.sqrt(eff_nn*(1-eff_nn) / n_den)

    return eff_comb, eff_bdt, eff_nn, n_comb, n_bdt, n_nn, n_den, err_comb, err_bdt, err_nn

def make_fig_ax(title, bins, x, xlabel, ylabel, eff_comb, eff_bdt, eff_nn, err_comb, err_bdt, err_nn):

    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_ylim(0,1.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    for i,w in enumerate(eff_comb):
        ax.plot([bins[i], bins[i+1]], [w,w], color='blue', label='Combined')
        ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='blue')
    for i,w in enumerate(eff_bdt):
        ax.plot([bins[i], bins[i+1]], [w,w], color='black', label='BDT')
        ax.plot([x[i], x[i]], [w-err_bdt[i], w+err_bdt[i]], color='black')
    for i,w in enumerate(eff_nn):
        ax.plot([bins[i], bins[i+1]], [w,w], color='red',  label='NN')
        ax.plot([x[i], x[i]], [w-err_nn[i], w+err_nn[i]], color='red')

    custom_lines = [Line2D([0], [0], color='blue', lw=2),
                    Line2D([0], [0], color='black', lw=2),
                    Line2D([0], [0], color='red', lw=2)]

    ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
    fig.savefig(f'{folder}{title}.pdf', bbox_inches=None)

ptbins = np.linspace(0, 60, 30)
x_pt = (ptbins[1:] + ptbins[:-1]) / 2

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

filelist = [filename1, filename2, filename3, filename4]

def make_dict():
    return {
        'bdt':{'low_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'med_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'high_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'all_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}}, 

        'nn':{'low_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'med_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'high_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'all_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}}, 
        
        'both':{'low_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'med_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'high_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}, 
            'all_eta':{'pt':{'eff':0, 'n':0, 'err':0}, 
                        'eta':{'eff':0, 'n':0, 'err':0}, 
                        'd0':{'eff':0, 'n':0, 'err':0}}},
                        
        'denom':{'low_eta':{'pt':{'n':0}, 
                        'eta':{'n':0}, 
                        'd0':{'n':0}}, 
            'med_eta':{'pt':{'n':0}, 
                        'eta':{'n':0}, 
                        'd0':{'n':0}}, 
            'high_eta':{'pt':{'n':0}, 
                        'eta':{'n':0}, 
                        'd0':{'n':0}}, 
            'all_eta':{'pt':{'n':0}, 
                        'eta':{'n':0}, 
                        'd0':{'n':0}}}}

dm_dict = make_dict()

for filename in filelist:
    info("Loading file")
    info(filename)

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

        matched_gen_pt      = ak.flatten(table['matched_gen_pt'])
        matched_gen_phi     = ak.flatten(table['matched_gen_phi'])
        matched_gen_eta     = ak.flatten(table['matched_gen_etaStar'])
        matched_gen_phiStar = ak.flatten(table['matched_gen_phiStar'])
        # matched_gen_etaStar = ak.flatten(table['matched_gen_etaStar'])
        matched_gen_q       = ak.flatten(table['matched_gen_q'])
        matched_gen_d0      = ak.flatten(table['matched_gen_d0'])
        matched_gen_vz      = ak.flatten(table['matched_gen_vz'])

        BDT_pt = ak.flatten(table['emtfTrack_pt'])
        NN_pt  = ak.flatten(table['emtfTrack_pt_dxy'])
        NN_d0  = ak.flatten(table['emtfTrack_dxy'])

        info("Making masks.")
        start_plots()

        info("Calculating low eta efficiency.")
        comb_n = dm_dict['both']['low_eta']['pt']['n']
        bdt_n = dm_dict['bdt']['low_eta']['pt']['n']
        nn_n = dm_dict['nn']['low_eta']['pt']['n']
        den_n = dm_dict['denom']['low_eta']['pt']['n']

        eff_comb, eff_bdt, eff_nn, n_comb, n_bdt, n_nn, n_den, err_comb, err_bdt, err_nn = get_n_eff_err(ptbins, matched_gen_pt, gen_pt, pt_low_eta_comb_mask, pt_low_eta_bdt_mask, pt_low_eta_nn_mask, pt_low_eta_denom_mask, comb_n, bdt_n, nn_n, den_n)
        
        dm_dict['both']['low_eta']['pt']['eff'] = eff_comb
        dm_dict['bdt']['low_eta']['pt']['eff'] = eff_bdt
        dm_dict['nn']['low_eta']['pt']['eff'] = eff_nn

        dm_dict['both']['low_eta']['pt']['n'] += n_comb
        dm_dict['bdt']['low_eta']['pt']['n'] += n_bdt
        dm_dict['nn']['low_eta']['pt']['n'] += n_nn
        dm_dict['denom']['low_eta']['pt']['n'] += n_den

        dm_dict['both']['low_eta']['pt']['err'] = err_comb
        dm_dict['bdt']['low_eta']['pt']['err'] = err_bdt
        dm_dict['nn']['low_eta']['pt']['err'] = err_nn

        print(eff_comb)
        print(n_comb)
        print(n_den)


        info("Calculating med eta efficiency.")
        comb_n = dm_dict['both']['med_eta']['pt']['n']
        bdt_n = dm_dict['bdt']['med_eta']['pt']['n']
        nn_n = dm_dict['nn']['med_eta']['pt']['n']
        den_n = dm_dict['denom']['med_eta']['pt']['n']

        eff_comb, eff_bdt, eff_nn, n_comb, n_bdt, n_nn, n_den, err_comb, err_bdt, err_nn = get_n_eff_err(ptbins, gen_pt, pt_med_eta_comb_mask, pt_med_eta_bdt_mask, pt_med_eta_nn_mask, pt_med_eta_denom_mask, comb_n, bdt_n, nn_n, den_n)
        
        dm_dict['both']['med_eta']['pt']['eff'] = eff_comb
        dm_dict['bdt']['med_eta']['pt']['eff'] = eff_bdt
        dm_dict['nn']['med_eta']['pt']['eff'] = eff_nn

        dm_dict['both']['med_eta']['pt']['n'] += n_comb
        dm_dict['bdt']['med_eta']['pt']['n'] += n_bdt
        dm_dict['nn']['med_eta']['pt']['n'] += n_nn

        dm_dict['denom']['med_eta']['pt']['n'] += n_den

        dm_dict['both']['med_eta']['pt']['err'] = err_comb
        dm_dict['bdt']['med_eta']['pt']['err'] = err_bdt
        dm_dict['nn']['med_eta']['pt']['err'] = err_nn


        info("Calculating high eta efficiency.")
        comb_n = dm_dict['both']['high_eta']['pt']['n']
        bdt_n = dm_dict['bdt']['high_eta']['pt']['n']
        nn_n = dm_dict['nn']['high_eta']['pt']['n']
        den_n = dm_dict['denom']['high_eta']['pt']['n']

        eff_comb, eff_bdt, eff_nn, n_comb, n_bdt, n_nn, n_den, err_comb, err_bdt, err_nn = get_n_eff_err(ptbins, gen_pt, pt_high_eta_comb_mask, pt_high_eta_bdt_mask, pt_high_eta_nn_mask, pt_high_eta_denom_mask, comb_n, bdt_n, nn_n, den_n)
        
        dm_dict['both']['high_eta']['pt']['eff'] = eff_comb
        dm_dict['bdt']['high_eta']['pt']['eff'] = eff_bdt
        dm_dict['nn']['high_eta']['pt']['eff'] = eff_nn

        dm_dict['both']['high_eta']['pt']['n'] += n_comb
        dm_dict['bdt']['high_eta']['pt']['n'] += n_bdt
        dm_dict['nn']['high_eta']['pt']['n'] += n_nn

        dm_dict['denom']['high_eta']['pt']['n'] += n_den

        dm_dict['both']['high_eta']['pt']['err'] = err_comb
        dm_dict['bdt']['high_eta']['pt']['err'] = err_bdt
        dm_dict['nn']['high_eta']['pt']['err'] = err_nn



        info("Calculating all eta efficiency.")
        comb_n = dm_dict['both']['all_eta']['pt']['n']
        bdt_n = dm_dict['bdt']['all_eta']['pt']['n']
        nn_n = dm_dict['nn']['all_eta']['pt']['n']
        den_n = dm_dict['denom']['all_eta']['pt']['n']

        eff_comb, eff_bdt, eff_nn, n_comb, n_bdt, n_nn, n_den, err_comb, err_bdt, err_nn = get_n_eff_err(ptbins, gen_pt, pt_comb_mask, pt_bdt_mask, pt_nn_mask, pt_denom_mask, comb_n, bdt_n, nn_n, den_n)
        
        dm_dict['both']['all_eta']['pt']['eff'] = eff_comb
        dm_dict['bdt']['all_eta']['pt']['eff'] = eff_bdt
        dm_dict['nn']['all_eta']['pt']['eff'] = eff_nn

        dm_dict['both']['all_eta']['pt']['n'] += n_comb
        dm_dict['bdt']['all_eta']['pt']['n'] += n_bdt
        dm_dict['nn']['all_eta']['pt']['n'] += n_nn

        dm_dict['denom']['all_eta']['pt']['n'] += n_den

        dm_dict['both']['all_eta']['pt']['err'] = err_comb
        dm_dict['bdt']['all_eta']['pt']['err'] = err_bdt
        dm_dict['nn']['all_eta']['pt']['err'] = err_nn
        

ylabel = 'L1T Efficiency'


eff_comb = dm_dict['both']['low_eta']['pt']['eff']
eff_bdt = dm_dict['bdt']['low_eta']['pt']['eff']
eff_nn = dm_dict['nn']['low_eta']['pt']['eff']

err_comb = dm_dict['both']['low_eta']['pt']['err']
err_bdt = dm_dict['bdt']['low_eta']['pt']['err']
err_nn = dm_dict['nn']['low_eta']['pt']['err']

xlabel = r'gen muon $p_T$ [GeV]'
title = 'pt_efficiency_low_eta'
make_fig_ax(title, ptbins, x_pt, xlabel, ylabel, eff_comb, eff_bdt, eff_nn, err_comb, err_bdt, err_nn)


eff_comb = dm_dict['both']['all_eta']['pt']['eff']
eff_bdt = dm_dict['bdt']['all_eta']['pt']['eff']
eff_nn = dm_dict['nn']['all_eta']['pt']['eff']

err_comb = dm_dict['both']['all_eta']['pt']['err']
err_bdt = dm_dict['bdt']['all_eta']['pt']['err']
err_nn = dm_dict['nn']['all_eta']['pt']['err']

xlabel = r'gen muon $p_T$ [GeV]'
title = 'pt_efficiency'
make_fig_ax(title, ptbins, x_pt, xlabel, ylabel, eff_comb, eff_bdt, eff_nn, err_comb, err_bdt, err_nn)





#############################################################
# pT from EMTF

# pt_n_comb, n_edges = np.histogram(ak.to_list(gen_pt[pt_med_eta_comb_mask]), bins=ptbins)
# pt_n_bdt, n_edges = np.histogram(ak.to_list(gen_pt[pt_med_eta_bdt_mask]), bins=ptbins)
# pt_n_nn, n_edges = np.histogram(ak.to_list(gen_pt[pt_med_eta_nn_mask]), bins=ptbins)
# pt_n_den, d_edges = np.histogram(ak.to_list(gen_pt[pt_med_eta_denom_mask]), bins=ptbins)
# pt_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_pt[pt_med_eta_bdt_denom_mask]), bins=ptbins)

# pt_eff_comb = np.where(pt_n_den == 0, 0, pt_n_comb/pt_n_den)
# pt_eff_bdt = np.where(pt_n_den == 0, 0, pt_n_bdt/pt_n_den_bdt)
# pt_eff_nn = np.where(pt_n_den == 0, 0, pt_n_nn/pt_n_den)

# ### pT

# x = (n_edges[1:] + n_edges[:-1]) / 2
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_ylim(0,1.05)
# ax.set_xlabel(r'gen muon $p_T$ [GeV]')
# ax.set_ylabel('L1T Efficiency')

# pt_err_comb = np.sqrt(pt_eff_comb*(1-pt_eff_comb) / pt_n_den)
# pt_err_bdt = np.sqrt(pt_eff_bdt*(1-pt_eff_bdt) / pt_n_den_bdt)
# pt_err_nn = np.sqrt(pt_eff_nn*(1-pt_eff_nn) / pt_n_den)

# for i,w in enumerate(pt_eff_comb):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='blue', label='Combined')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='blue')
# # for i,w in enumerate(eff_bdt):
# #     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
# #     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_bdt):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_nn):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='red',  label='NN')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='red')

# from matplotlib.lines import Line2D
# custom_lines = [Line2D([0], [0], color='blue', lw=2),
#                 Line2D([0], [0], color='black', lw=2),
#                 Line2D([0], [0], color='red', lw=2)]

# ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
# fig.savefig(f'{folder}pt_efficiency_med_eta.pdf', bbox_inches=None)

# print("Checkpoint 5")

# ### pT



# ptbins = np.linspace(0, 60, 30)
# pt_n_comb, n_edges = np.histogram(ak.to_list(gen_pt[pt_high_eta_comb_mask]), bins=ptbins)
# pt_n_bdt, n_edges = np.histogram(ak.to_list(gen_pt[pt_high_eta_bdt_mask]), bins=ptbins)
# pt_n_nn, n_edges = np.histogram(ak.to_list(gen_pt[pt_high_eta_nn_mask]), bins=ptbins)
# pt_n_den, d_edges = np.histogram(ak.to_list(gen_pt[pt_high_eta_denom_mask]), bins=ptbins)
# pt_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_pt[pt_high_eta_bdt_denom_mask]), bins=ptbins)

# pt_eff_comb = np.where(pt_n_den == 0, 0, pt_n_comb/pt_n_den)
# pt_eff_bdt = np.where(pt_n_den == 0, 0, pt_n_bdt/pt_n_den_bdt)
# pt_eff_nn = np.where(pt_n_den == 0, 0, pt_n_nn/pt_n_den)

# x = (n_edges[1:] + n_edges[:-1]) / 2
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_ylim(0,1.05)
# ax.set_xlabel(r'gen muon $p_T$ [GeV]')
# ax.set_ylabel('L1T Efficiency')

# pt_err_comb = np.sqrt(pt_eff_comb*(1-pt_eff_comb) / pt_n_den)
# pt_err_bdt = np.sqrt(pt_eff_bdt*(1-pt_eff_bdt) / pt_n_den_bdt)
# pt_err_nn = np.sqrt(pt_eff_nn*(1-pt_eff_nn) / pt_n_den)

# for i,w in enumerate(pt_eff_comb):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='blue', label='Combined')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='blue')
# # for i,w in enumerate(eff_bdt):
# #     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
# #     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_bdt):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_nn):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='red',  label='NN')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='red')

# from matplotlib.lines import Line2D
# custom_lines = [Line2D([0], [0], color='blue', lw=2),
#                 Line2D([0], [0], color='black', lw=2),
#                 Line2D([0], [0], color='red', lw=2)]

# ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
# fig.savefig(f'{folder}pt_efficiency_high_eta.pdf', bbox_inches=None)



# print("Checkpoint 6")


# ptbins = np.linspace(0, 60, 30)
# pt_n_comb, n_edges = np.histogram(ak.to_list(gen_pt[pt_comb_mask]), bins=ptbins)
# pt_n_bdt, n_edges = np.histogram(ak.to_list(gen_pt[pt_bdt_mask]), bins=ptbins)
# pt_n_nn, n_edges = np.histogram(ak.to_list(gen_pt[pt_nn_mask]), bins=ptbins)
# pt_n_den, d_edges = np.histogram(ak.to_list(gen_pt[pt_denom_mask]), bins=ptbins)
# pt_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_pt[pt_bdt_denom_mask]), bins=ptbins)

# pt_eff_comb = np.where(pt_n_den == 0, 0, pt_n_comb/pt_n_den)
# pt_eff_bdt = np.where(pt_n_den == 0, 0, pt_n_bdt/pt_n_den_bdt)
# pt_eff_nn = np.where(pt_n_den == 0, 0, pt_n_nn/pt_n_den)

# ### pT

# x = (n_edges[1:] + n_edges[:-1]) / 2
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_ylim(0,1.05)
# ax.set_xlabel(r'gen muon $p_T$ [GeV]')
# ax.set_ylabel('L1T Efficiency')

# pt_err_comb = np.sqrt(pt_eff_comb*(1-pt_eff_comb) / pt_n_den)
# pt_err_bdt = np.sqrt(pt_eff_bdt*(1-pt_eff_bdt) / pt_n_den_bdt)
# pt_err_nn = np.sqrt(pt_eff_nn*(1-pt_eff_nn) / pt_n_den)

# for i,w in enumerate(pt_eff_comb):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='blue', label='Combined')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='blue')
# # for i,w in enumerate(eff_bdt):
# #     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
# #     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_bdt):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_nn):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='red',  label='NN')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='red')

# from matplotlib.lines import Line2D
# custom_lines = [Line2D([0], [0], color='blue', lw=2),
#                 Line2D([0], [0], color='black', lw=2),
#                 Line2D([0], [0], color='red', lw=2)]

# ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
# fig.savefig(f'{folder}pt_efficiency.pdf', bbox_inches=None)


# print("Checkpoint 7")

# ### q/pT

# ptbins = np.linspace(-.5, .5, 30)
# pt_n_comb, n_edges = np.histogram(ak.to_list(gen_q[pt_comb_mask]/gen_pt[pt_comb_mask]), bins=ptbins)
# pt_n_bdt, n_edges = np.histogram(ak.to_list(gen_q[pt_bdt_mask]/gen_pt[pt_bdt_mask]), bins=ptbins)
# pt_n_nn, n_edges = np.histogram(ak.to_list(gen_q[pt_nn_mask]/gen_pt[pt_nn_mask]), bins=ptbins)
# pt_n_den, d_edges = np.histogram(ak.to_list(gen_q[pt_denom_mask]/gen_pt[pt_denom_mask]), bins=ptbins)
# pt_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_q[pt_bdt_denom_mask]/gen_pt[pt_bdt_denom_mask]), bins=ptbins)

# pt_eff_comb = np.where(pt_n_den == 0, 0, pt_n_comb/pt_n_den)
# pt_eff_bdt = np.where(pt_n_den == 0, 0, pt_n_bdt/pt_n_den_bdt)
# pt_eff_nn = np.where(pt_n_den == 0, 0, pt_n_nn/pt_n_den)

# x = (n_edges[1:] + n_edges[:-1]) / 2

# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_ylim(0,1.05)
# ax.set_xlabel(r'gen muon $p_T$ [GeV]')
# ax.set_ylabel('L1T Efficiency')

# pt_err_comb = np.sqrt(pt_eff_comb*(1-pt_eff_comb) / pt_n_den)
# pt_err_bdt = np.sqrt(pt_eff_bdt*(1-pt_eff_bdt) / pt_n_den_bdt)
# pt_err_nn = np.sqrt(pt_eff_nn*(1-pt_eff_nn) / pt_n_den)

# for i,w in enumerate(pt_eff_comb):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='blue', label='Combined')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='blue')
# # for i,w in enumerate(eff_bdt):
# #     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
# #     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_bdt):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='black')
# for i,w in enumerate(pt_eff_nn):
#     ax.plot([ptbins[i], ptbins[i+1]], [w,w], color='red',  label='NN')
#     ax.plot([x[i], x[i]], [w-pt_err_comb[i], w+pt_err_comb[i]], color='red')

# from matplotlib.lines import Line2D
# custom_lines = [Line2D([0], [0], color='blue', lw=2),
#                 Line2D([0], [0], color='black', lw=2),
#                 Line2D([0], [0], color='red', lw=2)]

# ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
# fig.savefig(f'{folder}invpt_efficiency.pdf', bbox_inches=None)


# print("Checkpoint 8")

# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'$p_T^\mathrm{NN} - p_T^\mathrm{GEN}$ [GeV]')
# ax.set_ylabel('Count')
# hist(ax, ak.to_list(NN_pt[pt_nn_mask]-gen_pt[pt_nn_mask]), bins=np.linspace(-60, 60, 60))
# fig.savefig(f'{folder}pt_nn_difference.pdf', bbox_inches=None)


# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'$\Delta p_T^\mathrm{NN}/p_T^\mathrm{GEN}$')
# ax.set_ylabel('Count')
# hist(ax, (NN_pt[pt_nn_mask]-gen_pt[pt_nn_mask])/gen_pt[pt_nn_mask], bins=np.linspace(-5, 15, 60))
# fig.savefig(f'{folder}pt_nn_diff_norm.pdf', bbox_inches=None)



# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'$p_T^\mathrm{BDT} - p_T^\mathrm{GEN}$ [GeV]')
# ax.set_ylabel('Count')
# hist(ax, BDT_pt[pt_bdt_mask]-gen_pt[pt_bdt_mask], bins=np.linspace(-60, 60, 60))
# fig.savefig(f'{folder}pt_bdt_difference.pdf', bbox_inches=None)


# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'$\Delta p_T^\mathrm{BDT}/p_T^\mathrm{GEN}$')
# ax.set_ylabel('Count')
# hist(ax, (BDT_pt[pt_bdt_mask]-gen_pt[pt_bdt_mask])/gen_pt[pt_bdt_mask], bins=np.linspace(-5, 15, 60))
# fig.savefig(f'{folder}pt_bdt_diff_norm.pdf', bbox_inches=None)



# ###                     ETA PLOTS


# etabins = [-2.5, -2.1, -1.6, -1.24]
# for eta in etabins[::-1]:
#     etabins.append(-eta)
# eta_n_comb, n_edges = np.histogram(ak.to_list(gen_eta[eta_comb_mask]), bins=etabins)
# eta_n_bdt, n_edges = np.histogram(ak.to_list(gen_eta[eta_bdt_mask]), bins=etabins)
# eta_n_nn, n_edges = np.histogram(ak.to_list(gen_eta[eta_nn_mask]), bins=etabins)
# eta_n_den, d_edges = np.histogram(ak.to_list(gen_eta[eta_denom_mask]), bins=etabins)
# eta_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_eta[eta_bdt_denom_mask]), bins=etabins)

# eta_eff_comb = np.where(eta_n_den == 0, 0, eta_n_comb/eta_n_den)
# eta_eff_bdt = np.where(eta_n_den == 0, 0, eta_n_bdt/eta_n_den_bdt)
# eta_eff_nn = np.where(eta_n_den == 0, 0, eta_n_nn/eta_n_den)

# x = (n_edges[1:] + n_edges[:-1]) / 2
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_ylim(0,1.05)
# ax.set_xlabel(r'gen muon $\eta$')
# ax.set_ylabel('L1T Efficiency')

# eta_err_comb = np.sqrt(eta_eff_comb*(1-eta_eff_comb) / eta_n_den)
# eta_err_bdt = np.sqrt(eta_eff_bdt*(1-eta_eff_bdt) / eta_n_den_bdt)
# eta_err_nn = np.sqrt(eta_eff_nn*(1-eta_eff_nn) / eta_n_den)

# for i,w in enumerate(eta_eff_comb):
#     if i == 3: continue
#     ax.plot([etabins[i], etabins[i+1]], [w,w], color='blue', label='Combined')
#     ax.plot([x[i], x[i]], [w-eta_err_comb[i], w+eta_err_comb[i]], color='blue')
# # for i,w in enumerate(eff_bdt):
# #     ax.plot([etabins[i], etabins[i+1]], [w,w], color='black', label='BDT')
# #     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
# for i,w in enumerate(eta_eff_bdt):
#     if i == 3: continue
#     ax.plot([etabins[i], etabins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-eta_err_comb[i], w+eta_err_comb[i]], color='black')
# for i,w in enumerate(eta_eff_nn):
#     if i == 3: continue
#     ax.plot([etabins[i], etabins[i+1]], [w,w], color='red',  label='NN')
#     ax.plot([x[i], x[i]], [w-eta_err_comb[i], w+eta_err_comb[i]], color='red')

# from matplotlib.lines import Line2D
# custom_lines = [Line2D([0], [0], color='blue', lw=2),
#                 Line2D([0], [0], color='black', lw=2),
#                 Line2D([0], [0], color='red', lw=2)]

# ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
# fig.savefig(f'{folder}eta_efficiency.pdf', bbox_inches=None)

# ###                     D0 PLOTS


# d0bins = np.linspace(0,100,25)
# d0_n_comb, n_edges = np.histogram(ak.to_list(gen_d0[d0_comb_mask]), bins=d0bins)
# d0_n_bdt, n_edges = np.histogram(ak.to_list(gen_d0[d0_bdt_mask]), bins=d0bins)
# d0_n_nn, n_edges = np.histogram(ak.to_list(gen_d0[d0_nn_mask]), bins=d0bins)
# d0_n_den, d_edges = np.histogram(ak.to_list(gen_d0[d0_denom_mask]), bins=d0bins)
# d0_n_den_bdt, d_edges = np.histogram(ak.to_list(gen_d0[d0_bdt_denom_mask]), bins=d0bins)

# d0_eff_comb = np.where(d0_n_den == 0, 0, d0_n_comb/d0_n_den)
# d0_eff_bdt = np.where(d0_n_den == 0, 0, d0_n_bdt/d0_n_den_bdt)
# d0_eff_nn = np.where(d0_n_den == 0, 0, d0_n_nn/d0_n_den)

# x = (n_edges[1:] + n_edges[:-1]) / 2
# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_ylim(0,1.05)
# ax.set_xlabel(r'gen muon $d_0$ [cm]')
# ax.set_ylabel('L1T Efficiency')

# d0_err_comb = np.sqrt(d0_eff_comb*(1-d0_eff_comb) / d0_n_den)
# d0_err_bdt = np.sqrt(d0_eff_bdt*(1-d0_eff_bdt) / d0_n_den_bdt)
# d0_err_nn = np.sqrt(d0_eff_nn*(1-d0_eff_nn) / d0_n_den)

# for i,w in enumerate(d0_eff_comb):
#     ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='blue', label='Combined')
#     ax.plot([x[i], x[i]], [w-d0_err_comb[i], w+d0_err_comb[i]], color='blue')
# # for i,w in enumerate(eff_bdt):
# #     ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='black', label='BDT')
# #     ax.plot([x[i], x[i]], [w-err_comb[i], w+err_comb[i]], color='black')
# for i,w in enumerate(d0_eff_bdt):
#     ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='black', label='BDT')
#     ax.plot([x[i], x[i]], [w-d0_err_comb[i], w+d0_err_comb[i]], color='black')
# for i,w in enumerate(d0_eff_nn):
#     ax.plot([d0bins[i], d0bins[i+1]], [w,w], color='red',  label='NN')
#     ax.plot([x[i], x[i]], [w-d0_err_comb[i], w+d0_err_comb[i]], color='red')

# from matplotlib.lines import Line2D
# custom_lines = [Line2D([0], [0], color='blue', lw=2),
#                 Line2D([0], [0], color='black', lw=2),
#                 Line2D([0], [0], color='red', lw=2)]

# ax.legend(custom_lines, ['Combined', 'BDT', 'NN'])
# fig.savefig(f'{folder}d0_efficiency.pdf', bbox_inches=None)


# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'$d_0^\mathrm{NN} - d_0^\mathrm{GEN}$ [cm]')
# ax.set_xlim(-250,250)
# ax.set_ylabel('Count')
# hist(ax, ak.to_list(NN_d0[d0_nn_mask]-gen_d0[d0_nn_mask]), bins=np.linspace(-60, 60, 60))
# fig.savefig(f'{folder}d0_nn_difference.pdf', bbox_inches=None)


# fig, ax = plt.subplots(figsize=(10,10))
# ax.set_xlabel(r'$\Delta d_0^\mathrm{NN}/d_0^\mathrm{GEN}$')
# ax.set_xlim(-10,10)
# ax.set_ylabel('Count')
# hist(ax, (NN_d0[d0_nn_mask]-gen_d0[d0_nn_mask])/gen_d0[d0_nn_mask], bins=np.linspace(-5, 15, 60))
# fig.savefig(f'{folder}d0_nn_diff_norm.pdf', bbox_inches=None)
