from configparser import ConfigParser
import awkward as ak
from icecream import ic
import numpy as np
np.seterr(all='ignore')
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import uproot
import sys

def get_n(arr, bins, mask):
    n, edges = np.histogram(ak.to_list(arr[mask]), bins=bins)
    return n

def plot_eff(n_bdt, n_nn, n_both, n_den, bins, title, xlabel, etaskip=False):
    global folder
    eff_both = np.where(n_den == 0, 0, n_both/n_den)
    eff_bdt = np.where(n_den == 0, 0, n_bdt/n_den)
    eff_nn = np.where(n_den == 0, 0, n_nn/n_den)

    # ic(title, eff_both, eff_bdt, eff_nn)

    err_both = np.sqrt(eff_both*(1-eff_both) / n_den)
    err_bdt = np.sqrt(eff_bdt*(1-eff_bdt) / n_den)
    err_nn = np.sqrt(eff_nn*(1-eff_nn) / n_den)

    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_ylim(0,1.05)

    x = (bins[1:] + bins[:-1]) / 2

    for i,w in enumerate(eff_both):
        if etaskip and i == 3: continue
        ax.plot([bins[i], bins[i+1]], [w,w], color='blue') # BDT or NN
        ax.plot([x[i], x[i]], [w-err_both[i], w+err_both[i]], color='blue')
    for i,w in enumerate(eff_bdt):
        if etaskip and i == 3: continue
        ax.plot([bins[i], bins[i+1]], [w,w], color='black') # BDT
        ax.plot([x[i], x[i]], [w-err_bdt[i], w+err_bdt[i]], color='black')
    for i,w in enumerate(eff_nn):
        if etaskip and i == 3: continue
        ax.plot([bins[i], bins[i+1]], [w,w], color='red') # NN
        ax.plot([x[i], x[i]], [w-err_nn[i], w+err_nn[i]], color='red')

    custom_lines = [Line2D([0], [0], color='blue', lw=2),
                    Line2D([0], [0], color='black', lw=2),
                    Line2D([0], [0], color='red', lw=2)]

    ax.legend(custom_lines, [r'L1 $p_T^\mathrm{BDT} > 10 \;\mathrm{GeV} \;||\; L1 p_T^\mathrm{NN} > 10 \;\mathrm{GeV}$',
                             r'L1 $p_T^\mathrm{BDT} > 10 \;\mathrm{GeV}$',
                             r'L1 $p_T^\mathrm{NN} > 10 \;\mathrm{GeV}$'], fontsize='xx-large')
    ax.set_xlabel(xlabel, fontsize='xx-large')
    ax.set_ylabel('Efficiency', fontsize='xx-large')
    ax.tick_params(axis='both', labelsize='xx-large')

    plt.minorticks_on()
    plt.grid(b=True, color='lightgrey', linestyle='--', linewidth=1, which='both')

    props = dict(boxstyle='round', facecolor='white')
    ax.text(0.8, 0.85, f'Denom: {sum(n_den)}', transform=ax.transAxes, bbox=props)

    fig.savefig(f'{folder}{title}.pdf', bbox_inches=None)

#############################################################
# Import file print from config file
config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

version = config['MuonGun']['version']
treename = config['MuonGun']['treename']
if version == '6':
    filename1 = f'samples/v{version}/' + config['MuonGun']['filename'] + "_1.root"
    filename2 = f'samples/v{version}/' + config['MuonGun']['filename'] + "_2.root"
    filename3 = f'samples/v{version}/' + config['MuonGun']['filename'] + "_3.root"
    filename4 = f'samples/v{version}/' + config['MuonGun']['filename'] + "_4.root"
    filelist = [filename1, filename2, filename3, filename4]
else:
    filename1 = f'samples/v{version}/' + config['MuonGun']['filename']
    filelist = [filename1]

folder = f'efficiencies/v{version}/'


matched_n_pt_num_bdt = 0
matched_n_pt_num_nn = 0
matched_n_pt_num_both = 0
matched_n_pt_denom = 0
n_pt_denom = 0
matched_n_pt_lo_eta_num_bdt = 0
matched_n_pt_lo_eta_num_nn = 0
matched_n_pt_lo_eta_num_both = 0
matched_n_pt_lo_eta_denom = 0
n_pt_lo_eta_denom = 0
matched_n_pt_md_eta_num_bdt = 0
matched_n_pt_md_eta_num_nn = 0
matched_n_pt_md_eta_num_both = 0
matched_n_pt_md_eta_denom = 0
n_pt_md_eta_denom = 0
matched_n_pt_hi_eta_num_bdt = 0
matched_n_pt_hi_eta_num_nn = 0
matched_n_pt_hi_eta_num_both = 0
matched_n_pt_hi_eta_denom = 0
n_pt_hi_eta_denom = 0
matched_n_pt_lo_d0_num_bdt = 0
matched_n_pt_lo_d0_num_nn = 0
matched_n_pt_lo_d0_num_both = 0
matched_n_pt_lo_d0_denom = 0
n_pt_lo_d0_denom = 0
matched_n_pt_md_d0_num_bdt = 0
matched_n_pt_md_d0_num_nn = 0
matched_n_pt_md_d0_num_both = 0
matched_n_pt_md_d0_denom = 0
n_pt_md_d0_denom = 0
matched_n_pt_hi_d0_num_bdt = 0
matched_n_pt_hi_d0_num_nn = 0
matched_n_pt_hi_d0_num_both = 0
matched_n_pt_hi_d0_denom = 0
n_pt_hi_d0_denom = 0

matched_n_invpt_num_bdt = 0
matched_n_invpt_num_nn = 0
matched_n_invpt_num_both = 0
matched_n_invpt_denom = 0
n_invpt_denom = 0

matched_n_eta_num_bdt = 0
matched_n_eta_num_nn = 0
matched_n_eta_num_both = 0
matched_n_eta_denom = 0
n_eta_denom = 0
matched_n_fine_eta_num_bdt = 0
matched_n_fine_eta_num_nn = 0
matched_n_fine_eta_num_both = 0
matched_n_fine_eta_denom = 0
n_fine_eta_denom = 0

matched_n_d0_num_bdt = 0
matched_n_d0_num_nn = 0
matched_n_d0_num_both = 0
matched_n_d0_denom = 0
n_d0_denom = 0
matched_n_d0_lo_eta_num_bdt = 0
matched_n_d0_lo_eta_num_nn = 0
matched_n_d0_lo_eta_num_both = 0
matched_n_d0_lo_eta_denom = 0
n_d0_lo_eta_denom = 0
matched_n_d0_md_eta_num_bdt = 0
matched_n_d0_md_eta_num_nn = 0
matched_n_d0_md_eta_num_both = 0
matched_n_d0_md_eta_denom = 0
n_d0_md_eta_denom = 0
matched_n_d0_hi_eta_num_bdt = 0
matched_n_d0_hi_eta_num_nn = 0
matched_n_d0_hi_eta_num_both = 0
matched_n_d0_hi_eta_denom = 0
n_d0_hi_eta_denom = 0

gen_pt_array = np.array(())
bdt_pt_array = np.array(())
nn_pt_array  = np.array(())
gen_d0_array = np.array(())
nn_d0_array  = np.array(())
matched_array  = np.array(())
eta_array  = np.array(()) 

gen_invpt_array = np.array(())
bdt_invpt_array = np.array(())
nn_invpt_array  = np.array(())

ptbins = np.linspace(0, 60, 30)
invptbins = np.linspace(-0.02, 0.02, 30)
d0bins = np.linspace(0, 100, 60)
etabins = [-2.5, -2.1, -1.6, -1.24]
for eta in etabins[::-1]:
    etabins.append(-eta)
etabins = np.array((etabins))
ic(etabins)
fine_etabins = np.linspace(-2.5,2.5,50)

for filename in filelist:
    print("Loading file")
    ic(filename)

    with uproot.open(filename) as file:
        print("Loaded file.")

        tree = file['tree']
        table = tree.arrays()

        print("Extracting branches.")

        gen_pt      = ak.flatten(table['gen_pt'])
        gen_phi     = ak.flatten(table['gen_phi'])
        gen_eta     = ak.flatten(table['gen_eta'])
        # ic(np.sum(gen_eta) > 0)
        gen_phiStar = ak.flatten(table['gen_phiStar'])
        gen_etaStar = ak.flatten(table['gen_etaStar'])
        gen_q       = ak.flatten(table['gen_q'])
        gen_d0      = ak.flatten(table['gen_d0'])
        gen_vz      = ak.flatten(table['gen_vz'])
        gen_dR      = ak.flatten(table['gen_dR'])

        matched_gen_pt      = ak.flatten(table['matched_gen_pt'])
        matched_gen_phi     = ak.flatten(table['matched_gen_phi'])
        matched_gen_eta     = ak.flatten(table['matched_gen_eta'])
        matched_gen_phiStar = ak.flatten(table['matched_gen_phiStar'])
        matched_gen_etaStar = ak.flatten(table['matched_gen_etaStar'])
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
        match_mask = (emtf_dR_mask & gen_dR_mask)

        matched_array = np.append(matched_array, ak.to_numpy(match_mask))
        gen_pt_array = np.append(gen_pt_array, ak.to_numpy(matched_gen_pt))
        bdt_pt_array = np.append(bdt_pt_array, ak.to_numpy(BDT_pt))
        nn_pt_array  = np.append(nn_pt_array, ak.to_numpy(NN_pt))
        gen_d0_array = np.append(gen_d0_array, ak.to_numpy(matched_gen_d0))
        nn_d0_array  = np.append(nn_d0_array, ak.to_numpy(NN_d0))

        temp_invpt     = matched_gen_q/matched_gen_pt
        temp_invpt_bdt = matched_gen_q/BDT_pt
        temp_invpt_nn  = matched_gen_q/NN_pt
        gen_invpt_array = np.append(gen_invpt_array, ak.to_numpy(temp_invpt))
        bdt_invpt_array = np.append(bdt_invpt_array, ak.to_numpy(temp_invpt_bdt))
        nn_invpt_array  = np.append(nn_invpt_array, ak.to_numpy(temp_invpt_nn))
        eta_array = np.append(eta_array, ak.to_numpy(matched_gen_etaStar))

        denom_mask = (abs(gen_vz) < 100) & (abs(gen_d0) < 100)
        num_mask = (abs(matched_gen_vz) < 100) & (abs(matched_gen_d0) < 100) 

        denom_eta_mask    = (abs(gen_etaStar) > 1.2) & (abs(gen_etaStar) < 2.5)
        denom_lo_eta_mask = (abs(gen_etaStar) > 1.2) & (abs(gen_etaStar) < 1.6)
        denom_md_eta_mask = (abs(gen_etaStar) > 1.6) & (abs(gen_etaStar) < 2.1)
        denom_hi_eta_mask = (abs(gen_etaStar) > 2.1) & (abs(gen_etaStar) < 2.5)
        denom_lo_d0_mask = (abs(gen_d0) < 25)
        denom_md_d0_mask = (abs(gen_d0) > 25) & (abs(gen_d0) < 50)
        denom_hi_d0_mask = (abs(gen_d0) > 50)

        num_eta_mask    = (abs(matched_gen_etaStar) > 1.2) & (abs(matched_gen_etaStar) < 2.5)
        num_lo_eta_mask = (abs(matched_gen_etaStar) > 1.2) & (abs(matched_gen_etaStar) < 1.6)
        num_md_eta_mask = (abs(matched_gen_etaStar) > 1.6) & (abs(matched_gen_etaStar) < 2.1)
        num_hi_eta_mask = (abs(matched_gen_etaStar) > 2.1) & (abs(matched_gen_etaStar) < 2.5)
        num_lo_d0_mask = (abs(matched_gen_d0) < 25)
        num_md_d0_mask = (abs(matched_gen_d0) > 25) & (abs(matched_gen_d0) < 50)
        num_hi_d0_mask = (abs(matched_gen_d0) > 50)

        denom_pt_mask = gen_pt > 15
        num_pt_mask = matched_gen_pt > 15
        
        pt_mask_bdt = BDT_pt > 10
        pt_mask_nn = NN_pt > 10
        pt_mask_both = pt_mask_bdt | pt_mask_nn

        denom_d0_mask = gen_d0 < 20


        # ##############################################################################
        ## pT efficiency
        den_pt_mask  = denom_mask & denom_eta_mask
        den_pt_lo_eta_mask  = denom_mask & denom_lo_eta_mask
        den_pt_md_eta_mask  = denom_mask & denom_md_eta_mask
        den_pt_hi_eta_mask  = denom_mask & denom_hi_eta_mask
        den_pt_lo_d0_mask  = denom_mask & denom_eta_mask & denom_lo_d0_mask
        den_pt_md_d0_mask  = denom_mask & denom_eta_mask & denom_md_d0_mask
        den_pt_hi_d0_mask  = denom_mask & denom_eta_mask & denom_hi_d0_mask

        # matching eff, 0.6, 5
        matched_num_pt_bdt_mask  = num_mask & num_eta_mask & match_mask & pt_mask_bdt
        matched_num_pt_nn_mask   = num_mask & num_eta_mask & match_mask & pt_mask_nn
        matched_num_pt_both_mask = num_mask & num_eta_mask & match_mask & pt_mask_both
        matched_num_pt_lo_eta_bdt_mask  = num_mask & num_lo_eta_mask & match_mask & pt_mask_bdt
        matched_num_pt_lo_eta_nn_mask   = num_mask & num_lo_eta_mask & match_mask & pt_mask_nn
        matched_num_pt_lo_eta_both_mask = num_mask & num_lo_eta_mask & match_mask & pt_mask_both
        matched_num_pt_md_eta_bdt_mask  = num_mask & num_md_eta_mask & match_mask & pt_mask_bdt
        matched_num_pt_md_eta_nn_mask   = num_mask & num_md_eta_mask & match_mask & pt_mask_nn
        matched_num_pt_md_eta_both_mask = num_mask & num_md_eta_mask & match_mask & pt_mask_both
        matched_num_pt_hi_eta_bdt_mask  = num_mask & num_hi_eta_mask & match_mask & pt_mask_bdt
        matched_num_pt_hi_eta_nn_mask   = num_mask & num_hi_eta_mask & match_mask & pt_mask_nn
        matched_num_pt_hi_eta_both_mask = num_mask & num_hi_eta_mask & match_mask & pt_mask_both
        matched_num_pt_lo_d0_bdt_mask  = num_mask & num_eta_mask & num_lo_d0_mask & match_mask & pt_mask_bdt
        matched_num_pt_lo_d0_nn_mask   = num_mask & num_eta_mask & num_lo_d0_mask & match_mask & pt_mask_nn
        matched_num_pt_lo_d0_both_mask = num_mask & num_eta_mask & num_lo_d0_mask & match_mask & pt_mask_both
        matched_num_pt_md_d0_bdt_mask  = num_mask & num_eta_mask & num_md_d0_mask & match_mask & pt_mask_bdt
        matched_num_pt_md_d0_nn_mask   = num_mask & num_eta_mask & num_md_d0_mask & match_mask & pt_mask_nn
        matched_num_pt_md_d0_both_mask = num_mask & num_eta_mask & num_md_d0_mask & match_mask & pt_mask_both
        matched_num_pt_hi_d0_bdt_mask  = num_mask & num_eta_mask & num_hi_d0_mask & match_mask & pt_mask_bdt
        matched_num_pt_hi_d0_nn_mask   = num_mask & num_eta_mask & num_hi_d0_mask & match_mask & pt_mask_nn
        matched_num_pt_hi_d0_both_mask = num_mask & num_eta_mask & num_hi_d0_mask & match_mask & pt_mask_both

        matched_n_pt_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_bdt_mask)
        matched_n_pt_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_nn_mask)
        matched_n_pt_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_both_mask)
        matched_n_pt_denom    += get_n(gen_pt, ptbins, den_pt_mask)
        matched_n_pt_lo_eta_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_lo_eta_bdt_mask)
        matched_n_pt_lo_eta_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_lo_eta_nn_mask)
        matched_n_pt_lo_eta_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_lo_eta_both_mask)
        matched_n_pt_lo_eta_denom    += get_n(gen_pt, ptbins, den_pt_lo_eta_mask)
        matched_n_pt_md_eta_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_md_eta_bdt_mask)
        matched_n_pt_md_eta_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_md_eta_nn_mask)
        matched_n_pt_md_eta_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_md_eta_both_mask)
        matched_n_pt_md_eta_denom    += get_n(gen_pt, ptbins, den_pt_md_eta_mask)
        matched_n_pt_hi_eta_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_hi_eta_bdt_mask)
        matched_n_pt_hi_eta_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_hi_eta_nn_mask)
        matched_n_pt_hi_eta_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_hi_eta_both_mask)
        matched_n_pt_hi_eta_denom    += get_n(gen_pt, ptbins, den_pt_hi_eta_mask)
        matched_n_pt_lo_d0_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_lo_d0_bdt_mask)
        matched_n_pt_lo_d0_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_lo_d0_nn_mask)
        matched_n_pt_lo_d0_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_lo_d0_both_mask)
        matched_n_pt_lo_d0_denom    += get_n(gen_pt, ptbins, den_pt_lo_d0_mask)
        matched_n_pt_md_d0_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_md_d0_bdt_mask)
        matched_n_pt_md_d0_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_md_d0_nn_mask)
        matched_n_pt_md_d0_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_md_d0_both_mask)
        matched_n_pt_md_d0_denom    += get_n(gen_pt, ptbins, den_pt_md_d0_mask)
        matched_n_pt_hi_d0_num_bdt  += get_n(matched_gen_pt, ptbins, matched_num_pt_hi_d0_bdt_mask)
        matched_n_pt_hi_d0_num_nn   += get_n(matched_gen_pt, ptbins, matched_num_pt_hi_d0_nn_mask)
        matched_n_pt_hi_d0_num_both += get_n(matched_gen_pt, ptbins, matched_num_pt_hi_d0_both_mask)
        matched_n_pt_hi_d0_denom    += get_n(gen_pt, ptbins, den_pt_hi_d0_mask)


        ##############################################################################
        ## 1/pT efficiency
        den_invpt_mask  = denom_mask & denom_eta_mask

        # matching eff, 0.6, 5
        matched_num_invpt_bdt_mask  = num_mask & num_eta_mask & match_mask & pt_mask_bdt
        matched_num_invpt_nn_mask   = num_mask & num_eta_mask & match_mask & pt_mask_nn
        matched_num_invpt_both_mask = num_mask & num_eta_mask & match_mask & pt_mask_both

        matched_n_invpt_num_bdt  += get_n(matched_gen_q/matched_gen_pt, invptbins, matched_num_invpt_bdt_mask)
        matched_n_invpt_num_nn   += get_n(matched_gen_q/matched_gen_pt, invptbins, matched_num_invpt_nn_mask)
        matched_n_invpt_num_both += get_n(matched_gen_q/matched_gen_pt, invptbins, matched_num_invpt_both_mask)
        matched_n_invpt_denom    += get_n(gen_q/gen_pt, invptbins, den_invpt_mask)

        ##############################################################################
        ## eta efficiency
        den_eta_mask  = denom_mask & denom_eta_mask & denom_pt_mask

        # matching eff, 0.6, 5
        matched_num_eta_bdt_mask  = num_mask & num_eta_mask & num_pt_mask & match_mask & pt_mask_bdt
        matched_num_eta_nn_mask   = num_mask & num_eta_mask & num_pt_mask & match_mask & pt_mask_nn
        matched_num_eta_both_mask = num_mask & num_eta_mask & num_pt_mask & match_mask & pt_mask_both

        matched_n_eta_num_bdt  += get_n(matched_gen_etaStar, etabins, matched_num_eta_bdt_mask)
        matched_n_eta_num_nn   += get_n(matched_gen_etaStar, etabins, matched_num_eta_nn_mask)
        matched_n_eta_num_both += get_n(matched_gen_etaStar, etabins, matched_num_eta_both_mask)
        matched_n_eta_denom    += get_n(gen_etaStar, etabins, den_eta_mask)

        matched_n_fine_eta_num_bdt  += get_n(matched_gen_etaStar, fine_etabins, matched_num_eta_bdt_mask)
        matched_n_fine_eta_num_nn   += get_n(matched_gen_etaStar, fine_etabins, matched_num_eta_nn_mask)
        matched_n_fine_eta_num_both += get_n(matched_gen_etaStar, fine_etabins, matched_num_eta_both_mask)
        matched_n_fine_eta_denom    += get_n(gen_etaStar, fine_etabins, den_eta_mask)

        ##############################################################################
        ## d0 efficiency
        den_d0_mask  = denom_mask & denom_eta_mask & denom_pt_mask
        den_d0_lo_eta_mask  = denom_mask & denom_lo_eta_mask & denom_pt_mask
        den_d0_md_eta_mask  = denom_mask & denom_md_eta_mask & denom_pt_mask
        den_d0_hi_eta_mask  = denom_mask & denom_hi_eta_mask & denom_pt_mask

        # matching eff, 0.6, 5
        matched_num_d0_bdt_mask  = num_mask & num_eta_mask & num_pt_mask & match_mask & pt_mask_bdt
        matched_num_d0_nn_mask   = num_mask & num_eta_mask & num_pt_mask & match_mask & pt_mask_nn
        matched_num_d0_both_mask = num_mask & num_eta_mask & num_pt_mask & match_mask & pt_mask_both
        matched_num_d0_lo_eta_bdt_mask  = num_mask & num_lo_eta_mask & num_pt_mask & match_mask & pt_mask_bdt
        matched_num_d0_lo_eta_nn_mask   = num_mask & num_lo_eta_mask & num_pt_mask & match_mask & pt_mask_nn
        matched_num_d0_lo_eta_both_mask = num_mask & num_lo_eta_mask & num_pt_mask & match_mask & pt_mask_both
        matched_num_d0_md_eta_bdt_mask  = num_mask & num_md_eta_mask & num_pt_mask & match_mask & pt_mask_bdt
        matched_num_d0_md_eta_nn_mask   = num_mask & num_md_eta_mask & num_pt_mask & match_mask & pt_mask_nn
        matched_num_d0_md_eta_both_mask = num_mask & num_md_eta_mask & num_pt_mask & match_mask & pt_mask_both
        matched_num_d0_hi_eta_bdt_mask  = num_mask & num_hi_eta_mask & num_pt_mask & match_mask & pt_mask_bdt
        matched_num_d0_hi_eta_nn_mask   = num_mask & num_hi_eta_mask & num_pt_mask & match_mask & pt_mask_nn
        matched_num_d0_hi_eta_both_mask = num_mask & num_hi_eta_mask & num_pt_mask & match_mask & pt_mask_both

        matched_n_d0_num_bdt  += get_n(matched_gen_d0, d0bins, matched_num_d0_bdt_mask)
        matched_n_d0_num_nn   += get_n(matched_gen_d0, d0bins, matched_num_d0_nn_mask)
        matched_n_d0_num_both += get_n(matched_gen_d0, d0bins, matched_num_d0_both_mask)
        matched_n_d0_denom    += get_n(gen_d0, d0bins, den_d0_mask)
        matched_n_d0_lo_eta_num_bdt  += get_n(matched_gen_d0, d0bins, matched_num_d0_lo_eta_bdt_mask)
        matched_n_d0_lo_eta_num_nn   += get_n(matched_gen_d0, d0bins, matched_num_d0_lo_eta_nn_mask)
        matched_n_d0_lo_eta_num_both += get_n(matched_gen_d0, d0bins, matched_num_d0_lo_eta_both_mask)
        matched_n_d0_lo_eta_denom    += get_n(gen_d0, d0bins, den_d0_lo_eta_mask)
        matched_n_d0_md_eta_num_bdt  += get_n(matched_gen_d0, d0bins, matched_num_d0_md_eta_bdt_mask)
        matched_n_d0_md_eta_num_nn   += get_n(matched_gen_d0, d0bins, matched_num_d0_md_eta_nn_mask)
        matched_n_d0_md_eta_num_both += get_n(matched_gen_d0, d0bins, matched_num_d0_md_eta_both_mask)
        matched_n_d0_md_eta_denom    += get_n(gen_d0, d0bins, den_d0_md_eta_mask)
        matched_n_d0_hi_eta_num_bdt  += get_n(matched_gen_d0, d0bins, matched_num_d0_hi_eta_bdt_mask)
        matched_n_d0_hi_eta_num_nn   += get_n(matched_gen_d0, d0bins, matched_num_d0_hi_eta_nn_mask)
        matched_n_d0_hi_eta_num_both += get_n(matched_gen_d0, d0bins, matched_num_d0_hi_eta_both_mask)
        matched_n_d0_hi_eta_denom    += get_n(gen_d0, d0bins, den_d0_hi_eta_mask)



## pt plotting

plot_eff(matched_n_pt_num_bdt, matched_n_pt_num_nn, matched_n_pt_num_both, matched_n_pt_denom, ptbins, r'pt efficiency, all eta', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')
plot_eff(matched_n_pt_lo_eta_num_bdt, matched_n_pt_lo_eta_num_nn, matched_n_pt_lo_eta_num_both, matched_n_pt_lo_eta_denom, ptbins, r'pt efficiency, low eta', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')
plot_eff(matched_n_pt_md_eta_num_bdt, matched_n_pt_md_eta_num_nn, matched_n_pt_md_eta_num_both, matched_n_pt_md_eta_denom, ptbins, r'pt efficiency, med eta', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')
plot_eff(matched_n_pt_hi_eta_num_bdt, matched_n_pt_hi_eta_num_nn, matched_n_pt_hi_eta_num_both, matched_n_pt_hi_eta_denom, ptbins, r'pt efficiency, high eta', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')
plot_eff(matched_n_pt_lo_d0_num_bdt, matched_n_pt_lo_d0_num_nn, matched_n_pt_lo_d0_num_both, matched_n_pt_lo_d0_denom, ptbins, r'pt efficiency, low d0', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')
plot_eff(matched_n_pt_md_d0_num_bdt, matched_n_pt_md_d0_num_nn, matched_n_pt_md_d0_num_both, matched_n_pt_md_d0_denom, ptbins, r'pt efficiency, med d0', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')
plot_eff(matched_n_pt_hi_d0_num_bdt, matched_n_pt_hi_d0_num_nn, matched_n_pt_hi_d0_num_both, matched_n_pt_hi_d0_denom, ptbins, r'pt efficiency, high d0', r'$p_T^\mathrm{gen} [\;\mathrm{GeV}]$')

## 1/pt plotting
plot_eff(matched_n_invpt_num_bdt, matched_n_invpt_num_nn, matched_n_invpt_num_both, matched_n_invpt_denom, invptbins, r'invpt efficiency, all eta', r'$q^\mathrm{gen}/p_T^\mathrm{gen} [\;\mathrm{GeV}]$')

## eta plotting
matched_n_eta_num_bdt[3], matched_n_eta_num_nn[3], matched_n_eta_num_both[3], matched_n_eta_denom[3] = 0, 0, 0, 0
plot_eff(matched_n_eta_num_bdt, matched_n_eta_num_nn, matched_n_eta_num_both, matched_n_eta_denom, etabins, r'eta efficiency, all eta', r'$\eta^{*\mathrm{gen}}$')
plot_eff(matched_n_fine_eta_num_bdt, matched_n_fine_eta_num_nn, matched_n_fine_eta_num_both, matched_n_fine_eta_denom, fine_etabins, r'eta efficiency, fine eta', r'$\eta^{*\mathrm{gen}}$')

## d0 plotting
plot_eff(matched_n_d0_num_bdt, matched_n_d0_num_nn, matched_n_d0_num_both, matched_n_d0_denom, d0bins, r'd0 efficiency, all eta', r'$d_0^\mathrm{gen}$')
plot_eff(matched_n_d0_lo_eta_num_bdt, matched_n_d0_lo_eta_num_nn, matched_n_d0_lo_eta_num_both, matched_n_d0_lo_eta_denom, d0bins, r'd0 efficiency, low eta', r'$d_0^\mathrm{gen} [\;\mathrm{cm}]$')
plot_eff(matched_n_d0_md_eta_num_bdt, matched_n_d0_md_eta_num_nn, matched_n_d0_md_eta_num_both, matched_n_d0_md_eta_denom, d0bins, r'd0 efficiency, med eta', r'$d_0^\mathrm{gen} [\;\mathrm{cm}]$')
plot_eff(matched_n_d0_hi_eta_num_bdt, matched_n_d0_hi_eta_num_nn, matched_n_d0_hi_eta_num_both, matched_n_d0_hi_eta_denom, d0bins, r'd0 efficiency, high eta', r'$d_0^\mathrm{gen} [\;\mathrm{cm}]$')


assert len(gen_pt_array) == len(bdt_pt_array), ic(len(gen_pt_array), len(bdt_pt_array))
assert len(gen_pt_array) == len(nn_pt_array),  ic(len(gen_pt_array), len(nn_pt_array))
assert len(gen_d0_array) == len(nn_d0_array),  ic(len(gen_d0_array), len(nn_d0_array))

npz_filename = f"dm_pt_v{version}.npz"
np.savez(npz_filename, gen_pt=gen_pt_array, bdt_pt=bdt_pt_array, nn_pt=nn_pt_array, gen_invpt=gen_invpt_array, bdt_invpt=bdt_invpt_array, nn_invpt=nn_invpt_array, gen_d0=gen_d0_array, nn_d0=nn_d0_array, match_mask=matched_array, gen_eta=eta_array)

print(f"File saved to {npz_filename}")