"""
This script matches GEN muons to GMT muons.  The output is a gen_mask (applied to gen muons to find those matched to GMT muons) and a gmt_mask (applied to gmt muons to find those matched to GEN muons).

Awkward arrays are used to mask matching muons without worrying about order.

Numpy masks are used to see gen, gmt relationships. (The masks will make sure the GEN and GMT  muons are in the same order when plotting.)
"""

from configparser import ConfigParser
import numpy as np
import awkward0 as awk
import os
import itertools

from myuproot import open_up
from kinematics import calcDeltaR, calcStar, convert_emtf
from logger import info, error
from colors import CYAN,W

def match_muons(evt_deltaR, threshold, evt):
    
    rows2cols_dict = {}

    if evt_deltaR.ndim > 1:
    #     # Sort dR by columns, save by index and value
        sorted_rows = np.argsort(evt_deltaR, axis=1) # sorted by col index
        sorted_cols = np.argsort(evt_deltaR, axis=0) # sorted by row index
        sorted_dR = np.sort(evt_deltaR, axis=1)

        nrows = evt_deltaR.shape[0]
        ncols = evt_deltaR.shape[1]

        # print(evt_deltaR)

        unq, count = np.unique(evt_deltaR, axis=0, return_counts=True)
        repeated_groups = unq[count > 1]

        repeats = np.array(())
        for repeated_group in repeated_groups:
            repeated_idx = np.argwhere(np.all(evt_deltaR == repeated_group, axis=1))
            repeats =  np.append(repeats, repeated_idx[1:])
        
        singles = np.array(([i for i in range(nrows) if i not in repeats]))

        if nrows <= ncols:
            col_ind = np.repeat(0, nrows)
            col_mask = sorted_rows[:,0]
            while len(set(col_mask)) <  len(col_mask[singles]):
                for i,c in enumerate(col_mask): # i increments by 1, c is the matched column index
                    if i == len(col_mask): continue
                    for j in range(i + 1,nrows): # j = i + 1
                        if j in repeats: continue
                        if c == col_mask[j]:
                            if sorted_dR[i,col_ind[i]] < sorted_dR[j,col_ind[j]]:
                                col_ind[j] += 1
                                col_mask[j] = sorted_rows[j,col_ind[j]]
                            else:
                                col_ind[i] += 1
                                col_mask[i] = sorted_rows[i,col_ind[i]]
                                c = col_mask[i]
            
            columns = sorted_cols[0,:][col_mask]
            rows = sorted_rows[:,0]
        else:
            row_ind = np.repeat(0, ncols)
            row_mask = sorted_cols[0,:]
            # print(row_mask)
            while len(set(row_mask)) <  len(row_mask):
                for i,r in enumerate(row_mask): # i increments by 1, r is the matched row index
                    if i == len(row_mask): continue
                    for j in range(i + 1,ncols): # j = i + 1
                        if r == row_mask[j]:
                            if sorted_dR[row_ind[i],i] < sorted_dR[row_ind[j],j]:
                                row_ind[j] += 1
                                row_mask[j] = sorted_cols[row_ind[j],j]
                            else:
                                row_ind[i] += 1
                                row_mask[i] = sorted_cols[row_ind[i],i]
                                r = row_mask[i]
                        # print(row_mask)

            columns = sorted_cols[0,:]
            rows = sorted_rows[:,0][row_mask]

        for row,col in zip(rows, columns):
            if evt_deltaR[col, row] < threshold and row not in repeats:
                rows2cols_dict[row] = col # Index of smallest dR for each gmt to its closest gen

    else:
        sorted_cols = evt_deltaR.argsort()
        sorted_dR = np.sort(evt_deltaR)
        ind = sorted_cols[0]
        if sorted_dR[0] < threshold:
                rows2cols_dict[ind] = 0

    return rows2cols_dict

threshold = 0.6

config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

sample = 'HtoLL'
sample = 'MuonGun'

filename    = config[sample]['filename']
treename = config[sample]['treename']

info(f"Importing ROOT file from {CYAN+filename+W}")
tree, table, nptab = open_up(filename, treename)

nevents = len(nptab['genPart_pt'])
muon = (abs(nptab['genPart_ID']) == 13) & (nptab['genPart_parentID'] == 6000113)
muons = nptab['genPart_pt'][muon].flatten()
nmuons = len(muons)

tot_gen = 0
tot_gmt = 0
tot_emtf = 0

flat_gen_indices = []
flat_emtf_indices =  []
gen_indices = []
emtf_indices =  []

ngmtmatch = 0
nemtfmatch = 0

info(f"Beginning matching process with gen-gmt dR < {threshold} and gmt-emtf dR < 5")

for evt in range(nevents):
    if evt%10000 == 0:
        info(f"Analyzing event {evt}/{nevents}")

    # if evt > 1: break
    ngmt = len(nptab['gmtMuon_pt'][evt])
    evt_muons = muon[evt]
    ngen = len(nptab['genPart_pt'][evt][evt_muons])
    nemtf = len(nptab['emtfTrack_pt'][evt])

    begin_gen = tot_gen
    begin_emtf = tot_emtf
    begin_gmt = tot_gmt

    tot_gen += ngen
    tot_gmt += ngmt
    tot_emtf += nemtf

    if ngmt == 0 or nemtf == 0:
        gen_indices.append([])
        emtf_indices.append([])
        continue

    # gen muons
    mu_eta = nptab['genPart_eta'][evt][evt_muons]
    mu_phi = nptab['genPart_phi'][evt][evt_muons]
    mu_vx = nptab['genPart_vx'][evt][evt_muons]
    mu_vy = nptab['genPart_vy'][evt][evt_muons]
    mu_vz = nptab['genPart_vz'][evt][evt_muons]
    mu_etastar, mu_phistar = calcStar(mu_eta, mu_phi, mu_vx, mu_vy, mu_vz)

    # gmt muons
    gmt_evt_eta = nptab['gmtMuon_eta'][evt]
    gmt_evt_phi = nptab['gmtMuon_phi'][evt]

    emtf_evt_eta = nptab['emtfTrack_eta'][evt]
    globPhi = (nptab['emtfTrack_sector'][evt] - 1) * 96 + nptab['emtfTrack_GMT_phi'][evt]
    globPhi = (globPhi + 600) % 576
    emtf_evt_phi = globPhi * 0.010908

    gmt_gen_deltaR = np.array(())
    for i in range(ngmt):
        # Loop over gmt muons
        # Calculate dR between gmt muon and all gen muons
        gmt_eta = np.repeat(gmt_evt_eta[i], ngen)
        gmt_phi = np.repeat(gmt_evt_phi[i], ngen)
        deltaR = calcDeltaR(gmt_eta, mu_etastar, gmt_phi, mu_phistar)
        # Append dR between gmt muon and all gen muons
        gmt_gen_deltaR = np.append(gmt_gen_deltaR, deltaR)
    # Reshape dR into ngmt rows, ngen cols
    if ngmt > 1:
        gmt_gen_deltaR = gmt_gen_deltaR.reshape(ngmt, ngen)

    emtf_gmt_deltaR = np.array(())
    for i in range(nemtf):
        emtf_eta = np.repeat(emtf_evt_eta[i], ngmt)
        emtf_phi = np.repeat(emtf_evt_phi[i], ngmt)
        deltaR = calcDeltaR(emtf_eta, gmt_evt_eta, emtf_phi, gmt_evt_phi)
        emtf_gmt_deltaR = np.append(emtf_gmt_deltaR, deltaR)
    if nemtf > 1:
        emtf_gmt_deltaR = emtf_gmt_deltaR.reshape(nemtf, ngmt)

    ###  MATCHING
    gen_gmt_dict = match_muons(gmt_gen_deltaR, threshold, evt) #  matching occurs here

    ###  MATCHING
    gmt_emtf_dict = match_muons(emtf_gmt_deltaR, 5, evt) #  matching occurs here

    # if evt < 31:
    #     print("Event",evt)
    #     print(gen_gmt_dict)


    if not bool(gen_gmt_dict) or not bool(gmt_emtf_dict):
        gen_indices.append([])
        emtf_indices.append([])
        continue
    else:
        gen2emtf_index = []
        gen_ind = []
        emtf_ind =  []
        for gen_mu in np.arange(ngen):
            try:
                gmt_match = gen_gmt_dict[gen_mu]
                ngmtmatch += 1
                try:
                    emtf_match = gmt_emtf_dict[gmt_match]
                    nemtfmatch += 1
                    flat_gen_indices.append(gen_mu + begin_gen)
                    flat_emtf_indices.append(emtf_match + begin_emtf)
                    gen_ind.append(gen_mu)
                    emtf_ind.append(emtf_match)
                except:
                    gen2emtf_index.append(-1)
            except:
                gen2emtf_index.append(-1)
        gen_indices.append(gen_ind)
        emtf_indices.append(emtf_ind)
    # if evt > 10: break

gen_indices = awk.fromiter(gen_indices)
emtf_indices = awk.fromiter(emtf_indices)

info("Event analysis completed.")

print("Num gmt  matched:",ngmtmatch)
print("Num emtf matched:",nemtfmatch)

try:
    awk.persist.save(f"matching/{sample}_gen_ind.awkd", gen_indices)
    awk.persist.save(f"matching/{sample}_emtf_ind.awkd", emtf_indices)
except:
    error(f"JaggedArray files must be deleted before they can be written! - gen_ind.awkd - emtf_ind.awkd")

flat_emtf_indices = np.array((flat_emtf_indices), dtype=int)
flat_gen_indices = np.array((flat_gen_indices), dtype=int)

print(flat_emtf_indices)
print(flat_gen_indices)

nemtf = len(nptab['emtfTrack_pt'].flatten())

assert nmuons == tot_gen, f"nmuons = {nmuons}, tot_gen = {tot_gen}"
assert nemtf == tot_emtf, f"nemtf = {nemtf}, tot_emtf = {tot_emtf}"

before_decimal = int(np.floor(threshold))
after_decimal = int((threshold - before_decimal)*10)

assert(len(flat_emtf_indices) == len(flat_gen_indices))

filename = f'matching/gen_gmt/{sample}_gen_gmt_matching_masks_{before_decimal}pt{after_decimal}.npz'
np.savez(filename, emtf_indices=flat_emtf_indices, gen_indices=flat_gen_indices)
print(filename)