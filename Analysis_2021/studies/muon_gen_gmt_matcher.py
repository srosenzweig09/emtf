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

from uproot_open import get_uproot_Table
from kinematics import calcDeltaR, calcStar, convert_emtf
from logger import info

def match_muons(evt_deltaR, evt, threshold):
    inds = evt_deltaR.argsort(axis=1)
    sorted_evt_deltaR = np.sort(evt_deltaR, axis=1)

    # Create dictionary where keys = GEN and values = GMT
    gen2gmt_dict = {} # dict[gen_ind] = gmt_ind
    gmt_repeats = []
    # Loop through gmt_ind (n = length of inds[:,0]) and see which gen muons
    # had the lowest DeltaR value to each muon. Multiples may occur.
    for gmt_ind, gen_ind in enumerate(inds[:,0]):
        # Don't keep any pairs with DeltaR below the matching threshold
        if evt_deltaR[gmt_ind, gen_ind] > threshold:
            continue
        elif gmt_ind not in gen2gmt_dict.values():
            # If gmt muon isn't assigned to a gen muon, go ahead and assign it as a pair
            gen2gmt_dict[gen_ind] = gmt_ind
        elif gmt_ind in gen2gmt_dict.values():
            # If gmt muon is already paired to a gen muon, assign it as a repeat
            gmt_repeats += [gmt_ind]
    
    # To handle multiple gmt muons being assigned to gen muons,
    # go through the list of repeats and find the muon with the lower DeltaR
    # and assign that gmt muon as the pair to the gen muon
    # Assign the duplicate to its next best gen muon match and
    # check if the DeltaR is below the threshold.
    for gmt_ind in set(gmt_repeats):
        for gen_ind, dR in enumerate(sorted_evt_deltaR[gmt_ind, :]):
            if dR >= threshold:
                continue
            elif gmt_ind not in gen2gmt_dict.values():
                gen2gmt_dict[gen_ind] = gmt_ind

    gmt_indices = np.array((), dtype=int)
    gen_indices  = np.array((), dtype=int)

    # assert len(gen_indices) == len(gmt_indices), f"Event {evt}"

    return gen2gmt_dict

threshold = 0.6


config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

HtoLL    = config['HtoLL']['filename']
treename = config['HtoLL']['treename']

tree = get_uproot_Table(HtoLL, treename)

nevents = len(tree['genPart_pt'])
muon = (abs(tree['genPart_ID']) == 13) & (tree['genPart_parentID'] == 6000113)
muons = tree['genPart_pt'][muon].flatten()
nmuons = len(muons)

print(nmuons)

tot_gen = 0
tot_gmt = 0
tot_emtf = 0

gen_mask = []
emtf_mask = []
gen2emtf_indices = []
emtf_indices = []
flat_gen_indices = []
flat_emtf_indices =  []
flat_gmt_indices =  []

for evt in range(nevents):
    if evt%10000 == 0:
        info(f"Analyzing event {evt}/{nevents} with index {evt}")

    # assert(len(gen_index) == len(gmt_index))

    ngmt = len(tree['gmtMuon_pt'][evt])
    evt_muons = muon[evt]
    ngen = len(tree['genPart_pt'][evt][evt_muons])
    nemtf = len(tree['emtfTrack_pt'][evt])

    begin_gen = tot_gen
    begin_emtf = tot_emtf
    begin_gmt = tot_gmt

    tot_gen += ngen
    tot_gmt += ngmt
    tot_emtf += nemtf

    if ngmt == 0 or nemtf == 0:
        gen2emtf_indices.append([-1 for _ in np.arange(ngen)])
        gen_mask.append([False for _ in np.arange(ngen)])
        emtf_mask.append([False for _ in np.arange(nemtf)])
        continue

    mu_eta = tree['genPart_eta'][evt][evt_muons]
    mu_phi = tree['genPart_phi'][evt][evt_muons]
    mu_vx = tree['genPart_vx'][evt][evt_muons]
    mu_vy = tree['genPart_vy'][evt][evt_muons]
    mu_vz = tree['genPart_vz'][evt][evt_muons]

    mu_etastar, mu_phistar = calcStar(mu_eta, mu_phi, mu_vx, mu_vy, mu_vz, is2darray=False)

    gmt_evt_pt, gmt_evt_eta, gmt_evt_phi = tree['gmtMuon_pt'][evt],  tree['gmtMuon_eta'][evt], tree['gmtMuon_phi'][evt]
    emtf_evt_pt, emtf_evt_eta, emtf_evt_phi = convert_emtf(tree['emtfTrack_pt'][evt],  tree['emtfTrack_eta'][evt], tree['emtfTrack_phi'][evt], is2darray=False)

    gmt_gen_deltaR = np.array(())
    for i in range(ngmt):
        gmt_eta = np.repeat(gmt_evt_eta[i], ngen)
        gmt_phi = np.repeat(gmt_evt_phi[i], ngen)
        deltaR = calcDeltaR(gmt_eta, mu_etastar, gmt_phi, mu_phistar)
        gmt_gen_deltaR = np.append(gmt_gen_deltaR, deltaR)

    emtf_gmt_deltaR = np.array(())
    for i in range(nemtf):
        emtf_eta = np.repeat(emtf_evt_eta[i], ngmt)
        emtf_phi = np.repeat(emtf_evt_phi[i], ngmt)
        deltaR = calcDeltaR(emtf_eta, gmt_evt_eta, emtf_phi, gmt_evt_phi)
        emtf_gmt_deltaR = np.append(emtf_gmt_deltaR, deltaR)

    gmt_gen_deltaR = gmt_gen_deltaR.reshape(ngmt, ngen)
    gen_gmt_dict = match_muons(gmt_gen_deltaR, evt,  threshold) #  matching occurs here
    gen_gmt_indices = gen_gmt_dict.keys()
    gmt_gen_indices = gen_gmt_dict.values()

    emtf_gmt_deltaR = emtf_gmt_deltaR.reshape(nemtf, ngmt)
    gmt_emtf_dict = match_muons(emtf_gmt_deltaR, evt,  5) #  matching occurs here
    gmt_emtf_indices = gmt_emtf_dict.keys()
    emtf_gmt_indices = gmt_emtf_dict.values()    
    
    if len(gen_gmt_indices) < 1:
        gen2emtf_indices.append([-1 for _ in np.arange(ngen)])
        gen_mask.append([False for _ in np.arange(ngen)])
        emtf_mask.append([False for _ in np.arange(nemtf)])
        continue
    else:
        
        evt_gen_mask = []
        emtf_index = []
        gen2emtf_index = []
        for gen_mu in np.arange(ngen):
            try:
                gmt_match = gen_gmt_dict[gen_mu]
                try:
                    emtf_match = gmt_emtf_dict[gmt_match]
                    gen2emtf_index.append(gen_mu)
                    emtf_index.append(emtf_match)
                    flat_emtf_indices.append(emtf_match + begin_emtf)
                    flat_gmt_indices.append(gmt_match + begin_gmt)
                    evt_gen_mask.append(True)
                except:
                    gen2emtf_index.append(-1)
                    emtf_index.append(-1)
                    evt_gen_mask.append(False)
            except:
                gen2emtf_index.append(-1)
                emtf_index.append(-1)
                evt_gen_mask.append(False)
            
        evt_emtf_mask = [True if i in emtf_index else False for i in np.arange(nemtf)]
        gen2emtf_indices.append(gen2emtf_index)
        gen_mask.append(evt_gen_mask)
        emtf_mask.append(evt_emtf_mask)
        
        flat_gen_indices.append([i + begin_gen for i in gen2emtf_index if i!=-1])
        emtf_indices.append(emtf_index)
        

        # if evt > 10: break

info("Event analysis completed.")

flat_emtf_indices = np.array((flat_emtf_indices))
flat_gmt_indices = np.array((flat_gmt_indices))
flat_gen_indices = np.array(([item for sublist in flat_gen_indices for item in sublist]))

# Awakward persist refuses to overwrite so we must remove them before we can use them.
if os.path.isfile("matching/awk_gen_mask.awkd"):
    os.remove("matching/awk_gen_mask.awkd")
if os.path.isfile("matching/awk_gen_index.awkd"):
    os.remove("matching/awk_gen_index.awkd")
if os.path.isfile("matching/awk_emtf_mask.awkd"):
    os.remove("matching/awk_emtf_mask.awkd")
if os.path.isfile("matching/awk_emtf_index.awkd"):
    os.remove("matching/awk_emtf_index.awkd")    

gen_mask = awk.fromiter(gen_mask)
awk.persist.save("matching/awk_gen_mask", gen_mask)
awk_gen_index = awk.fromiter(gen2emtf_indices)
awk.persist.save("matching/awk_gen_index", awk_gen_index)

emtf_mask = awk.fromiter(emtf_mask)
awk.persist.save("matching/awk_emtf_mask", emtf_mask)
awk_emtf_index = awk.fromiter(emtf_indices)
awk.persist.save("matching/awk_emtf_index", awk_emtf_index)

nemtf = len(tree['emtfTrack_pt'].flatten())

assert nmuons == tot_gen, f"nmuons = {nmuons}, tot_gen = {tot_gen}"
assert nemtf == tot_emtf, f"nemtf = {nemtf}, tot_emtf = {tot_emtf}"

before_decimal = int(np.floor(threshold))
after_decimal = int((threshold - before_decimal)*10)

assert(len(flat_emtf_indices) == len(flat_gen_indices))

filename = f'matching/gen_gmt/gen_gmt_matching_masks_{before_decimal}pt{after_decimal}_pre8.npz'
np.savez(filename, emtf_indices=flat_emtf_indices, gen_indices=flat_gen_indices, gmt_indices=flat_gmt_indices)
print(filename)
