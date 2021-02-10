from configparser import ConfigParser
import numpy as np
import awkward0 as awk

from uproot_open import get_uproot_Table
from kinematics import calcDeltaR, calcStar, convert_emtf
from logger import info

def muon_mask(gen_array):
    return gen_array[(abs(tree['genPart_ID']) == 13) & (tree['genPart_parentID'] == 6000113)]

def mask(gen_array):
    muons = muon_mask(gen_array)
    return muons[matching]

def match_muons(evt_deltaR, evt, threshold):

    inds = evt_deltaR.argsort(axis=1)
    sorted_evt_deltaR = np.sort(evt_deltaR, axis=1)

    gen2emtf_dict = {} # dict[gen_ind] = reco_ind
    gen_repeats = []
    emtf_repeats = []
    for emtf_ind, gen_ind in enumerate(inds[:,0]):
        if evt_deltaR[emtf_ind, gen_ind] > threshold:
            continue
        elif emtf_ind not in gen2emtf_dict.values():
            gen2emtf_dict[gen_ind] = emtf_ind
        elif emtf_ind in gen2emtf_dict.values():
            emtf_repeats += [emtf_ind]

    for i,gen_ind in enumerate(set(gen_repeats)):
        gen2emtf_dict.pop(gen_ind, None)
        dR = np.array(([val if inds[emtf_ind,0] == gen_ind else 9999 for emtf_ind,val in enumerate(sorted_evt_deltaR[:,0]) ]))
        emtf_ind_with_lowest_dR = np.argmin(dR)
        if np.min(dR) < threshold:
            gen2emtf_dict[gen_ind] = emtf_ind_with_lowest_dR
    
    for emtf_ind in set(emtf_repeats):
        for gen_ind, dR in enumerate(sorted_evt_deltaR[emtf_ind, :]):
            if dR >= threshold:
                continue
            elif emtf_ind not in gen2emtf_dict.values():
                gen2emtf_dict[gen_ind] = emtf_ind

    emtf_indices = np.array((), dtype=int)
    gen_indices  = np.array((), dtype=int)
    for key in gen2emtf_dict.keys():
        emtf_indices = np.append(emtf_indices, gen2emtf_dict[key])
        gen_indices  = np.append(gen_indices, key)

    assert len(gen_indices) == len(emtf_indices), f"Event {evt}"

    return emtf_indices, gen_indices

threshold = 1


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

gen_index = np.array(())
emtf_index = np.array(())
tot_gen = 0
tot_emtf = 0

for evt in range(nevents):
    if evt%10000 == 0:
        info(f"Analyzing event {evt+1}/{nevents} with index {evt}")

    assert(len(gen_index) == len(emtf_index))
    
    # print(f"\nEvent  {evt}\n")

    nL1 = len(tree['emtfTrack_pt'][evt])
    evt_muons = muon[evt]
    ngen = len(tree['genPart_pt'][evt][evt_muons])

    begin_gen = tot_gen
    begin_emtf = tot_emtf

    tot_gen += ngen
    tot_emtf += nL1

    if nL1 == 0:
        continue

    mu_eta = tree['genPart_eta'][evt][evt_muons]
    mu_phi = tree['genPart_phi'][evt][evt_muons]
    mu_vx = tree['genPart_vx'][evt][evt_muons]
    mu_vy = tree['genPart_vy'][evt][evt_muons]
    mu_vz = tree['genPart_vz'][evt][evt_muons]

    mu_etastar, mu_phistar = calcStar(mu_eta, mu_phi, mu_vx, mu_vy, mu_vz, is2darray=False)

    L1_evt_pt, L1_evt_eta, L1_evt_phi = convert_emtf(tree['emtfTrack_pt'][evt],  tree['emtfTrack_eta'][evt], tree['emtfTrack_phi'][evt], is2darray=False)

    evt_deltaR = np.array(())
    for i in range(nL1):
        L1_eta = np.repeat(L1_evt_eta[i], ngen)
        L1_phi = np.repeat(L1_evt_phi[i], ngen)
        deltaR = calcDeltaR(L1_eta, mu_etastar, L1_phi, mu_phistar)
        evt_deltaR = np.append(evt_deltaR, deltaR)

    evt_deltaR = evt_deltaR.reshape(nL1, ngen)

    matched = match_muons(evt_deltaR, evt,  threshold)
    

    if len(matched[0]) < 1 or len(matched[1]) < 1:
        continue
    else:
        emtf_indices = matched[0]
        gen_indices  = matched[1]

        if len(gen_index) > 0:
            for emtf, gen in zip(emtf_indices, gen_indices):
                if gen > ngen - 1:
                    print(f"!! ERROR !! {evt}")
                    print(emtf_indices, gen_indices)
                    print(ngen)
                if emtf + begin_emtf in emtf_index:
                    print(f"!! Event {evt} !! emtf")
                if gen + begin_gen in gen_index:
                    print(f"!! Event {evt} !! gen")

        emtf_indices = emtf_indices + begin_emtf
        gen_indices  = gen_indices + begin_gen

        gen_index = np.append(gen_index, gen_indices)
        emtf_index = np.append(emtf_index, emtf_indices)


info("Event analysis completed.")

assert(nmuons == tot_gen)

print(gen_index)
print(emtf_index)

print(len(set(gen_index)))
print(len(gen_index))

gen_index =  gen_index.astype(int)
emtf_index = emtf_index.astype(int)

assert(len(gen_index) == len(emtf_index))
assert(-1 not in gen_index)
assert(-1 not in emtf_index)
assert(len(set(emtf_index)) == len(emtf_index))
assert(len(set(gen_index)) == len(gen_index))


np.savez('gen_emtf_matching_masks.npz', gen_mask=gen_index, emtf_index=emtf_index, gen_muons=muon)