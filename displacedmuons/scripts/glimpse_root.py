from configparser import ConfigParser
import numpy as np
import awkward0 as awk

from myuproot import open_up
from kinematics import calcDeltaR, calcStar, convert_emtf
from logger import info

config_file = 'config/input_files.cfg'
config = ConfigParser()
config.optionxform = str
config.read(config_file)

HtoLL    = config['HtoLL']['filename']
treename = config['HtoLL']['treename']

tree = open_up(HtoLL, treename)

nevents = len(tree['genPart_pt'])
muon = (abs(tree['genPart_ID']) == 13) & (tree['genPart_parentID'] == 6000113)
nmuons = len(tree['genPart_pt'][muon].flatten())

# print(tree.columns)

# with open("branches.txt", "w") as f:
#     for branches in tree.columns:
#         f.writelines(branches + '\n')



efe_file = 'matchedNtuple_HTo2LLTo4Mu_combined_cmssw_11_0_2_fwImplementation_NNv5.root'
efe_tree = open_up(efe_file, 'tree')

print(efe_tree.columns)

for evt in np.arange(10):
    print(efe_tree['l1_pt'][evt])


masks = np.load("matching/gen_gmt/gen_gmt_matching_masks_0pt6.npz")


