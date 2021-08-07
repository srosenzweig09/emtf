# from configparser import ConfigParser
# import awkward as ak
# from icecream import ic
import numpy as np
np.seterr(all='ignore')
# import matplotlib as mpl
# from matplotlib.lines import Line2D
# import matplotlib.pyplot as plt
from myuproot import open_up
# import sys

N_bx = 2736 # IP5, CMS (Run 3)
# https://indico.cern.ch/event/751857/contributions/3259413/attachments/1781638/3257666/ExperimentsInRun3_proceedings.pdf
f_LHC = 11.246; # khZ

filename = '/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/crabOut/Nu_E10-pythia8-gun/NuGun_11_3_0_pre5_NNv6_5M/210421_154733/NuGun_11_3_0_pre5_NNv6_5M.root'

tree, ak_table, np_table = open_up(filename, 'EMTFNtuple/tree')

