from configparser import ConfigParser
import awkward as ak
from icecream import ic
import numpy as np
np.seterr(all='ignore')
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from myuproot import open_up
import sys

N_bx = 2736 # IP5, CMS
# https://indico.cern.ch/event/751857/contributions/3259413/attachments/1781638/3257666/ExperimentsInRun3_proceedings.pdf
f_LHC = 

filename = '/eos/cms/store/user/eyigitba/emtf/L1Ntuples/Run3/L1TMenuStudies/Run3Rates/NuGun_11.05.21/nuGun_11May21.root'

tree, ak_table, np_table = open_up(filename, 'tree')

