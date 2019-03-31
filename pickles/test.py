import pynbody                                                                              
import numpy as np
import pylab as plt 
import pynbody.analysis.profile as profile
import pynbody.filt as filt
import glob, os, pickle
import math
import pickle

galDict = pickle.load(open('agnQuenchingProperties.pkl', 'rb'))

print(galDict.keys())
