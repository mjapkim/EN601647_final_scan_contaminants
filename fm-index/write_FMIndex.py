from CG_FMIndex import *
import pickle

with open('fm_indexed.pkl', 'wb') as output:
    fm = FmIndex(t)
    pickle.dump(fm, output, pickle.HIGHEST_PROTOCOL)

del fm