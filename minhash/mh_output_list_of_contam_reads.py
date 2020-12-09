import pickle
import time
from operator import itemgetter 

def getList(dict): 
      
    return list(map(itemgetter(0), dict.items())) 

input1 = open('pkl/simulated_minhash_classify_contam.pkl', 'rb')
res = pickle.load(input1)

clean = {}
contam = {}
for key, jaccard_sim in res.items():
    if jaccard_sim != 1:
        contam[key] = jaccard_sim
    else:
        clean[key] = jaccard_sim

# print("# of contam reads: ", len(contam))
# print("# of clean reads: ", len(clean))

print(getList(contam)) 