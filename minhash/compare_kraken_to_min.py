import csv
import sys
from statistics import mean
import time
import matplotlib.pyplot as plt
import statistics
import pandas as pd
import seaborn as sns
from collections import Counter

kraken_file = "data/kraken_standard_db_results_044fasta_ci70"
min_file = "../sra_data_NO_contam_reads.txt" #"contamination.sra_data.shortread.mycoplasma.fasta.k41.skips6.m13.fasta"

kraken_names={}
min_read=[]
pos, negs = 0, 0
#Save all of the kraken file info 
with open(kraken_file, "r") as f:
	for line1 in zip(f):
		kraken_names[''.join(line1).strip().split('\t')[1]]=''.join(line1).strip().split('\t')[2]
		if ''.join(line1).strip().split('\t')[2] == "2104":
			negs += 1
		else:
			pos += 1

tps, fps = 0, 0
#Get the title of the fasta file, parse out the id that is important 
with open(min_file, "r") as f:
	c = 0
	for line1 in zip(f):
		name = ''.join(line1).strip().split(' ')[0]
		# if float(line1.strip().split(":")[1]) < 30:
		# 	name = line1.strip().split(">")[1].split(" ")[0]
		if name in kraken_names.keys():
			if kraken_names[name] != "2104":
				tps += 1
			else:
				fps += 1

print("TPS: " + str(tps))
print("FPS: " + str(fps))
print("FNS: " + str(pos-tps))
print("TNS: " + str(negs-fps))


