#import csv
import sys
from statistics import mean
import time
import statistics
import pandas as pd

#Example command line argument 
#python ./find_contam_with_min_real_data.py sra_data.shortread.fasta 13 1 10 

#Take in a bunch of arguments. 
fasta_file = sys.argv[1]
m=int(sys.argv[2])
k=int(sys.argv[3])
it=int(sys.argv[4])
print("M: " + str(m))
#The prefix will make naming things easier. 
prefix=fasta_file + ".k" + str(k) + ".skips" + str(it) + ".m" + str(m)

#We're going to save a lot of output files. Most for diagnostic purposes. 
out_db = open("db_size_" + prefix + ".txt","a")
target_out = open("target." + prefix + ".fasta","w")
contam_out = open("contamination." + prefix + ".fasta","w")
#potential_contam_out = open("potential.contamination." + prefix + ".fasta","w")
query_time_out = open("query_time_" + prefix + ".txt","a")
db_time_out = open("db_time_" + prefix + ".txt","a")

#This isn't as relevant anymore, but I left the kmer stuff in here for reference. 
def update_kmerdb(kmer, kmerDB):
    if kmer in kmerDB.keys():
        count = kmerDB[kmer][0] + 1
    else:
        count = 1
    kmerDB[kmer] = [count]
    return kmerDB, count
#Update the minimizer database. 
#Basically, just keep of a count of every time a min value occurs. 
def update_mindb(mins, minDB):
	if mins in minDB.keys():
		count = minDB[mins][0] + 1
	else:
		count = 1
	minDB[mins] = [count]
	return minDB

#Self explainatory 
def reverse_comp(s):
	r=s[::-1]
	rev=""
	alphabet = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	for i in r:
		rev=rev+alphabet[i]
	return rev

#Just returns all of the kmers in a string (k is taken as an argument.)
def find_kmers(line):
    kmers = []
    for i in range(0, len(line)+1-k):
        kmers.append(line[i:i+k])
        #Got to get the reverse as well 
        kmers.append(line[len(line)-k-i:len(line)-i])
    return kmers

def string_to_kmers(s, length):
	#This is adapted from 
    #https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/Minimizers.ipynb
    return [s[i:i+length] for i in range(len(s)-length+1)]

def minimizer(line, m):
    #This is adapted from 
    #https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/Minimizers.ipynb
    assert m <= len(line)
    return min(string_to_kmers(line, m))

def find_min(line):
	mins = []
	for i in range(0, len(line)+1-k, it):
		mini='vvvvoyager'
		l = k-m+1
		for j in range(0, l):
			forward = line[i:i+k][j:j+m]
			reverse = line[len(line)-k-i:len(line)-i][j:j+m]
			mini=min(mini,forward,reverse)		
		mins.append(mini)
	return mins


minDB = {}
all_reads = []
all_names = []
#Load in the reads
with open(fasta_file, "r") as f:
	for line1, line2 in zip(f,f):
		all_reads.append(line2.strip())
		all_names.append(line1.strip())
print("Read in the fasta file ")
all_mins = []

#Make the database 
t0=time.time()
for ele in all_reads:
	line_mins = find_min(ele)
	all_mins.append(line_mins)
	for mins in line_mins:
		minDB = update_mindb(mins, minDB)
t1=time.time()

#Save the time to a file 
db_time_out.write(str(t1-t0) + "\t" +  str(m) + "\t" + str(it) + sys.argv[1] +"\n")
print("Made db in: " + str(t1-t0) + "\tM: " + str(m) + "\tSkips: " + str(it))
f.close()

print("Finished making DB")
#There simply must be a more efficient way to figure out the size of the database
df=pd.DataFrame(minDB)
out_db.write(str(df.memory_usage(index=True).sum()) + "\t" + str(it) + "\t" + str(k) + "\t" + str(m) + "\n")
out_db.close()


all_counts = []
all_means = []

i = 0
t0=time.time()
for read in all_reads:
	if len(read) >= k:
		read_min = all_mins[i]
		this_read = []
		for mins in read_min:
			this=minDB[mins][0]
			this_read.append(this)
		all_means.append(statistics.median(this_read)) #this is the value that we use to classify the read 
	i += 1
t1=time.time()
#Write the time it took to query 
query_time_out.write(str(t1-t0) + "\t" +  str(m) + "\t" + str(it) + sys.argv[1] +"\n")
print("Query db in: " + str(t1-t0) + "\tM: " + str(m) + "\tSkips: " + str(it))

x=[]
classification = []
threshold = 50 
#Just playing around with the idea of a different threshold. 
#Didnt actually end up using the "potential threshold". 
potential_threshold = statistics.median(all_means)/10
print("Threshold: " + str(threshold))
print("Potential threshold: " + str(potential_threshold))

#Send the reads to different output files, depending on their relation to the threshold. 
positives = 0
negatives = 0
for i in range(0,len(all_means)):
	x.append(i)
	av = all_means[i]
	if av < threshold:
		contam_out.write(all_names[i] + ":" + str(av) + "\n")
		contam_out.write(all_reads[i] + "\n")
		positives += 1
	else:
		target_out.write(all_names[i] + ":" + str(av) + "\n")
		target_out.write(all_reads[i] + "\n")
		negatives += 1

print("Finished")


query_time_out.close()
db_time_out.close()
#




