import csv
import sys
import random
import statistics
import numpy as np
import scipy.stats as stats
from math import floor
from collections import defaultdict

def random_string_exclude(char):
	#When we want to pick a random letter (ACTG) but it can't be the same character as provided
	new_char = char
	while new_char == char:
		new_char = random.choice("ACTG")
	return new_char

def weighted_string(length):
	string = ""
	for i in range(length):
		nuc = random.choice("ACTG")
		string += nuc
	return string

def add_het(first_genome, het_level):
	nums = random.sample(range(0,len(first_genome)), int((het_level/100)*len(first_genome)))
	nums.sort()
	new_genome = list(first_genome)
	for i in nums:
		new_genome[i] = random_string_exclude(new_genome[i]) 
	return "".join(new_genome)

def add_error(string, error_rate):
	read_w_errors=list(string)
	nums = random.sample(range(0,len(string)), int((error_rate/100)*len(string)))
	nums.sort()
	for i in nums:
		read_w_errors[i] = random_string_exclude(read_w_errors[i])
	return "".join(read_w_errors)



first_gs = 100000 #size of the target genome 
second_gs = 100000 #Size of the contamination genome
perc_cont = 0 #Percent of contamination
first_perc_het = 0.5 #Percent of het in first genome
second_perc_het = 0.5 #Percent of het in the contamintion genome
cov = 30 #Coverage 
read_size = 10000 #Size of the reads 
perc_err = 1 #Percent of error in the sequening 


out=open("simulation.gs"+str(first_gs)+".cov"+str(cov)+".het"+str(first_perc_het)+".err"+str(perc_err)+".contamination"+str(perc_cont)+".fasta", "w")
out_mat_target=open("target_mat_genome"+".fasta", "w")
out_pat_target=open("target_pat_genome"+".fasta", "w")
out_mat_contam=open("contamination_mat_genome"+".fasta", "w")
out_pat_contam=open("contamination_pat_genome"+".fasta", "w")
#Make the first and contamination genomes 

m_genome_target = weighted_string(first_gs)
p_genome_target = add_het(m_genome_target, first_perc_het)
out_mat_target.write( ">Mat Target" + "\n")
out_mat_target.write( m_genome_target+ "\n")
out_pat_target.write( ">Pat Target" + "\n")
out_pat_target.write( p_genome_target+ "\n")

m_genome_contam = weighted_string(second_gs)
p_genome_contam = add_het(m_genome_contam, second_perc_het)
out_mat_contam.write( ">Mat contamination" + "\n")
out_mat_contam.write( m_genome_contam+ "\n")
out_pat_contam.write( ">Pat contamination" + "\n")
out_pat_contam.write( p_genome_contam+ "\n")

#Now we have our genomes. Next we need to make the reads 
num_reads = int(( first_gs * cov ) / read_size)
num_reads_target = int(num_reads * ((100-perc_cont)/100) )
num_reads_contam = int(num_reads * (perc_cont/100))

for i in range(0, num_reads_target):
	#Find a random spot in the genome, and make a read 
	random_number = random.randint(0,first_gs)
	if (i%2) == 0: #Use mat 
		this_read = m_genome_target[random_number:(random_number + read_size)]
		which_genome="mat_target"
	else: # Use pat 
		this_read = p_genome_target[random_number:(random_number + read_size)]
		which_genome = "pat_target"
	this_read = add_error(this_read, perc_err)
	out.write(">read_"+str(i)+"_"+which_genome + "\n")
	out.write(this_read + "\n")

for i in range(0, num_reads_contam):
	#Find a random spot in the genome and make a read 
	random_number = random.randint(0,second_gs)
	if (i%2) == 0: #Use mat 
		this_read = m_genome_contam[random_number:(random_number + read_size)]
		which_genome="mat_contamination"
	else: # Use pat 
		this_read = p_genome_contam[random_number:(random_number + read_size)]
		which_genome = "pat_contamination"
	this_read = add_error(this_read, perc_err)
	out.write(">read_"+str(i)+"_"+which_genome + "\n")
	out.write(this_read + "\n")

out.close()
out_mat_target.close()
out_pat_target.close()
out_mat_contam.close()
out_pat_contam.close()




