# Bloom filter for contamination detection

## Files included:

genome_preproc.py....preprocessing step for reference genomes to be indexed

bloomfilter.py.......class file for making bloom filters

mycoplasma_test.py...produces output for mycoplasma tests

synthetic_test.py....produces output for synthetic data

### genome_preproc.py
#### ran on reference genomes to get in proper format for bloom filter build

Takes the contiguous reference genome and splits it into all the possible kmers of length k. Where k is the length of the reads that we tested our implementation on.

For Mycoplasma: k = 101 bp
For synthetic: k = 150 bp

#### To run:

	python genome_proc.py
* to change reference genome, change the variable 'fasta' at the beginning of the main part of the file 
and the name of the output file

### bloomfilter.py
#### class for bloom filter objects, based on the geeksforgeeks implementation

The class takes in n, the number of items to be added to the filter, and p, the false positive probibility for the membership testing.

Used by importing the class into the testing file

### mycoplasma_test.py
#### runs the contamination detection using the reference mycoplasma bloom filter index

Takes the Mycoplasma reference genome and makes a bloom filter object based on it. Then, it reads in the sequencing reads of a sample and tests for membership, aka detects contamination.

#### To run: 

	python mycoplasma_test.py

* can be run with different p values and different input files for reference and sample

### synthetic_test.py
#### same as mycoplasma_test.py but for the synthetic data set

Ran with both the contamination reference genome and the target reference genome as the basis for the bloom filter.

For contaminant genome based bloom filter:

	if bloomf.check(word):

		this word is a contaminant

	else:

		this word is not in the set and probably not a contaminant

For target genome based bloom filter:
	if bloomf.check(word):

		this word is 'clean,' it is in the target genome

	else:

		this word is probably a contaminant, it is NOT in the target genome 
