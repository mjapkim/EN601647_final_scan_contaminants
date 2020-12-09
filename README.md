## EN.601.647 Contamination Scanning/Detection via Genome Indexing

### Completed during the Fall 2020 semester by Natalia Rincon, April Kim, Katie Jenike, Jessica Bonnie

Identification of contamination of sequencing reads by foreign biological samples remains an active challenge within the field of computational genomics. For the purposes of this project, we selected 3 different indexing methods to examine their efficacy on this problem: Kmer Profiling, FM-Indexing, and Bloom Filter Indexing. Each approach was evaluated for performance on both real world and synthetic data and adapted from its original implementation through parameter selection or space reduction techniques to improve its accuracy and/or efficiency. Our analyses have shown that, while it runs in a reasonable amount of time, kmer profiling produces conservative results; current implementations of Bloom filter indexing cannot be used as a standalone method for contamination detection; and the resources required to utilize unmodified FM-Indexing on genomic data are enormous. Future work on these and other methods will be required to conquer this challenge.

Python implementations unique to each method are saved in the respective folders (bloom_filter, minimizers, minhash, and fm-index), with README files unique to each method that contains file explanations and how to run the codes. 

Prior to implementing/analyzing either of the four methods, please run the install_python_packages.sh script (./install_python_packages.sh) to ensure all packages necessary are installed in your environment.
