# Libraries
import sys
import gzip # Provides support for *.gz (compressed) files
import time # Used for runtime measurement
time_start = time.time() # start timer

print("\nProgram: [R]apid [R]esistance [F]inder")
print("Tool for identifying which resistance genes are in the sample.")
print("Version: 1.0")
print("\nUsage: RRF.py <FASTA_resistance_genes> <FASTQ_read_file1.gz> <FASTQ_read_file2.gz> <length_of_kmer>\n")

# Check user input
if len(sys.argv) == 1: # no command line argument
    print("No command line arguments.")
    try:
        ResFile = open(str(input("Path to FASTA resistance gene file: ")), "r")
    except IOError as error:
        print("Can't open file, reason:", str(error))
        sys.exit(1)
    Reads1_path = str(input("Path to FASTQ (*.gz) read 1 file: ")) # Unknown3_raw_reads_1.txt.gz
    Reads2_path = str(input("Path to FASTQ (*.gz) read 2 file: ")) # Unknown3_raw_reads_2.txt.gz
    try:
        Reads1File = gzip.open(Reads1_path, "rt")
        Reads2File = gzip.open(Reads2_path, "rt")
    except IOError as error:
        print("Can't open file, reason:", str(error))
        sys.exit(1)
    try:
        kmer_length = int(input("kmer length as integer (recommended is 19): "))
    except ValueError as error:
        sys.stderr.write("Kmer is not an integer. Exiting.")
        sys.exit(1)
elif len(sys.argv) == 5: # user provided all command line arguments
    try:
        ResFile = open(sys.argv[1], "r")
        Reads1File = gzip.open(sys.argv[2], "rt")
        Reads2File = gzip.open(sys.argv[3], "rt")
    except IOError as error:
        print("Can't open file, reason:", str(error))
        sys.exit(1)
    try:
        kmer_length = int(sys.argv[4])
    except ValueError as error:
        sys.stderr.write("Kmer is not an integer. Exiting.")
        sys.exit(1)
else:
    sys.stderr.write("Usage: RRF.py <FASTA_resistance_genes> <FASTQ_read_file1.gz> <FASTQ_read_file2.gz> <length_of_kmer>\n")
    sys.exit(1)

# Functions
def process_res_gene_file(ResFile):
    """Iterate the resistance gene file.
    Store gene header and its sequence in a list
    Generate a kmer set and a dictionary consisting of kmers and value of 0"""
    fasta = [] # Will contain resistance gene sequence temporarily

    for line in ResFile:
        line = line.rstrip()
        if line.startswith(">"):
            if fasta:
                full_seq = "".join(fasta) # Store entire resistance gene sequence
                store_fasta_seqs.append(header)
                store_fasta_seqs.append(full_seq)

                # Generate kmers of the resistance gene,
                # add them to a kmer set, and
                # add them to a kmer dictionary with value 0
                for i in range(0, len(full_seq) - kmer_length + 1, 1):
                    kmer = full_seq[i:i+kmer_length]
                    ResKmerSet.add(kmer)
                    if kmer not in ResKmerDict.keys():
                        ResKmerDict[kmer] = 0
                fasta = []
            header = line
        else:
            sequence = line
            fasta.append(sequence)

    # Process last resistance gene
    if fasta:
        full_seq = "".join(fasta) # Store entire resistance gene sequence
        store_fasta_seqs.append(header)
        store_fasta_seqs.append(full_seq)

        # Generate kmers of the resistance gene,
        # add them to a kmer set, and
        # add them to a kmer dictionary with value 0
        for i in range(0, len(full_seq) - kmer_length + 1, 1):
            kmer = full_seq[i:i+kmer_length]
            ResKmerSet.add(kmer)
            if kmer not in ResKmerDict.keys():
                ResKmerDict[kmer] = 0

    ResFile.close()

def reverse_complement(DNA_sequence):
    """Reverse complement a given DNA sequence"""
    DNA_translation_table = str.maketrans("ACGT", "TGCA")
    rev_compliment = DNA_sequence.translate(DNA_translation_table)
    return rev_compliment[::-1]

def read_in_seq_reads(read, kmer_length):
    """Go through sequencing reads and check if at least 2 out of 3 kmer match
    the resistance gene kmer set. If so, generate all kmers for that read and
    increase the count for matching kmers. Also consider the reverse complement."""

    # Variables
    read_seq = ""
    rev_read_seq = ""
    read_kmer = ""
    line_count = 3

    for read_line in read:
        if (line_count % 4 == 0):
            read_seq = read_line.rstrip() # arrived at read sequence

            # Generate 3 read kmers and check if at least 2 match to the
            # resistance gene dictionary. If so, generate all kmers for
            # that read, otherwise try the reverse complement. If the
            # reverse complement does not have at least 2 matches, ignore
            # this read.
            read_kmer_set = set()
            read_kmer_set.add(read_seq[1:1+kmer_length])
            read_kmer_set.add(read_seq[41:41+kmer_length])
            read_kmer_set.add(read_seq[81:81+kmer_length])

            if len(read_kmer_set.intersection(ResKmerSet)) > 1:
                for j in range(0, len(read_seq)  - kmer_length + 1, 1):
                    read_kmer = read_seq[j:j+kmer_length]
                    if read_kmer in ResKmerDict.keys():
                        ResKmerDict[read_kmer] += 1
            else:
                rev_read_seq = reverse_complement(read_seq) # reverse complement

                read_kmer_set = set()
                read_kmer_set.add(rev_read_seq[1:1+kmer_length])
                read_kmer_set.add(rev_read_seq[41:41+kmer_length])
                read_kmer_set.add(rev_read_seq[81:81+kmer_length])

                if len(read_kmer_set.intersection(ResKmerSet)) > 1:
                    for j in range(0, len(rev_read_seq) - kmer_length + 1, 1):
                        read_kmer = rev_read_seq[j:j+kmer_length]
                        if read_kmer in ResKmerDict.keys():
                            ResKmerDict[read_kmer] += 1
        line_count += 1

def filter_res_genes(store_fasta_seqs):
    counter = 0
    sequence = ""

    ### plots
    import matplotlib.pyplot as plt
    import numpy as np
    plt_counter = 1
    ### end plots

    for i in range(len(store_fasta_seqs)):
        if (counter % 2 == 0):
            # arrived at gene id
            skip_sequence = True
            header = store_fasta_seqs[i]
        else:
            # arrived at gene sequence
            temp = list()
            sequence = store_fasta_seqs[i]
            len_of_sequence = len(sequence)

            # Generate kmers of the sequence.
            # Add the depth to a kmer depth list, and check for coverage.
            # As soon as the coverage drops below 95%, ignore the sequence.
            # (= does not fulfill the requirements).
            for j in range(0, len(sequence) - kmer_length + 1, 1):
                read_kmer = sequence[j:j+kmer_length]
                temp_depth = ResKmerDict[read_kmer] >= 10

                if (temp_depth):
                    depth_of_kmer = ResKmerDict[read_kmer]
                    temp.append(depth_of_kmer)
                    #temp_coverage = 1 - temp.count(0) / len_of_sequence
                    print(header)
                    print(len(temp) / (len(sequence) - kmer_length + 1))
            gene_coverage = len(temp) / (len(sequence) - kmer_length + 1)
            if (gene_coverage >= 0.95): # check coverage, stop once below 95%
                # start plots
                #x = np.linspace(0, len(temp), len(temp))
                #plt.plot(x, temp, "-b")
                #plt.ylim(0, max(temp))
                #plt.title(header)
                #plt.ylabel("depth")
                #plt.xlabel("gene position")
                #plt.savefig(fname="depth_coverage_" + str(plt_counter) + ".png", format="png", dpi=300)
                #plt.clf()
                #plt_counter += 1
                # end plots

                # Once we processed every kmer of the sequence and the coverage
                # is still >95%, we check for minimum depth (>=10)
                final_result.append([header, gene_coverage, min(temp)])
        counter += 1

    # Sort results according to coverage and then minimum depth. Print results.
    final_result.sort(key=lambda x: (x[1], x[2]), reverse=True)

print("Processing resistance genes...")
store_fasta_seqs = [] # List of header and sequence of resistance file
ResKmerSet  = set()   # Resistance genes kmer set
ResKmerDict = dict()  # Unique resistance genes kmer dictionary {kmer1:0, kmer2:0, ...}
process_res_gene_file(ResFile)
print("Finished processing resistance gene file.")

print("Processing read file 1...")
read_in_seq_reads(Reads1File, kmer_length)
Reads1File.close()

print("Processing read file 2...")
read_in_seq_reads(Reads2File, kmer_length)
Reads2File.close()

print("Filtering resistance genes that do not match the requirements (>=95% coverage, >=10 depth)...")
final_result = list()
filter_res_genes(store_fasta_seqs)

# Print resistance genes that passed coverage and depth criteria
print("Coverage [%]\tMinimum depth\tGene")
for i in range(len(final_result)):
    print("{:.2f}\t\t{:d}\t\t{}".format(final_result[i][1] * 100, final_result[i][2], final_result[i][0]))

time_end = time.time()
time_elapsed = (time_end - time_start) / 60
print("Elapsed time in minutes: {:.2f}".format(time_elapsed))
