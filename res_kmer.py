# Libraries
import sys
import gzip # Provides support for *.gzip files
import time # Used for runtime measurement
time_start = time.time() # start timer

print("\nProgram: Tool for identifying which resistance genes are in the sample.")
print("Version: 1.0")
print("\nUsage: res_kmer.py <FASTA_resistance_genes> <FASTQ_read_file1.gz> <FASTQ_read_file2.gz> <length_of_kmer>\n")

if len(sys.argv) == 1: # no command line argument
    print("No command line arguments.")
    try:
        ResFile = open(str(input("Path to FASTA resistance gene file: ")), "r")
    except IOError as error:
        print("Can't open file, reason:", str(error))
        sys.exit(1)
    Reads1_path = str(input("Path to FASTQ (*.gz) read 1 file: ")) # Unknown3_raw_reads_1.txt.gz
    Reads1_path_split = Reads1_path.rsplit("_", 1)
    Reads2_path = Reads1_path_split[0] + "_" + Reads1_path_split[1].translate(str.maketrans("1", "2"))
    print("Automatically determined path to FASTQ (*.gz) read 2 file:", Reads2_path)
    try:
        Reads1File = gzip.open(Reads1_path, "rt")
        Reads2File = gzip.open(Reads2_path, "rt")
    except IOError as error:
        print("Can't open file, reason:", str(error))
        sys.exit(1)
    try:
        kmer_length = int(input("kmer length as integer (the lower the slower, recommended is 19): "))
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
    sys.stderr.write("Usage: res_kmer.py <FASTA_resistance_genes> <FASTQ_read_file1.gz> <FASTQ_read_file2.gz> <length_of_kmer>\n")
    sys.exit(1)

# Functions
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

            if len(read_kmer_set.intersection(ResKmerSet)) >= 2:
                for j in range(0, len(read_seq), 1):
                    if (j < len(read_seq) - kmer_length + 1):
                        read_kmer = read_seq[j:j+kmer_length]
                        if read_kmer in ResKmerDict.keys():
                            ResKmerDict[read_kmer] += 1
            else:
                rev_read_seq = reverse_complement(read_seq) # reverse complement

                read_kmer_set = set()
                read_kmer_set.add(rev_read_seq[1:1+kmer_length])
                read_kmer_set.add(rev_read_seq[41:41+kmer_length])
                read_kmer_set.add(rev_read_seq[81:81+kmer_length])

                if len(read_kmer_set.intersection(ResKmerSet)) >= 2:
                    for j in range(0, len(rev_read_seq), 1):
                        if (j < len(rev_read_seq) - kmer_length + 1):
                            read_kmer = rev_read_seq[j:j+kmer_length]
                            if read_kmer in ResKmerDict.keys():
                                ResKmerDict[read_kmer] += 1
        line_count += 1

print("Processing resistance genes...")
fasta = []            # Will contain resistance gene sequence temporarily
store_fasta_seqs = [] # List of header and sequence of resistance file
ResKmerSet  = set()   # Resistance genes kmer set
ResKmerDict = dict()  # Unique resistance genes kmer dictionary {kmer1:0, kmer2:0, ...}

ResFile = open("data/resistance_genes.fsa", "r")
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
            for i in range(0, len(full_seq), 1):
                kmer = full_seq[i:i+kmer_length]
                if (i < len(full_seq) - kmer_length + 1):
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
    for i in range(0, len(full_seq), 1):
        kmer = full_seq[i:i+kmer_length]
        if (i < len(full_seq) - kmer_length + 1):
            ResKmerSet.add(kmer)
        if kmer not in ResKmerDict.keys():
            ResKmerDict[kmer] = 0

ResFile.close()

print("Finished processing resistance gene file.")

print("Processing read file 1...")
read_in_seq_reads(Reads1File, kmer_length)
Reads1File.close()

print("Processing read file 2...")
read_in_seq_reads(Reads2File, kmer_length)
Reads2File.close()

print("Filtering resistance genes that do not match the requirements (>95% coverage, >=10 depth)...")

skip_sequence = True
counter = 0
sequence = ""
final_result = list()

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
        for j in range(0, len(sequence), 1):
            if (j < len(sequence) - kmer_length + 1):
                read_kmer = sequence[j:j+kmer_length]
                depth_of_kmer = ResKmerDict[read_kmer]
                temp.append(depth_of_kmer)
                temp_coverage = 1 - temp.count(0) / len_of_sequence
                if (temp_coverage < 0.95): # check coverage, stop once below 95%
                    skip_sequence = False
                    break
        # Once we processed every kmer of the sequence and the coverage
        # is still >95%, we check for minimum depth (>=10)
        if skip_sequence and (min(temp) >= 10):
            final_result.append([header, temp_coverage, min(temp)])
    counter += 1

# Sort results according to coverage and then minimum depth. Print results.
final_result.sort(key=lambda x: (x[1], x[2]), reverse=True)
print("Coverage [%]\tMinimum depth\tGene")
for i in range(len(final_result)):
    print("{:.2f}\t\t{:d}\t\t{}".format(final_result[i][1] * 100, final_result[i][2], final_result[i][0]))

time_end = time.time()
time_elapsed = time_end - time_start
time_elapsed
