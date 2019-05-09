
#from functions2 import *

import time, gzip       # time measurement
start = time.time() # time measurement

### my files  - VLad
infile_read1 = "file1.txt.gz"
infile_read2 = "file2.txt.gz"
infile_res = "resistance_genes.fsa" 

def f_kmer_generator(sequence, kmer_length):
    """k-mer generator using list"""
    kmer_list = list()
    for i in range(0, len(sequence)-kmer_length+1):
        kmer_list.append(sequence[i:i+kmer_length])
    return(kmer_list)
    
def reverse_complement(DNA_sequence):
    """Reverse complement a given DNA sequence"""
    DNA_translation_table = str.maketrans("ACGT", "TGCA")
    rev_compliment = DNA_sequence.translate(DNA_translation_table)
    return rev_compliment[::-1]    

def f_ResGeneDict_Kmers(infile):
    ID_list = []
    SEQ_list = []
    RES_KMER_set = set ()
    ResGene_Dict = {}
    seq = ""
    """ Create a dictionary of all resistace genes (keys) and their sequences (values,
    and a nested dictionary of positions (keys) and counts (as values)"""
    print("Current function f_ResGeneDict")
    ### Needs a TRY statement_________________________________________________
    
    # Kmer creation from Resistance Genes file
    ResFile = open(infile, "r")
    ### _____________________________________________________________________
 
    for line in ResFile:                                                        # When ">" , assemble the sequence
        if line.startswith('>'):                                                # Read Line by line
            if seq:   # push to seqs
                ResGene_Dict[info] = {}
                ResGene_Dict[info][seq]= {}                
                SEQ_list.append(seq)
                for i in f_kmer_generator(seq, 19):
                    RES_KMER_set.add(i)
                    ResGene_Dict[info][seq][i]=0
                seq = ''
            info = line[:-1]
            ID_list.append(info)      # append the ID      
        else:
            seq += line.strip()   
                                                                                  
    ResGene_Dict[info] = {}                                                     # Last entry push
    ResGene_Dict[info][seq]= {}
    
    ### Last entry
    #ID_list.append(line[:-1])                                                      # Last entry push
    SEQ_list.append(seq)
    for i in f_kmer_generator(seq, 19):
        RES_KMER_set.add(i)
        ResGene_Dict[info][seq][i]=0
    print("f_ResGeneDict_Kmers of %s is complete!" % infile) 
    return RES_KMER_set, ResGene_Dict

def f_ResGeneDict_Kmers_old(infile):
    """ Create a dictionary of all resistace genes (keys) and their sequences (values,
    and a nested dictionary of positions (keys) and counts (as values)"""
    print("Current function f_ResGeneDict")
    ### Needs a TRY statement_________________________________________________
    
    # Kmer creation from Resistance Genes file
    ResFile = open(infile, "r")
    ### _____________________________________________________________________
    
    info = [] # info line
    ResGene_Dict = {}
    data = ''    
    for line in ResFile:                                                        # When ">" , assemble the sequence
        if line.startswith('>'):                                                # Read Line by line
            if data:   # push to seqs
                ResGene_Dict[info] = {}
                ResGene_Dict[info][data]= {}
                for i in f_kmer_generator(data, 19):
                    ResGene_Dict[info][data][i]=0
                data = ''
            info = line[:-1]           
        else:
            data += line.strip()                                                                                     
    ResGene_Dict[info] = {}                                                     # Last entry push
    ResGene_Dict[info][data]= {}
    for i in f_kmer_generator(data, 19):
        ResGene_Dict[info][data][i]=0
    print("f_ResGeneDict_Kmers of %s is complete!" % infile) 
    return ResGene_Dict

def f_gzip_extract_reads(filename1, filename2):
    """ Opening a gzip file and extracting the reads"""
    print("Initialized: sequence exraction from files <%s , %s>" %(filename1,filename2))
    Kmer_Dict = dict() # TRy DICT from teh start
    seq_counter=0 
    line_count=2
    with gzip.open(filename1, 'rt') as f1, gzip.open(filename2, "rt") as  f2:
        while True:
            
            line_1 = f1.readline()
            line_2 = f2.readline()
            line_count += 1        
            if line_count % 4 == 0:                                             # grab the sequence line
                
                read1 = line_1.rstrip()
                read2 = line_2.rstrip()
                
                read_reverse1 = reverse_complement(read1)  
                read_reverse2 = reverse_complement(read2)     
                       
                f_kmer_generator(Kmer_Dict, read1,19)
                f_kmer_generator(Kmer_Dict, read2,19)
                
                f_kmer_generator(Kmer_Dict, read_reverse1,19)
                f_kmer_generator(Kmer_Dict, read_reverse2,19)
                
            if (line_count-2)*2 % 1000000== 0:                                      # Counter , two files, started from 2
                print("# of sequences processed: ",line_count-2) 
                break
            if not line_1 and not line_2:
                break     

    f1.close()
    f2.close()
    Kmer_Dict = {key:val for key, val in Kmer_Dict.items() if val >= 10}
    print("Reading of the following files is complete! \n %s \n %s"  % (infile_read1, infile_read2))
    return Kmer_Dict




def f_reads_subselect(infile_read1,infile_read2):  
    print("Initialized: sequence exraction from files <%s , %s>" %(infile_read1,infile_read1))
    READS_list =  []
    line_track = [] 
    FIND_list = []
    seq_counter=0 
    line_count=2
    line1,line2 =0,0
    with gzip.open(infile_read1, 'rt') as f1, gzip.open(infile_read2, "rt") as  f2:
        while True:
            
            line_1 = f1.readline()
            line_2 = f2.readline()
            line_count += 1  
            line1+=1
            line2 +=1
            if line_count % 4 == 0:                                             # grab the sequence line
                
                    read1 = line_1.rstrip()
                    read1_RC = reverse_complement(read1)
                    read2 = line_2.rstrip()
                    read2_RC = reverse_complement(read2)
                    
                    if read1[:19] in RES_KMER_set or read1[-19:] in RES_KMER_set:
                        READS_list.append(read1)
                        line_track.append(line_count)
                        
                        
                    elif read1_RC[:19] in RES_KMER_set or read1_RC[-19:] in RES_KMER_set:
                        READS_list.append(read1)
                        line_track.append(line_count)
                        
                    if read2[:19] in RES_KMER_set or read2[-19:] in RES_KMER_set:
                        READS_list.append(read2)
                        line_track.append(line_count)
                        
                    elif read2_RC[:19]in RES_KMER_set or read2_RC[-19:] in RES_KMER_set:
                        READS_list.append(read2)
                        line_track.append(line_count)
    
            if (line_count-2)*2 % 100000== 0:                                      # Counter , two files, started from 2
                print("# of sequences processed: ",line_count-2) 
                
            if not line_1 and not line_2:
                break     
    
        f1.close()
        f2.close()
        return READS_list

def f_Dict_crosscheck2(Res_Dict, Reads_subset):
    ### So here i can quickly check for the presence of the kmers in both dics
    ### If we construct a loop to go over each kmer for reference gene in turn and if complete
    ### create a file with the gene ID , coverage depth value ( which is set for the subselection)
    ### That could be a good output as a start??? 
    
    """ Cross check reads and resistance genes dictionaries, 
    check if the kmers are present in both and assign the count from 
    the reads dictionary"""
    results = []
    try: 
        gene_counter = 0
        #inner = 0
        

        for key1 in Res_Dict:                                                 # for ID
            for key2 in Res_Dict[key1]:                                       # for Sequence
                for key3 in Res_Dict[key1][key2]:      
                    if key3 in Reads_subset.keys():
                        newRes_Dict[key1][key2][key3]=newRes_Dict[key1][key2][key3]+Reads_subset[key3]
                        
                        
                        
                        kmer_count +=1                                           # if kmer matched +=1
                        newRes_Dict[key1][key2][key3]=Reads_subset[key3]
                        if min_depth is None or Reads_subset[key3] < min_depth:
                            min_depth = Reads_subset[key3]
                gene_coverage = int(kmer_count*100/seq_len)  
                if gene_coverage >=95:
                    print("%s : \nCoverage of %s %% and minimum depth of %s\n" % (key1,gene_coverage, min_depth))
                    results.append((key1,gene_coverage, min_depth))
            if gene_counter % 100 == 0:
                print("Number of genes processed: ",gene_counter)
    except TypeError as e:
        print(e)        
        print(key1, key2, key3,)

    return results   
 
def f_Dict_crosscheck2(newRes_Dict, Reads_subset):
    ### So here i can quickly check for the presence of the kmers in both dics
    ### If we construct a loop to go over each kmer for reference gene in turn and if complete
    ### create a file with the gene ID , coverage depth value ( which is set for the subselection)
    ### That could be a good output as a start??? 
    
    """ Cross check reads and resistance genes dictionaries, 
    check if the kmers are present in both and assign the count from 
    the reads dictionary"""
    results = []
    try: 
        gene_counter = 0
        #inner = 0
        

        for key1 in newRes_Dict: 
                                                # for ID
            gene_counter +=1
            
            for key2 in newRes_Dict[key1]:                                       # for Sequence
                kmer_count = 0                                                   # how many kmers has been matched
                seq_len = len(key2)-18                                           # how many kmers per gene (Seq_len - kmer_size_len+1 (0-indexing))
                min_depth = None
                
                for key3 in newRes_Dict[key1][key2]:
                                               # the inner key aka kmers
                    if key3 in Reads_subset:
                        kmer_count +=1                                           # if kmer matched +=1
                        newRes_Dict[key1][key2][key3]=Reads_subset[key3]
                        if min_depth is None or Reads_subset[key3] < min_depth:
                            min_depth = Reads_subset[key3]
                gene_coverage = int(kmer_count*100/seq_len)  
                if gene_coverage >=95:
                    print("%s : \nCoverage of %s %% and minimum depth of %s\n" % (key1,gene_coverage, min_depth))
                    results.append((key1,gene_coverage, min_depth))
            if gene_counter % 100 == 0:
                print("Number of genes processed: ",gene_counter)
    except TypeError as e:
        print(e)        
    return results  

def f_Kmer_depth_rev_com(Reads_list,min_depth):
    """ Create all kmers from reads and subselect 
    for the depth(aka number of occurance"""
    print("Curent function f_Kmer_depth")
    test_counter=0    
    Kmer_Dict = dict()
    for read in Reads_list:
        test_counter +=1
        
        kmer_list = f_kmer_generator(read,19)
        rev_read = reverse_complement(read)                                       # Reverse complement the sequence 
        rev_kmer_list = f_kmer_generator(rev_read,19)
        kmer_list.extend(rev_kmer_list)
        for kmer in kmer_list:
            if kmer in Kmer_Dict:
                Kmer_Dict[kmer] += 1
            else:
                Kmer_Dict[kmer] = 1
        if test_counter%100000 ==0:
            print("f_Kmer_depth, counter state: ", test_counter ) 
    Reads_subset = {key:val for key, val in Kmer_Dict.items() if val > min_depth} # Create a subset based on depth specifications
    print("Done") 
    return Reads_subset # this can be removed if we choose to subselect (will free up the memory )
 ### currently with these files it is .... 123 445 562 kmers in total, subset is 11 138 671 (x100 less if cut to min depth 10)
  
 
RES_KMER_set, ResGene_Dict = f_ResGeneDict_Kmers(infile_res)
READS_list=f_reads_subselect(infile_read1,infile_read2)
READS_byDEPTH_dict = f_Kmer_depth_rev_com(READS_list,10)
results = f_Dict_crosscheck2(ResGene_Dict, READS_byDEPTH_dict)

end = time.time() # time measurement

time_elapced= (end-start)/60
print(time_elapced)
# Write the results
outfile = open("ResFinder_older.txt", "w")
results.sort(reverse = True,key=lambda x: (x[1], x[2]))
for i in range(len(results)):
    print("{:.2f}\t\t{:.2f}\t\t{}".format(results[i][1], results[i][2], results[i][0]), file = outfile)
print("Processing time elapsed :",int((end-start)/60), "min", file = outfile)
outfile.close()    

print("Processing time elapsed :",int((end-start)/60), "min")

end = time.time() # time measurement

time_elapced= (end-start)/60
print(time_elapced)
