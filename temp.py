test_sequence="AGGCTTGACCTTGTGCGGCAGCCGGCAGCCGCCGGCCGGCAGGTTGCCGTACGGATCCTCGTCGGCCGCCGTCGCGGGGGCCGCCGGCATCACGGGGGCGGTCGACGCGACCAGCATCGGGGGGAACGGCAACGCATGCACCCGCTCGCCGGTCAGCACGAACGCACCGTTGGCCACCGCCGGCGCGATCGGCGGCACGCCGGCGTCGCTCAACCCGGTCGGCTGGGCGTCCGACGGCACGAAGAAGACATCCACCGGCGGCGCTTCCTGCATGCGTATCGGCGAATAGTCGGCGAAACCGGCGTTGCGGACCGCGCCATGGTCGACGTCGATCGCGAAGCCGGGCTTCGTCGTCGCGAGACCGAACAGCGCGCCGCCCTGGATCTGCGCTTGCGCGCCGGTCGGGTTGACGATGCGGCCCGCATACACGCCGGCCGTCACGCGATGCACGCGCGGTTGTTGCGCTTCGATCGACACTTCCGTCACGTACGCGACGACCGAGCCGGCCGTTTCGTGCATCGCGACGCCCCACGCGTGCCCGGCCGGCAGCGTGCGCGCGCCGTAGCCGGACTTGTCGACGGCCAGCGCGAGCGCCTGCCGATGCGCGGCGTGCTCGGGGCCGGCCAGCCGCGTCATCCGGTAGGCGACCGGATCCTGCCGCGCCGAGTGCGCGAGCTCGTCGACCAGCGTTTCCATCACGAACGCCGTATGCGAGTTGCCGCCCGAGCGCCACGTCTGGACCGGCACGTCGGCCTCGGTCTGATGAACCGATACCTGCATCGGGAAGCCGTACGGGCTGTTCGTCACGCCTTCGGTCAGGCTCGGATCGGTGCCGCGCTTGAGCATCGTCGTGCGCTCGAGCGGCGAGCCCTTCAGCACAGACTGGCCGACGACCACGTGCTGCCAGTCGCGCACGGCGCCGCTGCCGTCCACGCCGATGTCGACGCGATGCAGCACCATCGGGCGGTAATAGCCGCCGCGCAGATCGTCCTCGCGCGTCCAGATCGTCTTGACGGGGCCGAGATGGCCGGCCGCGAGGTACGCGGCGGACACGTGGGCGGCTTCGACCACGTAGTCCGACGTCGGCGTCGAGCGCCGGCCATAGTCGCCGCCCGAGGTCAGCGTGAAGATCTGGACTTTCTCCGGGGCGACGCCGAGCGCCTTCGCGACCGCCGCGCGGTCGGTCGTC"

def frame_orfs(frame):
    #define the start and stop codons
    start_codon = "ATG"
    print(len(frame))
    stop_codon_tuple = ("TAA", "TAG", "TGA")

    #trimming the frame (from the end) to number of nucleotides in multiples of 3
    trimmed_frame = ""
    for i in range((len(frame) - (len(frame)%3))):
        trimmed_frame = trimmed_frame + frame[i]

    #now get all the positions of start and stop codons
    start_pos_list = []
    stop_pos_list = []
    for i in range((len(trimmed_frame) - 3)):
        codon = trimmed_frame[i:i+3]

        #storing the start positions 
        if codon == start_codon:
            start_pos_list.append(i)

        #storing the stop positions 
        if codon in stop_codon_tuple:
            stop_pos_list.append(i)

        #incrementing i by 3
        i = i + 3

    print(start_pos_list)
    print(stop_pos_list)

    #Now getting the ORFs.     
    orf_dict = {}
    for start_position in start_pos_list:
        i = 1 #initializing a variable to store the ORF number
        # a new list of stop positions is to be built which will contain stop positions after the start position
        new_stop_pos_list = []
        for stop_position in stop_pos_list:
            if stop_position > start_position:
                new_stop_pos_list.append(stop_position)

        print(new_stop_pos_list)

        for stop_position in new_stop_pos_list:
            #print(stop_position)
            orf_key = "ORF"+str(i)
            orf_dict[orf_key] = [[],
                                      start_position,
                                      stop_position]
            for j in range(start_position, stop_position):
                orf_dict[orf_key][0].append(trimmed_frame[j:j+3])

        #updating ORF number 
        i = i + 1
    return orf_dict

print(frame_orfs(test_sequence))
