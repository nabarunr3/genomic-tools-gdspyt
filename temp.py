def frame_orfs(frame):
    #define the start and stop codons
    start_codon = "ATG"
    # print("Calculating all ORFs in frame of length", len(frame), "...")
    stop_codon_tuple = ("TAA", "TAG", "TGA")

    #trimming the frame (from the end) to number of nucleotides in multiples of 3
    trimmed_frame = ""
    for i in range((len(frame) - (len(frame)%3))):
        trimmed_frame = trimmed_frame + frame[i]

    #now get all the positions of start and stop codons
    start_pos_list = []
    stop_pos_list = []

    #initialize the while loop variable
    h = 0
    while (h <= (len(trimmed_frame) - 3)):
        codon = trimmed_frame[h:h+3]

        #storing the start positions 
        if codon == start_codon:
            start_pos_list.append(h)

        #storing the stop positions 
        if codon in stop_codon_tuple:
            # print(codon)
            stop_pos_list.append(h)

        #incrementing i by 3
        h = h + 3

    # print(start_pos_list)
    # print(stop_pos_list)

    #Now getting the ORFs.     
    orf_dict = {}
    i = 1 #initializing a variable to store the ORF number
    for start_position in start_pos_list:
        # a new list of stop positions is to be built which will contain stop positions after the start position
        new_stop_pos_list = []
        for stop_position in stop_pos_list:
            if stop_position > start_position:
                new_stop_pos_list.append(stop_position)

        #print(new_stop_pos_list)
        for stop_posi in new_stop_pos_list:
            # print(trimmed_frame[stop_posi:stop_posi+3])
            orf_key = "ORF"+str(i)
            orf_dict[orf_key] = [[],
                                 start_position + 1, #+1 to start indexing from 1
                                 stop_posi + 3]
            for j in range(start_position, stop_posi+1):
                orf_dict[orf_key][0].append(trimmed_frame[j:j+3])
            #updating ORF number 
            i = i + 1

    return orf_dict

#This is a simple function to reverse a string
def sequence_reverse(sequence):
    """
    This program finds the reverse of the sequence entered. 
    """
    #first we initialize the string variable which will contain our reversed string
    reversed_sequence = ''
    #next we initialize our loop variable to the value of the last index of 'sequence'
    i=len(sequence)-1
    #we use a while loop, which continues till i is greater than or equal to zero 
    while i>=0:
        #we concatenate each letter from the end of the sequence entered to form a new string which has the order reversed
        reversed_sequence = reversed_sequence+sequence[i]
        #in each iteration we decrease the value of the variable i to maintain the loop condition
        i=i-1

    #finally, exit from the function and return the reversed sequence    
    return reversed_sequence

#this function returns the complement of a DNA sequence entered, NOT the reverse complement
def dna_complement(sequence):
    """
    This function returns the complementary sequence of a DNA sequence entered.

    """
    #first we define a dictionary containing all the bases and their complements, even in lowecase
    bases_dict = {"A":"T", "a":"t", "C":"G", "c":"g", "G":"C", "g":"c", "T":"A", "t":"a", "N":"N", "n":"n"}

    #we initialize a string variable to store the complementary sequence
    complement = ''

    #next we loop over our sequence and make a new complementary string
    for i in sequence:
        complement = complement + bases_dict[i]

    #finally, return the complement sequence and exit from the function
    return complement

#a simple function which calls the reverse string and sequence complement functions to calculate the reverse complements
def reverse_complement(sequence):
    return(dna_complement(sequence_reverse(sequence)))

def frame_extract(sequence):
    """This function extracts all the 6 ORFs from a given sequence."""

    #initializing the dictionary which will store all the sequences of all the frames 
    frame_dictionary = {}

    #this variable stores the frame numbers 
    frame_no = 1

    #loop extracting frames 1, 2 and 3
    for i in range(3):
        frame_name = "Frame" + str(frame_no)
        frame_dictionary[frame_name] = sequence[i:]
        frame_no = frame_no + 1

    # print(frame_no)

    sequence_reverse_complement = reverse_complement(sequence)

    #loop extracting frames 4, 5, and 6
    for i in range(3):
        frame_name = "Frame" + str(frame_no)
        frame_dictionary[frame_name] = sequence_reverse_complement[i:]
        frame_no = frame_no + 1

    # print(frame_no)

    return(frame_dictionary)

def sequence_orfs(sequence):
    """This function finds and returns all the ORFs in all frames of a sequence."""
    #initializing the dictionary to store all the frame information of the sequence 
    sequence_orfs_dictionary = {}

    #calling the frame_extract() to find all the frames and store it in frames_dictionary variable 
    frames_dictionary = frame_extract(sequence)

    #this loop goes through all the frames in the dictionary, calls the frame_orfs() function for each frame and stores the output in frames_dictionary.    
    for frame, frame_sequence in frames_dictionary.items():
        sequence_orfs_dictionary[frame] = frame_orfs(frame_sequence)

    #for frame 5, 6, and 7, which are the reverse complement frames, we need to alter the start and stop positions of the frame, so that they are with respect to the forward sequence

    reverse_frames_tuple = ("Frame4", "Frame5", "Frame6")

    #calculating the length of the sequence
    seqlen = len(sequence)

    # now altering the start and stop positions 
    for frame, ORF_dict in sequence_orfs_dictionary.items():
        if frame in reverse_frames_tuple:
            for orfs, info in ORF_dict.items():
                info[1] = seqlen - info[1]
                info[2] = seqlen - info[2]

    return sequence_orfs_dictionary

def print_sequence_orfs(sequence):
    """This function prints the ORF information of a sequence in a presentable format."""

    seq_orf_info = sequence_orfs(sequence)

    for frame, frame_dict in seq_orf_info.items():
        print("** _ _ _ _ _ _", frame, "_ _ _ _ _ _")
        print("Retrieving", len(frame_dict), "sequences in", frame, "...")
        for orf, orf_dict in frame_dict.items():
            print("***", orf, "\n", orf_dict[0][0], orf_dict[0][1], "...", orf_dict[0][-2], orf_dict[0][-1])
            print(" ", orf_dict[1], "           ", orf_dict[2])

    return 0

def compare_seq_orfs(sequence):
    """This function compares all the ORFs in a sequence."""

    #getting the ORF information of the sequences 
    seq_orf_info = sequence_orfs(sequence)

    #a dictionary to store the lengths of all the ORFs defined in seq_orf_info
    seq_orf_len_dict = {}

    #this loop builds the seq_orf_len_dict as well as finds the maximum length
    max_len = 0
    for frame, frame_dict in seq_orf_info.items():
        seq_orf_len_dict[frame] = {}
        for orf, orf_list in frame_dict.items():
            #we subtract the starting and ending positions of the ORFs
            orf_len = abs(orf_list[1] - orf_list[2])
            #initialize a new dictionary inside of seq_orf_len_dict, corresponding to a key, frame
            seq_orf_len_dict[frame][orf] = orf_len
            if orf_len >= max_len:
                max_len = orf_len

    #now we find the minimum length 
    min_len = max_len 
    for frame, frame_dict in seq_orf_len_dict.items():
        for orf, orf_len in frame_dict.items():
            if orf_len <= min_len:
                min_len = orf_len


    #since there can be more than one sequences with the maximum or minimum lengths, we initialize a dictionary which will store all frames and ORFs corresponting to the maximum and minimum lengths 
    max_len_dict = {}
    min_len_dict = {}
    for frame, frame_dict in seq_orf_len_dict.items():
        for orf, orf_len in frame_dict.items():
            if orf_len == max_len:
                max_len_dict[frame] = {}
                max_len_dict[frame][orf] = [seq_orf_info[frame][orf][1], seq_orf_info[frame][orf][2]]
            elif orf_len == min_len:
                min_len_dict[frame] = {}
                min_len_dict[frame][orf] = [seq_orf_info[frame][orf][1], seq_orf_info[frame][orf][2]]

    return (max_len, max_len_dict, min_len, min_len_dict, seq_orf_len_dict)
