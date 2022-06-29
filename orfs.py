def frame_orfs(frame):
    """
    This function extracts all non-overlapping ORFs from a sequence, 
    that is, ignoring start codons that fall between 
    a previous start codon and stop codon.

    """

    #define the start and stop codons
    start_codon = "ATG"
    stop_codon_tuple = ("TAA", "TAG", "TGA")

    #initialize the dictionary to store the orfs, in the form of
    # orf_dictionary = {
    #     ORF1:[["start_codon, ...., stop_codon"]
    #           starting_position,
    #           stopping_position],
    #     ORF2:[["start_codon, ...., stop_codon"]
    #           starting_position,
    #           stopping_position],
    #     .
    #     .
    #     .
    #     ORFn:[["start_codon, ...., stop_codon"]
    #           starting_position,
    #           stopping_position]}
    orf_dict = {}

    #trimming the frame (from the end) to number of nucleotides in multiples of 3
    trimmed_frame = ""
    for i in range((len(frame) - (len(frame)%3))):
        trimmed_frame = trimmed_frame + frame[i]

    #this will keep track of the orf number
    orf_count = 0  
    h = 0
    while (h <= (len(trimmed_frame) - 3)):
        codon = trimmed_frame[h:h+3]

        if codon == start_codon:
            start_position = h 

            #we move to check from the next codon, for a stop codon 
            h = h + 3
            #we now enter the stop codon searching loop 
            while (h <= (len(trimmed_frame) - 3)):
                codon = trimmed_frame[h:h+3]

                if codon in stop_codon_tuple:
                    stop_position = h 
                    orf_count = orf_count + 1 

                    #to assign a name to the ORF, which will act as a key in the orf_dict. 
                    orf_name = "ORF" + str(orf_count) 
                    orf_dict[orf_name] = []

                    #now we store the ORF information between start and stop positions in orf_dict
                    #first we store the list of codons in a nested list in the first position of the list corresponding to the orf_number 
                    orf_dict[orf_name].append([])
                    j = start_position
                    while(j <= stop_position):
                        codon = trimmed_frame[j:j+3]
                        orf_dict[orf_name][0].append(codon)
                        j = j + 3
                    #next we store the start and stop positions in the second and third positions of the orf_name list, respectively 
                    orf_dict[orf_name].append(start_position)
                    #we add 3 to the stop position to consider the 3rd nucleotide of the stop codon as the stop position
                    orf_dict[orf_name].append((stop_position + 3)) 

                    #we now break out of the stop codon searching loop 
                    break

                #if the codon is not a stop codon, we increment h by 3 to search for another codon 
                h = h + 3

        #if the codon is not a start codon, to move to the next test codon, have increment h by 3
        h = h + 3 

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

def print_orfs_of_file(fasta_dict):
    """This function takes as input a dictionary containing all sequences in a file in the format "identifier:sequence" and prints all the ORFs of all reading frames of all sequences of a file."""

    for identifier, sequence_orf_info in fasta_dict.items():
        print("* ++++++", identifier, "++++++")
        print_sequence_orfs(sequence_orf_info)

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

def get_all_ORFs_of_length(fasta_dict, x):
    """This function gets all ORFs of length x from a file whose sequence information has been stored in fasta_dict"""
    all_ORFs_dict = {}

    file_ORF_lengths_dict = {}
    for identifier, sequence in fasta_dict.items():
        seq_info = compare_seq_orfs(sequence)
        file_ORF_lengths_dict[identifier] = seq_info[4]


    #initializing the list containing ORF of length x
    dict_file_orfs_length_x = []
    for identifier, identifier_info in file_ORF_lengths_dict.items():
        for frame, frame_ORF_info in identifier_info.items():
            for ORF, length in frame_ORF_info.items():
                if length == x:
                    ORF_coordinates = (identifier, frame, ORF)
                    dict_file_orfs_length_x.append(ORF_coordinates)


    return dict_file_orfs_length_x

def print_longest_shortest_ORFs(fasta_dict):
    """This function prints all the longest and shortest ORFs, considering all the frames of all sequences in a file """

    #first we call the file_ORF_compare_lengths() to store the maximum and minimum ORF lengths of the file in a tuple returned by the function 
    max_min_tuple = file_ORF_compare_lengths(fasta_dict)
    #then we store the all the ORFs having the maximum and minimum lengths
    max_ORF_list = get_all_ORFs_of_length(fasta_dict, max_min_tuple[0])
    min_ORF_list = get_all_ORFs_of_length(fasta_dict, max_min_tuple[1])

    #Finally to print them out in a nicely formatted way 
    print("\n\nThe ORF(s) with the maximum lengths of", max_min_tuple[0], "are:")
    print("\nIdentifier \t\t\t Frame \t\t ORF")
    print("__________ \t\t\t _____ \t\t ___")
    for i in range(len(max_ORF_list)):
        print(max_ORF_list[i][0], "\t", max_ORF_list[i][1], "\t", max_ORF_list[i][2])
    print("\n\nThe ORF(s) with the minimum lengths of", max_min_tuple[1], "are:")
    print("\nIdentifier \t\t\t Frame \t\t ORF")
    print("__________ \t\t\t _____ \t\t ___")
    for i in range(len(min_ORF_list)):
        print(min_ORF_list[i][0], "\t", min_ORF_list[i][1], "\t", min_ORF_list[i][2])

    return 0
