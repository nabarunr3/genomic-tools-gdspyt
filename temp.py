test_sequence="AGGCTTGACCTTGTGCGGCAGCCGGCAGCCGCCGGCCGGCAGGTTGCCGTACGGATCCTCGTCGGCCGCCGTCGCGGGGGCCGCCGGCATCACGGGGGCGGTCGACGCGACCAGCATCGGGGGGAACGGCAACGCATGCACCCGCTCGCCGGTCAGCACGAACGCACCGTTGGCCACCGCCGGCGCGATCGGCGGCACGCCGGCGTCGCTCAACCCGGTCGGCTGGGCGTCCGACGGCACGAAGAAGACATCCACCGGCGGCGCTTCCTGCATGCGTATCGGCGAATAGTCGGCGAAACCGGCGTTGCGGACCGCGCCATGGTCGACGTCGATCGCGAAGCCGGGCTTCGTCGTCGCGAGACCGAACAGCGCGCCGCCCTGGATCTGCGCTTGCGCGCCGGTCGGGTTGACGATGCGGCCCGCATACACGCCGGCCGTCACGCGATGCACGCGCGGTTGTTGCGCTTCGATCGACACTTCCGTCACGTACGCGACGACCGAGCCGGCCGTTTCGTGCATCGCGACGCCCCACGCGTGCCCGGCCGGCAGCGTGCGCGCGCCGTAGCCGGACTTGTCGACGGCCAGCGCGAGCGCCTGCCGATGCGCGGCGTGCTCGGGGCCGGCCAGCCGCGTCATCCGGTAGGCGACCGGATCCTGCCGCGCCGAGTGCGCGAGCTCGTCGACCAGCGTTTCCATCACGAACGCCGTATGCGAGTTGCCGCCCGAGCGCCACGTCTGGACCGGCACGTCGGCCTCGGTCTGATGAACCGATACCTGCATCGGGAAGCCGTACGGGCTGTTCGTCACGCCTTCGGTCAGGCTCGGATCGGTGCCGCGCTTGAGCATCGTCGTGCGCTCGAGCGGCGAGCCCTTCAGCACAGACTGGCCGACGACCACGTGCTGCCAGTCGCGCACGGCGCCGCTGCCGTCCACGCCGATGTCGACGCGATGCAGCACCATCGGGCGGTAATAGCCGCCGCGCAGATCGTCCTCGCGCGTCCAGATCGTCTTGACGGGGCCGAGATGGCCGGCCGCGAGGTACGCGGCGGACACGTGGGCGGCTTCGACCACGTAGTCCGACGTCGGCGTCGAGCGCCGGCCATAGTCGCCGCCCGAGGTCAGCGTGAAGATCTGGACTTTCTCCGGGGCGACGCCGAGCGCCTTCGCGACCGCCGCGCGGTCGGTCGTC"

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

frm_dic = frame_extract(test_sequence)
for keys, items in frm_dic.items():
    print(keys, "\n", items, "\n")
