def file_ORF_compare_lengths(fasta_dict):
 """This function compares the maximum and minimum lengths of all sequences in a file and returns the lengths of the longest and shortest ORF in the file."""
 #initializing the dictionary which will store the ORF lengths of all the sequences 
 file_ORF_len_dict = {}

 #building the file_ORF_len_dict dictionary 
 for identifier, sequence in fasta_dict.items():
  file_ORF_len_dict[identifier] = compare_seq_orfs(sequence)

 #now we will go over the entries of the dictionary and get the values of the maximum ORF
 max_len_file = 0
 for identifier, sequence_info in file_ORF_len_dict.items():
  #the first entry in the tuple returned by compare_seq_ORFs contains the maximum length of the sequence 
  if sequence_info[0] > max_len_file:
   max_len_file = sequence_info[0]

 #getting the value of the minimum length of the ORF
 min_len_file = max_len_file
 for identifier, sequence_info in file_ORF_len_dict.items():
  #the third entry in the tuple returned by compare_seq_ORFs contains the minimum length of the sequence 
  if sequence_info[2] < min_len_file:
   min_len_file = sequence_info[2]

 return (max_len_file, min_len_file)

#!/usr/bin/python3
import sys
def usage():
    print
    """
    The program allows the user to ....

    Basic Usage:

    python3 fastinfo.py <filename.FASTA> [option1] [option2...]

    """
    return 0

#now we open the file, which is stored in the second position of the sys.argv[] list
try:
    f = open(sys.argv[1], 'r')
    print("FASTA File found. Processing...")
except IOError:
    print("The file does not exist! Please recheck the filename.")
    sys.exit()

from parsefasta import *
fast_dict = multi_FASTA_to_dict(f)

print("\nThere are a total of", len(fast_dict), "sequences in the file.")

print("\nThe following identifiers have been retrieved:")
print("Identifier \t\t\t Length")
for keys, sequence in fast_dict.items():
    print(keys, "\t", len(sequence))

from seqlen_compare import *
seq_len_compare(fast_dict)

from orfs import *

print_orfs_of_file(fast_dict)

print_longest_shortest_ORFs(fast_dict)

from repeats import *

def repeat_times_file(fasta_dict, repeat):
    """This returns the the frequency of occurance of a repeat in the file."""

    frequency = 0
    #this loop simply adds the number of times the repeat occurs in each sequence of fasa_dict
    for ID, sequence in fasta_dict.items():
        frequency = frequency + get_repeat_frequency(sequence, repeat)

    return frequency

input_repeat_unit = 'CGCGCTCG'
print("The frequency of", input_repeat_unit, "is", repeat_times_file(fast_dict, input_repeat_unit))

def print_max_repeat(fasta_dict, n):
    """This function prints the maximum number of times of occurance of repeat units of lengths n and all such repeat units with occur with the maximum frequency"""

    #first we get the max repeats list by calling the max_occurance_n function. Remember that the first position of the list returned by the function contains the maximum frequency
    max_repeats_list = max_occurance_n(fasta_dict, n)

    print("The repeat of length", n, "occurs a max of", max_repeats_list[0], "times in the file.")
    print("The following repeats occur at the maximum frequency:")
    i = 1 #because we want to print from the second position of the list
    list_len = len(max_repeats_list)
    while(i < list_len):
        print(max_repeats_list[i])
        i = i + 1 

    return 0

print_max_repeat(fast_dict, 8)
