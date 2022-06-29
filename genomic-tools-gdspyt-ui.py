from orfs import *
comparison_tuple = compare_seq_orfs(test_sequence, 2)
print("Max length:", comparison_tuple[0])
print("Sequences with max length:", comparison_tuple[1])
print("Min length:", comparison_tuple[2])
print("Sequences with min length:", comparison_tuple[3])
print("All the ORFs with their starting and stop positions and their length:", comparison_tuple[4])

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

print("* Sequence information of the file:")
print("\nThere are a total of", len(fast_dict), "sequences in the file.")

print("\nThe following identifiers have been retrieved:\n")
print("Identifier\t\t\t\t Length")
print("__________\t\t\t\t ______")
for keys, sequence in fast_dict.items():
    print(keys, "\t", len(sequence))

from seqlen_compare import *
seq_len_compare(fast_dict)

from orfs import *

print("* ORFs in the file")
print_orfs_of_file(fast_dict)

print("* Longest and Shortest ORFs in the file.")
print_longest_shortest_ORFs(fast_dict)
print("* Longest and Shortest ORFs in individual frames.")
i = 1
while (i <= 3):
    print("** Frame"+str(i))
    print_longest_shortest_ORFs(fast_dict, i)
    i = i + 1

from repeats import *

def repeat_times_file(fasta_dict, repeat):
    """This returns the the frequency of occurance of a repeat in the file."""

    frequency = 0
    #this loop simply adds the number of times the repeat occurs in each sequence of fasa_dict
    for ID, sequence in fasta_dict.items():
        frequency = frequency + get_repeat_frequency(sequence, repeat)

    return frequency

input_repeat_unit = 'CGCGCTCG'
#print("The frequency of", input_repeat_unit, "is", repeat_times_file(fast_dict, input_repeat_unit))

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

#print_max_repeat(fast_dict, 8)
