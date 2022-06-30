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

print("* Identifier ORFs")
identifier = "gi|142022655|gb|EQ086233.1|43"
tester = fast_dict[identifier]
print("The identifier", identifier, "has ORFs with the following:\n")
print_seq_orf_info(tester)

from repeats import *

input_repeat_unit = 'GCCGCCG'
print("* Repeat Information")
print("\nThe repeat", input_repeat_unit, "occurs with a frequency of", repeat_times_file(fast_dict, input_repeat_unit))

print_max_repeat_n(fast_dict, 11)
