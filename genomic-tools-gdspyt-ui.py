#!/usr/bin/python3
#the following code reads from a FASTA file and puts the sequences in a dictionary
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

# Putting the sequences of the file in a dictionary

def multi_FASTA_to_dict(FASTA_file_object):
    """
    This function opens a FASTA file and puts the sequences in a dictionary. 
    """
    sequence_dictionary = {}
    for line in FASTA_file_object:
        line = line.rstrip() #removing the trailing newlines from the FASTA file
        #print(line)
        if line[0] == '>':
            #split the header into a list consisting of the sequence identifier and the description 
            header_word_list = line[1:].split() 
            #identifier would be the first word of the header word list
            identifier = header_word_list[0]
            #print(identifier)
            #initializing an empty entry corresponding to the sequence identifier
            sequence_dictionary[identifier] = ""
            #print(sequence_dictionary[identifier])
        else:
            sequence_dictionary[identifier] = sequence_dictionary[identifier]+line #we concatenate the subsequent lines to the entry corresponding to an identifier
            #print(sequence_dictionary[identifier])

    return sequence_dictionary

# Printing file information 
# Now we print the sequence identifiers in the dictionary and the number of entries in the multi-FASTA file. 

fast_dict = multi_FASTA_to_dict(f)
print("\nThe following identifiers have been retrieved:")
for keys in fast_dict:
    print(keys)
    #print(fast_dict[keys]) #prints all sequences, hence commented out
print("\nThere are a total of", len(fast_dict), "sequences in the file.")

#first we got to import the functions which store sequence lengths of a file in the form of a dictionary.
from temp import * 

def print_orfs_of_file(fasta_dict):
    """This function takes as input a dictionary containing all sequences in a file in the format "identifier:sequence" and prints all the ORFs of all reading frames of all sequences of a file."""

    for identifier, sequence_orf_info in fasta_dict.items():
        print("* ++++++", identifier, "++++++")
        print_sequence_orfs(sequence_orf_info)

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

from repeats import *
print("\n\nThe following repeats of length", 8, "have been identified.\n")
print(get_file_repeat_units(fast_dict, 8))

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

    print("The repeat of length", n, "occurs", max_repeats_list[0], "times in the file.")
    print("The following repeats occur at the maximum frequency:")
    i = 1 #because we want to print from the second position of the list
    list_len = len(max_repeats_list)
    while(i < listlen):
        print(max_repeats_list[i])

    return 0

print_max_repeat(fast_dict, 8)
