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

# Printing length information of sequences

from seqlen_compare import *
seq_len_compare(fast_dict)

#first we got to import the functions which store sequence lengths of a file in the form of a dictionary.
from temp import * 

def print_orfs_of_file(fasta_dict):
    """This function takes as input a dictionary containing all sequences in a file in the format "identifier:sequence" and prints all the ORFs of all reading frames of all sequences of a file."""

    for identifier, sequence_orf_info in fasta_dict.items():
        print("* ++++++", identifier, "++++++")
        print_sequence_orfs(sequence_orf_info)

print_orfs_of_file(fast_dict)

def file_ORF_compare(fasta_dict):
 #initializing the dictionary which will store the ORF lengths of all the sequences 
 file_ORF_len_dict = {}

 #building the file_ORF_len_dict dictionary 
 for identifier, sequence in fasta_dict.items():
  file_ORF_len_dict[identifier] = compare_seq_orfs{fasta_dict}

 #now we will go over the entries of the dictionary and get the values of the maximum ORF
 max_len_file = 0
 for identifier, sequence_info in file_ORF_len_dict.items():
  for i in range(len(sequence_info)):
   if sequence_info[0] > max_len_file:
    max_len_file = sequence_info[0]

 #getting the value of the minimum length of the ORF
 min_len_file = max_len_file
 for identifier, sequence_info in file_ORF_len_dict.items():
  for i in range(len(sequence_info)):
   if sequence_info[2] < min_len_file:
    min_len_file = sequence_info[2]

#storing ORFs with the maximum and minimum lengths 
max_len_dict = {}
min_len_dict = {}
 for identifier, sequence_info in file_ORF_len_dict.items():
  for i in range(len(sequence_info)):
    for frame, ORF_info in sequence_info[1].items():
     for ORF, positions in ORF_info.items():


 return 0
