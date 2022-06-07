#!/usr/bin/python3
#the following code snippet reads from a fasta file and puts the sequences in a dictionary
import sys
def usage():
    print
    """
    The program allows the user to ....

    Basic Usage:

    python3 fastinfo.py <filename.fasta> [option1] [option2...]

    """
    return 0

#now we try to open the file, which is stored in the second position of the sys.argv[] list
try:
    f = open(sys.argv[1], 'r')
    print("FASTA File found. Processing...")
except IOError:
    print("The file does not exist! Please recheck the filename.")
    sys.exit()

# Putting the sequences of the file in a dictionary

def multi_fasta_to_dict(fasta_file_object):
    """
    This function opens a FASTA file and puts the sequences in a dictionary. 
    """
    sequence_dictionary = {}
    for line in fasta_file_object:
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

# Printing the identifiers and the number of sequences
# Now we print the sequence identifiers in the dictionary and the number of entries in the multi-FASTA file. 

fast_dict = multi_fasta_to_dict(f)
print("\nThe following identifiers have been retrieved:")
for keys in fast_dict:
    print(keys)
    print(fast_dict[keys])
print("\nThere are a total of", len(fast_dict), "sequences in the file.")
