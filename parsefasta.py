def multi_FASTA_to_dict(FASTA_file_object):
    """
    This function takes a file object of a FASTA file and puts the sequences in a dictionary. 
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
