def seq_len_compare(seq_dict):
    """This function compares the length of DNA sequences in a dictionary
    """

    """first, we create a new sequence dictionary which stores lengths against identifiers
    """
    seq_len_dict = {}
    for keys,values in seq_dict.items():
        seq_len_dict[keys] = len(values)

    #for testing     
    # for keys,values in seq_len_dict.items():
    #     print(keys,values)


    """
    Now we compare the lengths of sequences in the dictionary.
    """

    #1. Maximum Length

    #first we initialize the maximum sequence length to 0
    max = 0
    #next we initialize a variable to store the identifier corresponding to the maximum sequence length
    max_identifier = ''
    #now we create a list which holds all identifiers which have the maximum value
    max_list = []
    for keys,values in seq_len_dict.items():
        if values>max:
            max = values
            max_identifier = keys

    #next we store all the identifiers in the dictionary which have the value corresponding to the maximum length
    for keys,values in seq_len_dict.items():
        if values==max:
            max_list.append(keys)

    #2. Minimum Length 
    #first we initialize the minimum sequence length to the maximum value obtained in the previous section
    min = max 
    #next we initialize a variable to store the identifier corresponding to the minimum sequence length
    min_identifier = ''
    #now we create a list which holds all identifiers which have the minimum value
    min_list = []
    for keys,values in seq_len_dict.items():
        if values<min:
            min = values
            min_identifier = keys

    #next we store all the identifiers in the dictionary which have the value corresponding to the minimum length
    for keys,values in seq_len_dict.items():
        if values==min:
            min_list.append(keys)

    """
    Now we print out the values obtained along with their identifiers.
    """
    #1. Print Maximum Length(s)
    print("\nThe following", len(max_list),  "sequence(s) has/have the maximum length of", max, '.')
    for x in range(len(max_list)):
        #print the identifiers in the max_list
        print("\n", max_list[x])

    #2. Print Minimum Length(s)
    print("\nThe following", len(min_list), "sequence(s) has/have the minimum length of", min,'.')
    for x in range(len(min_list)):
        #print the identifiers in the min_list
        print("\n", min_list[x])

    #to exit from the function
    return 0
