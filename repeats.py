def get_repeat_frequency(sequence, repeat_unit):
    """This function returns the number of times a repeat unit occurs in a sequence."""

    repeat_count = 0
    #we start from the first position of the sequence and take out test windows corresponding to length of repeat unit.
    repeat_len = len(repeat_unit)
    seq_len = len(sequence)
    i = 0
    #we want slide the window till the last position in the sequence which can accomodate the windo
    while (i <= (seq_len - repeat_len)):
        test_window = sequence[i:(i + repeat_len)]
        if test_window == repeat_unit:
            repeat_count = repeat_count + 1

        i = i + 1

    return repeat_count

def get_seq_repeat_info(sequence, n):
    """This function gets a dictionary containing all subsequences of length n which occur repeatedly, that is, two times or more, in a sequence."""

    sequence_repeats_dict = {}

    #we start from the first position and extract the substrings of length l. We stop at the last last position minus the length of the repeat. 
    i = 0
    seqlen = len(sequence)
    while (i <= (seqlen - n)):
        test_repeat_unit = sequence[i:(i + n)]
        #now we get the frequency of this unit in our sequence
        frequency = get_repeat_frequency(sequence, test_repeat_unit)
        if (frequency > 1):
            sequence_repeats_dict[test_repeat_unit] = frequency

        i = i + 1 

    return sequence_repeats_dict

def get_file_repeat_units(fasta_dict, n):
    #dictionary storing the non-unique repeat information of the file by calling get_seq_repeat_info()
    file_repeat_info = {}
    for uid, sequence in fasta_dict.items():
        file_repeat_info[uid] = get_seq_repeat_info(sequence, n)

    #dictionary to store the times of occurance of each unique repeat unit in the file 
    file_unit_times_dict = {}
    #to help us keep track of the unique repeat units
    unit_list = []
    for uid, sequence_repeat_info in file_repeat_info.items():
        for repeat_unit, times in sequence_repeat_info.items():
            if repeat_unit in unit_list:
                file_unit_times_dict[repeat_unit] = file_unit_times_dict[repeat_unit] + times
            else:
                file_unit_times_dict[repeat_unit] = times
                unit_list.append(repeat_unit)

    return file_unit_times_dict

def max_occurance_n(fasta_dict, n):
    """This function returns the unit(s) with the maximum occurance with the number of times it occured."""

    repeat_units_dict = get_file_repeat_units(fasta_dict, n)
    max_times = 0

    for repeat_unit, times in repeat_units_dict.items():
        if times > max_times:
            max_times = times

    #now we build the max_repeats list. A list is used here because there can be more than one repeat unit with the maximum occurance
    max_repeats_list = []
    #we assign the first position of the max_repeats_list with the max time of occurance 
    max_repeats_list.append(max_times)
    for repeat_unit, times in repeat_units_dict.items():
        if times == max_times:
            max_repeats_list.append(repeat_unit)

    return max_repeats_list
