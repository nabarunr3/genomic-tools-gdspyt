test_sequence="AGGCTTGACCTTGTGCGGCAGCCGGCAGCCGCCGGCCGGCAGGTTGCCGTACGGATCCTCGTCGGCCGCCGTCGCGGGGGCCGCCGGCATCACGGGGGCGGTCGACGCGACCAGCATCGGGGGGAACGGCAACGCATGCACCCGCTCGCCGGTCAGCACGAACGCACCGTTGGCCACCGCCGGCGCGATCGGCGGCACGCCGGCGTCGCTCAACCCGGTCGGCTGGGCGTCCGACGGCACGAAGAAGACATCCACCGGCGGCGCTTCCTGCATGCGTATCGGCGAATAGTCGGCGAAACCGGCGTTGCGGACCGCGCCATGGTCGACGTCGATCGCGAAGCCGGGCTTCGTCGTCGCGAGACCGAACAGCGCGCCGCCCTGGATCTGCGCTTGCGCGCCGGTCGGGTTGACGATGCGGCCCGCATACACGCCGGCCGTCACGCGATGCACGCGCGGTTGTTGCGCTTCGATCGACACTTCCGTCACGTACGCGACGACCGAGCCGGCCGTTTCGTGCATCGCGACGCCCCACGCGTGCCCGGCCGGCAGCGTGCGCGCGCCGTAGCCGGACTTGTCGACGGCCAGCGCGAGCGCCTGCCGATGCGCGGCGTGCTCGGGGCCGGCCAGCCGCGTCATCCGGTAGGCGACCGGATCCTGCCGCGCCGAGTGCGCGAGCTCGTCGACCAGCGTTTCCATCACGAACGCCGTATGCGAGTTGCCGCCCGAGCGCCACGTCTGGACCGGCACGTCGGCCTCGGTCTGATGAACCGATACCTGCATCGGGAAGCCGTACGGGCTGTTCGTCACGCCTTCGGTCAGGCTCGGATCGGTGCCGCGCTTGAGCATCGTCGTGCGCTCGAGCGGCGAGCCCTTCAGCACAGACTGGCCGACGACCACGTGCTGCCAGTCGCGCACGGCGCCGCTGCCGTCCACGCCGATGTCGACGCGATGCAGCACCATCGGGCGGTAATAGCCGCCGCGCAGATCGTCCTCGCGCGTCCAGATCGTCTTGACGGGGCCGAGATGGCCGGCCGCGAGGTACGCGGCGGACACGTGGGCGGCTTCGACCACGTAGTCCGACGTCGGCGTCGAGCGCCGGCCATAGTCGCCGCCCGAGGTCAGCGTGAAGATCTGGACTTTCTCCGGGGCGACGCCGAGCGCCTTCGCGACCGCCGCGCGGTCGGTCGTC"

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

def get_repeat_units(sequence, n):
    """This function gets a dictionary containing all subsequences of length n which occur repeatedly, that is, two times or more."""

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
