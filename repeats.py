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

print(get_repeat_frequency(test_sequence, "AGTAGGA"))
