def repeat_mapper(sequence, repeat_length):
    """This function goes over a sequence and finds all the repeats of length x and their occurances."""

    #we slide a window of 'length' nucleotides from position 1 to 1+length and then 2 to 2+length and so on, to get sequences of 'length' nucleotides.
    #This window will have to stop at end minus 
    i = 0
    while (i < = len(sequence) - repeat_length):
       repeat_candidate = sequence[i:(i + repeat_length)]
       j = 0
       #now we go over the sequence and find substrings matching the repeat sequence
       while (j <= (len(sequence) - repeat_length)):
           test_window = sequence[j:(j + repeat_length)]
           if test_window == repeat_candidate:
               new_match_pos = j
               repeat_enquirer = j + 1
               prev_match_pos = new_match_pos
               #for the occurance to be considered a repeat, the distance between the position of the present (new occurance) and that of the future occurance should not be greater than the length of the repeat 
               while((new_match_pos - (prev_match_pos + repeat_length)) <=0):
                   j = j + 1
                   test_window = sequence[j:(j + repeat_length)]
                   if test_window == repeat_candidate:
                       prev_match_pos = new_match_pos
                       new_match_pos = j 
                       repeat_counter = repeat_counter + 1




           j = j + 1

        i = i + 1
