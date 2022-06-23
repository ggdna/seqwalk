def filter_rc_3letter(library, k):
    """
    filter library to be RC free
    (Supplementary note X)

    Args:
        library: list of sequences
        k: SSM k value

    Returns:
        list of strings : filtered_library
             list of sequences without reverse complementary k-mers
    """

    assert (k % 2 == 1), "SSM k must be odd for RC filtering"

    to_remove = []
    middle = int((k+1)/2)

    for seq in library:
        for i in range(len(seq)-k):
            if sum([(s == "C" or s == "G") for s in seq[i:i+k]]) == 0 :
                if seq[i+middle-1] == "A":
                    to_remove.append(seq)

    return [seq for seq in library if seq not in to_remove]



def filter_gc(library, gc_min, gc_max):
    """
    filters library for sequences that have desired GC content

    Args:
        library: list of sequences in string representation
        gc_min: minimum number of GC bases (int)
        gc_max: maximimum number of GC bases (int)

    Returns:
        list of strings : filtered_library
            list of sequences in string representation
    """

    assert (gc_min <= gc_max), "gc_min cannot be greater than gc_max"
    assert (gc_max <= len(library[0])),  "gc_max cannot be greater than seq length"

    filtered_library = []

    for seq in library:

        gc = sum([(s == "C" or s == "G") for s in seq])

        if gc >= gc_min:
            if gc <= gc_max:
                filtered_library.append(seq)

    return filtered_library



def filter_pattern(library, pattern):
    """
    filters library to remove specific patterns

    Args:
        library: list of sequences in string representation
        pattern: sequence pattern to be prevented

    Returns:
        list of strings : filtered_library
            list of sequences in string representation
    """

    filtered_library = []

    for seq in library:
        if pattern not in seq:
            filtered_library.append(seq)

    return filtered_library


