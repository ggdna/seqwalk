def rc(seq):
    """
    reverse complement of DNA sequence

    Args:
        seq: string with letters in {A, C, G, T}

    Returns:
        string corresponding to reverse complement
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

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
        for i in range(len(seq)-k+1):
            if sum([(s == "C" or s == "G") for s in seq[i:i+k]]) == 0 :
                if seq[i+middle-1] == "A":
                    to_remove.append(seq)

    return [seq for seq in library if seq not in to_remove]

def rc_hash_filtering(library, k):
    """
    filter any library to be RC free, using simple hash approach
    could be slow for large libraries

    Args:
        library: list of sequences
        k: SSM k value

    Returns:
        list of strings : filtered_library
             list of sequences without reverse complementary k-mers
    """
    seen_kmers = set()
    to_remove = set()

    for seq in library:
        bad_seq = False
        for i in range(len(seq)-k+1):
            
            if rc(seq[i:i+k]) in seq:
                bad_seq = True
            if rc(seq[i:i+k]) in seen_kmers:
                bad_seq = True
        if not bad_seq:
            for i in range(len(seq)-k+1):
                seen_kmers.add(seq[i:i+k])
        else:
            to_remove.add(seq)

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


