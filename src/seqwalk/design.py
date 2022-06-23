from seqwalk.generation import *
from seqwalk.filtering import *
import math


def max_size(L, k, alphabet="ACT", RCfree=False, GClims=None,
             prevented_patterns=["AAAA", "CCCC", "GGGG", "TTTT"]):
    """
    design a max size library of length L sequences with SSM k

    Args:
        L: integer length of desired seqs
        k: SSM k value
        alphabet: string of allowable letters (default "ACT")
        RCfree: bool, True if orthogonality with RCs is required
        GClims: tuple of (GCmin, GCmax), allowable range of number of GC bases
        prevented_patterns: list of prevented patterns (default 4N)

    Returns:
        list of strings : seqs
            library of orthogonal sequences
    """

    assert (L > k), "L must greater than k"

    if RCfree and len(alphabet) == 4:
        seq = adapted_hierholzer(k, alphabet)
        seqs = partition_path(seq, L, k)
        if len(seqs) == 0:
            return seqs
        if GClims != None:
             GCmin, GCmax = GClims
             seqs = filter_gc(seqs, GCmin, GCmax)
        for pattern in prevented_patterns:
             seqs = filter_pattern(seqs, pattern)
        return seqs

    else:
         seq = "".join([alphabet[i-1] for i in simple_shift(k, len(alphabet))])
         seqs = partition_path(seq[:-1], L, k)
         if RCfree:
             seqs = filter_rc_3letter(seqs, k)
         if len(seqs) == 0:
             return seqs
         if GClims != None:
             GCmin, GCmax = GClims
             seqs = filter_gc(seqs, GCmin, GCmax)
         for pattern in prevented_patterns:
             seqs = filter_pattern(seqs, pattern)
         return seqs


def max_orthogonality(N, L, alphabet="ACT", RCfree=False, GClims=None,
                      prevented_patterns=["AAAA", "CCCC", "GGGG","TTTT"],
                      k_init=None):
    """
    design a maximally orthogonal library of N length L sequences

    Args:
        N: minimum number of sequences in library
        L: integer length of desired seqs
        alphabet: string of allowable letters (default "ACT")
        RCfree: bool, True if orthogonality with RCs is required
        GClims: tuple of (GCmin, GCmax), allowable range of number of GC bases
        prevented_patterns: list of prevented patterns (default 4N)
        k_init: initial guess for SSM k value
    
    Returns:
        list of strings : seqs
            library of orthogonal sequences

    """
    if k_init == None:
        k_init = max(int(math.log(N)/math.log(len(alphabet))), 2)
        if (k_init % 2 == 0) and RCfree:
            k_init += 1

    while True:

        library = max_size(L, k_init, alphabet, RCfree, GClims, prevented_patterns)

        if len(library)>N:
            print("Number of sequences: %d" % len(library))
            print("SSM k value: %d" % k_init)
            return library

        k_init += (1 + RCfree)
