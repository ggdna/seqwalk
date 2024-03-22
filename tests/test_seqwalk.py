from seqwalk import design

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

def check_ssm(library, k, RCfree=False):
    """
    check a library to verify SSM is satisfied for length-k
    RCfree is boolean that is False if RC-tolerant
    returns True if SSM is satisfied
    """

    seen_kmers = set()

    for seq in library:
        for i in range(len(seq)-k+1):

            if seq[i:i+k] in seen_kmers:
                return False
            seen_kmers.add(seq[i:i+k])

            if RCfree:
                if rc(seq[i:i+k]) in seen_kmers:
                    return False
                seen_kmers.add(rc(seq[i:i+k]))
    return True

def check_pattern_free(library, pattern):

    for seq in library:
        if pattern in seq:
            return False
    return True

def check_GC(library, GClims):
    
    for seq in library:
        GC_count = sum([i in ['G', 'C'] for i in seq])
        if GC_count < GClims[0]:
            return False
        if GC_count > GClims[1]:
            return False
    return True


def test_3letter_RC_tolerant():
    ## test 3 letter no RCfree
    library = design.max_size(25, 6, alphabet="ACT")
    assert check_ssm(library, 6, RCfree=False), "SSM failed, 4 letter RC tolerant"
    library = design.max_size(25, 6, alphabet="TAC")
    assert check_ssm(library, 6, RCfree=False), "SSM failed, 4 letter RC tolerant"
    library = design.max_size(25, 6, alphabet="CAT")
    assert check_ssm(library, 6, RCfree=False), "SSM failed, 4 letter RC tolerant"

def test_4letter_RC_tolerant():
    ## test 4 letter no RCfree
    library = design.max_size(25, 6, alphabet="ACGT")
    assert check_ssm(library, 6, RCfree=False), "SSM failed, 4 letter RC tolerant"


def test_3letter_RC_free():
    ## test 3 letter filtered
    library = design.max_size(25, 4, alphabet="ACT", RCfree=True) 
    assert check_ssm(library, 4, RCfree=True), "SSM failed, 3 letter RC free"
    library = design.max_size(25, 5, alphabet="ACT", RCfree=True) 
    assert check_ssm(library, 5, RCfree=True), "SSM failed, 3 letter RC free"

def test_4letter_RC_free():
    library = design.max_size(25, 5, alphabet="ACTG", RCfree=True) 
    assert check_ssm(library, 5, RCfree=True), "SSM failed, 3 letter RC free"
    library = design.max_size(25, 4, alphabet="ACT", RCfree=True) 
    assert check_ssm(library, 4, RCfree=True), "SSM failed, 3 letter RC free"

def test_pattern_free():
    library = design.max_size(25, 5, alphabet="ACTG", RCfree=True)
    assert check_pattern_free(library, "AAAA"), "Failed pattern filtering"

def test_GC_lims():
    library = design.max_size(25, 5, alphabet="ACTG", GClims=(10, 20))
    assert check_GC(library, (10, 20)), "Failed GC filtering"

def test_max_orthogonality():
    library = design.max_orthogonality(100, 10, RCfree=True)
    assert check_ssm(library, 6), "Failed max orthogonality"
