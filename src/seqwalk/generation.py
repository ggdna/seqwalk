from itertools import product
import random



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




def partition_path(seq, L, k):
    """
    partitions self-avoiding walk into appropriate length sequences

    Args:
        seq: self avoiding walk (string or list)
        L: length of desired sequences
        k: SSM k-value

    Returns:
        seqs: list of strings, each a length L seq
    """
    seqs = []
    for i in range(0, len(seq)-L+1, L-k+1):
        seqs.append(seq[i:i+L])
    return seqs


def is_necklace(seq):
    """
    computes if a sequence is a necklace (as defined in Wong 2017)

    Args:
        seq: typically list of ints
    """

    p = 1

    for i in range(1, len(seq)):
        if seq[i-p] < seq[i]:
            p = i + 1
        elif seq[i-p] > seq[i]:
            return False

    if (len(seq) % (p)) == 0:
        return True
    return False


def f(seq, k):
    """
    helper function defined in Wong 2017
    * note confusing notation switch

    Args:
       seq: list of ints
       k: alphabet length

    Returns:
       int corresponding to next element of seq
    """

    p = 1

    if seq[0] == k:
        if sum(seq[1:]) == len(seq) - 1:
            return 1
        for i in range(0, k-1):
            if not is_necklace(seq[1:] + [k-i]):
                return k-i
        return 1

    elif is_necklace(seq[1:] + [(seq[0] % k) + 1]):
        return (seq[0] % k) + 1

    else:
        return seq[0]

def simple_shift(n, k):
    """
    simple shift rule from Wong 2017
    * note confusing notation switch

    Args:
        n: SSM k value
        k: alphabet size

    Returns:
        list of integers corresponding to H. path
    """

    seq = [1]*n

    while True:
        seq.append(f(seq[-n:], k))
        if seq[-n:] == seq[:n]:
            return seq




def out_edges(v, alphabet):
    """
    lists outedges for a node in a kmer graph

    Args:
        v: string corresponding to node
        alphabet: string containing all valid letters

    Returns:
        list of outedges represented as strings
    """
    return [v + l for l in alphabet]


def adapted_hierholzer(k, alphabet):
    """
    finds RC-free path through 4-letter kmer graph using modified hierholzer
    (see Supplementary Note X)

    Args:
        k: SSM k (integer)
        alphabet: string containing all valid letters

    Returns:
        path: string corresponding to RC-free self-avoiding walk
    """

    assert (k % 2 == 1), "k must be odd for adapted Hierholzer"

    # dictionary to store visited nodes
    # list of all nodes
    marked = {"".join(l) : 0 for l in product(alphabet, repeat=k)}
    nodes = ["".join(l) for l in product(alphabet, repeat=(k-1))]

    # initialize stack with a random starting node
    v_stack = [random.choice(nodes)]
    path = ""

    while len(v_stack) != 0:
        v = v_stack.pop()
        unmarked = [s for s in out_edges(v, alphabet) if not marked[s]]

        if unmarked == []:
            if path == "":
                path = v
            else:
                path = v[0] + path
        else:
            n_edge = random.choice(unmarked)
            v_stack.append(v)
            v_stack.append(n_edge[1:])
            marked[n_edge] = 1
            marked[rc(n_edge)] = 1

    return path








