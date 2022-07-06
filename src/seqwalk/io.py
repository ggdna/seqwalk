from importlib import resources

def load_library(identifier):
    """
    load a library of prebuilt sequences

    Args:
        identifier: string identifier of prebuilt library. listed on usage page

    Returns:
        list of strings : seqs
            library of orthogonal sequences

    """

    with resources.path("seqwalk.prebuilt_libs", identifier+".txt") as f:
        seqs = [s.strip() for s in open(f, "r").readlines()]
    return seqs

def write_library(seqs, filename):
    """
    writes a list of sequences to file. contains no information beyond sequence

    Args:
        seqs: list of strings
        filename: string corresponding to filename to save to
    Returns:
        None
    """
    f = open(filename, "w+")
    f.writelines([s + "\n" for s in seqs])
    f.close()
    print("File written!")
