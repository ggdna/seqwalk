import csv
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

def write_csv(seq, filename):
    """
    Writes a list of sequences to a .csv file. 
    Copy the sequence column into an IDT order template

    Args:
        seqs: list of strings
        filename: string corresponding to filename to save to
    Returns:
        None
    """
    with open (filename + '.csv', 'w') as f:
        write = csv.writer(f, delimiter='\n')
        write.writerow(["Sequence"])
        write.writerow(seq)

