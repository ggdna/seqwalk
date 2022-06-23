def load_library(identifier):
    f = open("prebuilt_libs/%s.txt" % identifier, "r")
    seqs = [s.strip() for s in f.readlines()]
    f.close()
    return seqs

def write_library(seqs, filename):
    f = open(filename, "w+")
    f.writelines(seqs)
    f.close()
    print("File written!")
