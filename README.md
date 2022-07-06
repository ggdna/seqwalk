# seqwalk

`seqwalk` is a package for designing orthogonal DNA sequence libraries.  If you want to design DNA barcodes (for sequencing, multiplexed imaging, molecular programming, etc.) `seqwalk` is for you! It can efficiently generate libraries of maximal size or maximal predicted orthogonality based on sequence symmetry. `seqwalk` additionally includes off-the-shelf orthogonal sequence libraries, as well as tools for analyzing orthogonal sequence libraries.

A code-free, interactive version of `seqwalk` can be found [here](https://colab.research.google.com/drive/1eVbcn_b5EE5FcL9NL5EyxeFAqNoImNSa?usp=sharing)

For more details, see our preprint (coming soon!).

## Installation

```bash
$ pip install seqwalk
```

## Usage

### Designing a set of barcodes with maximum orthogonality

If you want a certain number of barcodes with maximum orthogonality, you can use the `max_orthogonality` function from the `design` module. You must specify the length of desired sequences (L) and the number of desired sequences (N). Optionally, specify the prevention of reverse complementary sequences, GC content limits, allowable alphabet, and specific prevented patterns. By default, reverse complementary sequences are allowed, there are no GC content constraints, a 3 letter (A/C/T, no G) code is used and any 4N sequence is prevented.

For example, if you want 100 barcodes with length 25, with prevented reverse complements, and a 4 letter alphabet, and between 10 and 15 G/C bases, you can use the following code:

```python
from seqwalk import design

library = design.max_orthogonality(100, 25, alphabet="ACGT", RCfree=True, GClims=(10, 15))
```

### Designing a set of orthogonal barcodes with maximum size

If you have an orthogonality constraint in mind, you can use the `max_size` function from the `design` module to generate a maximally sized library. Orthogonality constraints must be sequence symmetry minimization k values. That is, the shortest k for which no substring of length k appears twice.

If you want sequences that satisfy SSM for k=12, and you want barcodes of length 25, without considering reverse complementarity, and using a 4 letter alphabet, with no GC constraints, you can use the following code:

```python
from seqwalk import design

library = design.max_size(25, 12, alphabet="ACGT")
```

### Importing "off-the-shelf" experimentally characterized libraries

The `io` module provides the ability to import libraries that have been previously experimentally characterized, using code of the following format.

```python
from seqwalk import io

PERprimers = io.load_library("kishi2018")
```

We provide the following libraries, accessible with the identifier tag.

| identifier | # of seqs | seq length | original use case | ref |
|------------|-----------|------------|-------------------|-----|
| `kishi2018` | 50 | 9nt | PER primers | [Kishi et al, 2018](https://www.nature.com/articles/nchem.2872) |

If you have an orthogonal library you would like to add, please submit a PR!

### Quality control using pairwise comparisons

Once you have a library in the form of a list of sequences, you can use the `analysis` module to perform additional quality control. For example, we provide a function to compute pairwise Hamming distances.

```python
from seqwalk import analysis

h_crosstalk = analysis.hamming_matrix(seqs)
```

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`seqwalk` is licensed under the terms of the MIT license.

## Credits

`seqwalk` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
