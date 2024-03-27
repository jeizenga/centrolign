# centrolign

### About

`centrolign` is an multiple sequence alignment algorithm that is intended to be used for long tandem repeat nucleotide sequences. The original motivation for this algorithm was constructing pangenome graphs for the human centromeres (hence the name). However, it is also suitable for other tandem arrays of a similar scale (roughly 100 kbp-10 Mbp). `centrolign` has been successfully used to align >50 such sequences on a high memory compute server.

### Installation

`centrolign` supports macOS and Linux operating systems. The developers regularly build on macOS Monterey 12.3 and Ubuntu 22.04. Windows builds are not supported. 

#### Dependencies

* `cmake` ≥ 3.10
* `gcc` ≥ 4.8 or `clang` ≥ 3.3

#### Building and installing `centrolign`

```
git clone https://github.com/jeizenga/centrolign.git
cd centrolign
mkdir build
cd build
cmake .. # to choose install location, add '-DCMAKE_INSTALL_PREFIX=/path/to/install/'

# make executable binary in the build directory
make -j 8

# OPTIONAL: install the executable, library, and headers 
make install
```

### Running `centrolign`

Two primary inputs are required to run `centrolign` as a command line utility:

1. The sequences to align in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
2. A guide tree for progressive alignment in [Newick format](https://en.wikipedia.org/wiki/Newick_format).

**Important:** The sequence names from the FASTA file (i.e. following the `>` up to the first whitespace) must match the sample names in the Newick file identically.

The output of `centrolign` is a sequence graph that represents the multiple sequence alignment in [GFA format v1.0](https://gfa-spec.github.io/GFA-spec/GFA1.html). The syntax is:

```
centrolign -T guide_tree.nwk sequences.fasta > msa.gfa
```

Some notes:

* The guide tree (`-T`) is not strictly necessary, although highly recommended. If it is not provided, the sequences will be aligned in the order they are provided.
* The alignments of each progressive alignment subproblem can optionally be saved as GFA files by providing a prefix (`-S`).
* If `centrolign` was run with the `-S` parameter, a failed run can be restarted mid-execution using `-R`. The `-S` parameter must be the same in both runs.

#### Using `centrolign` for pairwise alignment

If a FASTA containing only two sequences is provided as input, `centrolign` outputs a [CIGAR string](https://en.wikipedia.org/wiki/Sequence_alignment#Representations) instead of a GFA file. Also, the guide tree (`-T`) is unnecessary for pairwise alignment.

```
centrolign two_sequences.fasta > cigar.txt
```

#### Using `centrolign` as a library

While it is primarily intended as a command line utility, the build process for `centrolign` also creates (and optionally installs) a library. To incorporate this library into another project, it will be necessary to replicate the basic I/O and parameter setting code from `centrolign`'s [main function](https://github.com/jeizenga/centrolign/blob/main/src/main.cpp).

### Limitations and known issues

* `centrolign` performs only global, co-linear alignment. Accordingly, the sequence graph outputs are all acyclic and they lack inversions. If cyclic or inverting motifs are necessary to align your sequences, you will need to build additional layers around it to generate global, acyclic alignment problems.


### Citation and credit

The design of `centrolign` was been substantially influenced by the pairwise alignment algorithm [UniAligner](https://github.com/seryrzu/unialigner).

There is currently no preprint or publication associated with the `centrolign` project. For the time being, cite this GitHub repository.