# centrolign

### About

`centrolign` is an multiple sequence alignment algorithm that is intended to be used for long tandem repeat nucleotide sequences. The original motivation for this algorithm was constructing pangenome graphs for the human centromeres (hence the name). However, it is also suitable for other tandem arrays of a similar scale (roughly 100 kbp-10 Mbp).

### Installation

`centrolign` supports building on macOS and Linux operating systems. The developers regularly build on macOS Monterey 12.3 and Ubuntu 22.04. Windows builds are not supported. 

##### Dependencies

* `cmake` ≥ 3.10
* `gcc` ≥ 4.8 or `clang` ≥ 3.3

##### Building and installing `centrolign`



```
git clone https://github.com/jeizenga/centrolign.git
cd centrolign
mkdir build
cd build
cmake .. # to choose install location, add '-DCMAKE_INSTALL_PREFIX=/path/to/install/'

# makes executable binary in the build directory
make -j 8

# OPTIONAL: installs executable, library, and headers 
make install
```

### Running `centrolign`

Two primary inputs are required to run `centrolign` as a command line utility:

1. The sequences to align in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).
2. A guide tree for the alignment in [Newick format](https://en.wikipedia.org/wiki/Newick_format).




### Limitations and known issues

