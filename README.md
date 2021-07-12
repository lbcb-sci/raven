# Raven

[![Latest GitHub release](https://img.shields.io/github/release/lbcb-sci/raven.svg)](https://github.com/lbcb-sci/raven/releases/latest)
![Build status for gcc/clang](https://github.com/lbcb-sci/raven/actions/workflows/raven.yml/badge.svg)
[![Published in Nature Computational Science](https://img.shields.io/badge/published%20in-Nature%20Computational%20Science-blue)](https://www.nature.com/articles/s43588-021-00073-4)

Raven is a de novo genome assembler for long uncorrected reads.

## Usage
To build raven run the following commands (< 30s):

```bash
git clone https://github.com/lbcb-sci/raven && cd raven && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

which will create raven executable and unit tests (running `make install` will install the executable to your system). Running the executable will display the following usage:

```bash
usage: raven [options ...] <sequences> [<sequences> ...]

  # default output is to stdout in FASTA format
  <sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    --weaken
      use larger (k, w) when assembling highly accurate sequences
    -p, --polishing-rounds <int>
      default: 2
      number of times racon is invoked
    -m, --match <int>
      default: 3
      score for matching bases
    -n, --mismatch <int>
      default: -5
      score for mismatching bases
    -g, --gap <int>
      default: -4
      gap penalty (must be negative)
    --graphical-fragment-assembly <string>
      prints the assembly graph in GFA format
    --resume
      resume previous run from last checkpoint
    --disable-checkpoints
      disable checkpoint file creation
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage

  only available when built with CUDA:
    -c, --cuda-poa-batches <int>
      default: 0
      number of batches for CUDA accelerated polishing
    -b, --cuda-banded-alignment
      use banding approximation for polishing on GPU
      (only applicable when -c is used)
    -a, --cuda-alignment-batches <int>
      default: 0
      number of batches for CUDA accelerated alignment
```

#### Build options
- `raven_build_tests`: build unit tests
- `racon_enable_cuda`: build with NVIDIA CUDA support

#### Dependencies
- gcc 4.8+ | clang 4.0+
- cmake 3.11+
- zlib 1.2.8+

###### Hidden
- lbcb-sci/racon/tree/library 3.0.2
- rvaser/bioparser 3.0.13
- (racon_test) google/googletest 1.10.0

### Other options

#### Brew
Install [Linuxbrew](https://docs.brew.sh/Homebrew-on-Linux) and run the following command:

```bash
brew install brewsci/bio/raven-assember
```

#### Conda
Install [conda](https://conda.io/en/latest/miniconda.html) and run the following command:
```bash
conda install -c bioconda raven-assembler
```

## Acknowledgment
This work has been supported in part by the Genome Institute of Singapore (A\*STAR), by the Croatian Science Foundation under projects Algorithms for genome sequence analysis (UIP-11-2013-7353) and Single genome and metagenome assembly (IP-2018-01-5886), and in part by the European Regional Development Fund under grant KK.01.1.1.01.0009 (DATACROSS).
