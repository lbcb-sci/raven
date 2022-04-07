# Raven

[![Latest GitHub release](https://img.shields.io/github/release/lbcb-sci/raven.svg)](https://github.com/lbcb-sci/raven/releases/latest)
![Build status for gcc/clang](https://github.com/lbcb-sci/raven/actions/workflows/raven.yml/badge.svg)
[![Published in Nature Computational Science](https://img.shields.io/badge/published%20in-Nature%20Computational%20Science-blue)](https://www.nature.com/articles/s43588-021-00073-4)

Raven is a de novo genome assembler for long uncorrected reads.

## Usage
To build raven executable run the following commands:

```bash
git clone https://github.com/lbcb-sci/raven && cd raven
cmake -S ./ -B./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

For faster build times optionally use [ninja](https://ninja-build.org/) and enable threading in cmake.
Eg.
```bash
git clone https://github.com/lbcb-sci/raven && cd raven
cmake -S ./ -B ./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release -G Ninja
cmake --build build -j 4
```

To install the raven executable after build run:

```bash
cmake --install ./build
```

To install python bindings run the following:
```bash
pip install git+git://github.com/lbcb-sci/raven.git@master
```

Python example can be found at `PythonLib/example.py`

```bash
usage: raven [options ...] <sequences> [<sequences> ...]

  # default output is to stdout in FASTA format
  <sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -k, --kmer-len <int>
      default: 15
      length of minimizers used to find overlaps
    -w, --window-len <int>
      default: 5
      length of sliding window from which minimizers are sampled
    -f, --frequency <double>
      default: 0.001
      threshold for ignoring most frequent minimizers
    -i, --identity <double>
      default: 0
      threshold for overlap between two reads in order to construct an edge between them
      if set to zero, this functionality is disabled
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

To use raven library component in your project, add the following to your cmake file:
```cmake
include(FetchContent)

FetchContent_Declare(
        raven
        GIT_REPOSITORY https://github.com/lbcb-sci/raven
        GIT_TAG v1.8.1)

FetchContent_GetProperties(raven)
if (NOT raven_POPULATED)
    FetchContent_Populate(raven)
    add_subdirectory(
            ${raven_SOURCE_DIR}
            ${raven_BINARY_DIR}
            EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(<YourTarget> <PRIVATE|PUBLIC|INTERFACE> raven)
```

#### Build options
- `RAVEN_BUILD_TESTS`: build unit tests
- `RAVEN_BUILD_PYTHON`: builds python module
- `RAVEN_BUILD_SHARED_LIBS: build raven lib and it's dependencies as shared libraries`
- `RAVEN_BUILD_EXE`: build raven executable
- `racon_enable_cuda`: build with NVIDIA CUDA support

#### Dependencies
- gcc 7.5+ | clang 8.0+
- cmake 3.11+
- zlib 1.2.8+

###### Hidden
- [pybind11](git@github.com:pybind/pybind11.git)
- lbcb-sci/racon/tree/library 3.0.2
- rvaser/bioparser 3.0.13
- (raven_test) google/googletest 1.10.0

### Other options

**NOTE**: not updated for 1.8 release

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
