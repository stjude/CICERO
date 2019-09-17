### Cicero 0.1.6

(Explain what this is and what paper to cite)

Invoke the Cicero wrapper as
```
Cicero.sh [-n ncores] -b bamfile -g genome -r refdir [-j junctions]
```

Where `ncores` is the number of cores to be run on (with [GNU parallel](https://www.gnu.org/software/parallel/),
`bamfile` is the input bamfile, `genome` is GRCh37-lite, `refdir` is the reference file directory specific to
Cicero, and `junctions` is the (optional) junctions file output from RNApeg.

### Dependencies

* [GNU parallel](https://www.gnu.org/software/parallel/)
* [Samtools 1.3.1](http://www.htslib.org/doc/samtools-1.3.1.html)
* [Cap3](https://www.ncbi.nlm.nih.gov/pubmed/10508846)
* [Blat](https://genome.ucsc.edu/goldenpath/help/blatSpec.html)
* Java 1.8.0
* Perl 5.10.1 with libraries:
    - base
    - Bio
    - Carp
    - Compress
    - Cwd
    - Data
    - DBI
    - diagnostics
    - Digest
    - English
    - enum
    - Exporter
    - File
    - FileHandle
    - List
    - POSIX
    - strict
    - Sys
    - Tree
    - warnings

### Downloading RNApeg

RNApeg is required to generate a junctions file for use by Cicero. You can get RNApeg from both docker and singularity.

Running RNApeg via Docker:
* docker run mnedmonson/public:rnapeg RNApeg.sh

Running RNApeg via Singularity:
* singularity run docker://mnedmonson/public:rnapeg
* singularity run docker://mnedmonson/public:rnapeg RNApeg.sh

### Downloading reference files

(Explain how to download the reference files here)
