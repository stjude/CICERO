# CICERO 0.2.0

<div align="center">

  [![Status](https://img.shields.io/badge/status-active-success.svg)](https://github.com/stjude/CICERO)
  [![GitHub issues](https://img.shields.io/github/issues/stjude/CICERO)](https://github.com/stjude/CICERO/issues)
  [![GitHub Pull Requests](https://img.shields.io/github/issues-pr/stjude/CICERO)](https://github.com/stjude/CICERO/pulls)

</div>

CICERO (Clipped-reads Extended for RNA Optimization) is an assembly-based algorithm to detect diverse classes
of driver gene fusions from RNA-seq.

## 📝 Table of Contents
- [Running CICERO](#running)
- [Dependencies](#dependencies)
- [Running with Docker](#docker)
- [Running with St. Jude Cloud](#cloud)
- [Generate Junctions](#junctions)
- [Reference Files](#reference)
- [Supported Genomes](#genomes)
- [Demo](#demo)
- [Citation](#citation)
- [License](#license)

## Running CICERO <a name="running"></a>

Invoke the CICERO wrapper as
```
Cicero.sh [-n ncores] -b bamfile -g genome -r refdir [-j junctions] [-o outdir] [-t threshold] [-s sc_cutoff] [-c sc_shift] [-p]
```

Where `ncores` is the number of cores to be run on (with [GNU parallel](https://www.gnu.org/software/parallel/)),
`bamfile` is the input bamfile, `genome` is GRCh37-lite, `refdir` is the reference file directory specific to
CICERO, and `junctions` is the (optional) junctions file output from RNApeg.

Once you have the output from CICERO, use the following [guide](https://www.stjude.cloud/docs/guides/genomics-platform/analyzing-data/rapid-rnaseq/) to interpret the results.

## Dependencies <a name="dependencies"></a>

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

## Running with Docker <a name="docker"></a>

CICERO can be run with Docker. To begin, build the Docker image using the Dockerfile in this repository. 

```
docker build -t stjude/cicero:0.2.0 .
```

Then invoke the CICERO wrapper using Docker.

```
docker run -v reference:/reference stjude/cicero:0.2.0 [-n cores] -b <bam file path> -g <genome, e.g. GRCh37-lite> -r /reference [-j junctions file] [-o output directory] [-p] [-s int] [-t int] [-c int]
```

## Running with St. Jude Cloud <a name="cloud"></a>

CICERO is integrated in the St. Jude Cloud Rapid RNA-Seq workflow. To run CICERO in St. Jude Cloud, access the tool through the [platform page](https://platform.stjude.cloud/workflows/rapid_rna-seq). Documentation
for running and interpreting results is available in the [user guide](https://stjudecloud.github.io/docs/guides/genomics-platform/analyzing-data/rapid-rnaseq/).

## Generate junctions file with RNApeg <a name="junctions"></a>

RNApeg is required to generate a junctions file for use by CICERO. You can get RNApeg from both Docker and Singularity.

Running RNApeg via Docker:
```
docker run mnedmonson/public:rnapeg RNApeg.sh -b bamfile -f fasta -r refflat [-rg refflat]
```
Running RNApeg via Singularity:
```
singularity run docker://mnedmonson/public:rnapeg RNApeg.sh -b bamfile -f fasta -r refflat [-rg refflat]
```

## Downloading reference files <a name="reference"></a>

Reference files are required to run CICERO. They can be found at the following location:
* https://doi.org/10.5281/zenodo.3817656

## Supported Genome Version <a name="genomes"></a>

CICERO currently only supports `GRCh37-lite`. We are working towards support for `GRCh38` in the future. 

## Demo <a name="demo"></a>

A demo of CICERO can be found at the following location:
* https://www.stjuderesearch.org/site/lab/zhang/cicero

## Citation <a name="citation"></a>

Tian, L., Li, Y., Edmonson, M.N. et al. CICERO: a versatile method for detecting complex and diverse driver fusions using cancer RNA sequencing data. Genome Biol 21, 126 (2020). https://doi.org/10.1186/s13059-020-02043-x

## License <a name="license"></a>
Copyright 2019 St. Jude Children's Research Hospital

Licensed under a modified version of the Apache License, Version 2.0
(the "License") for academic research use only; you may not use this
file except in compliance with the License. To inquire about commercial
use, please contact the St. Jude Office of Technology Licensing at
scott.elmer@stjude.org.
    
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
