<p align="center">

  <h1 align="center">
    CICERO [UNMAINTAINED]
  </h1>

  <p align="center">
   <a href="https://github.com/stjude/CICERO" target="_blank">
     <img alt="Status" src="https://img.shields.io/badge/status-unmaintained-red">
   </a>
   <a href="https://github.com/stjude/CICERO/issues" target="_blank">
     <img alt="Github Issues"
          src="https://img.shields.io/github/issues/stjude/CICERO" />
   </a>
   <a href="https://github.com/stjude/CICERO/pulls" target="_blank">
     <img alt="Pull Requests"
          src="https://img.shields.io/github/issues-pr/stjude/CICERO" />
   </a>
   <a href="https://actions-badge.atrox.dev/stjude/CICERO/goto" target="_blank">
     <img alt="Actions: CI Status"
          src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fstjude%2FCICERO%2Fbadge&style=flat" />
   </a>
   <a href="https://github.com/stjude/CICERO/blob/master/LICENSE.md" target="_blank">
     <img alt="License: MIT"
          src="https://img.shields.io/badge/License-MIT-blue.svg" />
   </a>
   <img alt="Maintenance" src="https://img.shields.io/maintenance/no/2025">
  </p>


  <p align="center">
    CICERO (Clipped-reads Extended for RNA Optimization) is an assembly-based algorithm to detect diverse classes
    of driver gene fusions from RNA-seq.
   <br />
   <a href="https://stjude.github.io/CICERO/"><strong>Explore the docs »</strong></a>
   <br />
   <a href="https://stjuderesearch.org/site/lab/zhang/cicero"><strong>Work with demo data »</strong></a>
   <br />
   <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02043-x"><strong>Read the paper »</strong></a>
   <br />
   <br />
   <a href="https://github.com/stjude/cicero/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
   <a href="https://github.com/stjude/cicero/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
   <br />
    ⭐ Consider starring the repo! ⭐
   <br />
  </p>
</p>

---

  <p align="center">
  To discover driver fusions beyond canonical exon-to-exon chimeric transcripts, we develop CICERO, a local assembly-based algorithm that integrates RNA-seq read support with extensive annotation for candidate ranking. CICERO outperforms commonly used methods, achieving a 95% detection rate for 184 independently validated driver fusions including internal tandem duplications and other non-canonical events in 170 pediatric cancer transcriptomes.
   <img alt="Overview of CICERO algorithm which consists of fusion detection through analysis of candidate SV breakpoints and splice junction, fusion annotation, and ranking." src="CICERO.png"/>
  </p>


<br />

## 📝 Table of Contents

- [Unmaintained](#unmaintained)
- [Running CICERO](#running)
- [Dependencies](#dependencies)
- [Running with Docker](#docker)
- [Running with St. Jude Cloud](#cloud)
- [Generate Junctions](#junctions)
- [Reference Files](#reference)
- [Supported Genomes](#genomes)
- [Demo](#demo)
- [Output Fields](#outputfields)
- [Citation](#citation)
- [License](#license)

## Unmaintained <a name="unmaintained"></a>

As of 2025, CICERO is no longer receiving regular maintenance or updates. It should generally remain usable and suitable for use. Any questions can be addressed to the senior author of the CICERO [manuscript](#citation).

## Running CICERO <a name="running"></a>

Add the `src/scripts` directory to your system `PATH` variable. Add the `src/perllib` and `dependencies/lib/perl` directories to your system `PERL5LIB` variable.

Then invoke the CICERO wrapper as
```
Cicero.sh [-h] [-n ncores] -b bamfile -g genome -r refdir [-j junctions] [-o outdir] [-t threshold] [-s sc_cutoff] [-c sc_shift] [-p] [-d]

-p - optimize CICERO, sets sc_cutoff=3 and sc_shift=10 [default true]
-s <num> - minimum number of soft clip support required [default=2]
-t <num> - threshold for enabling increased soft clip cutoff [default=200000]
-c <num> - clustering distance for grouping similar sites [default=3]
-j <file> - junctions file from RNApeg
-n <num> - number of cores to utilize with GNU parallel
-d - disable excluded regions file use
```

- `ncores` is the number of cores to be run on (with [GNU parallel](https://www.gnu.org/software/parallel/)).
- `bamfile` is the input bamfile mapped to human genome builds GRCh37-lite or GRCh38_no_alt. [Contact us](mailto:liqing.tian@stjude.org) if your bam is based on other reference version.
- `genome` is either `GRCh37-lite` or `GRCh38_no_alt`. CICERO only support the two human reference genome versions.
- `refdir` is the reference file directory specific to CICERO. Download [Reference Files](#reference) below. e.g. `-r /home/user/software/CICERO/reference_hg38/` or `-r /home/user/software/CICERO/reference_hg19/`
- `outdir` is user defined output file folder.
- `junctions` is the junctions file output from RNApeg. See [Generate Junctions](#junctions) below. CICERO can detect fusion by analysis of splice junction reads. If this option is omitted, fusions generated by small deletions may be missed as these events may lack the soft-clipped reads.
- `threshold` CICERO first detects all soft-clipped positions supported by >=2 reads from bam file. For sample with <=threshold (default 200,000) soft-clipped positions, CICERO will detect fusions based on these soft-clipped positons; otherwise, to speed-up CICERO running, CICERO will detect fusions based on soft-clipped positions supported by >=sc_cutoff (3, default for optimize mode, see below) reads. For sample with lots of soft-clipped positions, a smaller threshold will speed-up CICERO running, however, some fusion events (i.e. only supported by 2 reads) may be missed.
- `sc_cutoff` controls the number of soft clip reads required to support a putative site. The default is 2, but for samples with large numbers of soft clip reads, it may be desirable to require additional support to reduce the computational time required.
- `sc_shift` sets the threshold for considering events to be the same site.
- `optimize` defaults to ON. This sets `sc_cutoff` to 3 for samples where the number of soft clip sites exceeds 200,000. It also sets `sc_shift` to 10 which sets the distance to consider events the same.
- `-no-optimize` turns optimizations off. This can increase sensitivity, but increases the computational requirements.

The final CICERO fusion result file will be located at `<outdir>/CICERO_DATADIR/<sample name>/final_fusions.txt`. Use the following [guide](https://university.stjude.cloud/docs/genomics-platform/workflow-guides/rapid-rnaseq/#interpreting-results) to interpret the results.

To visualize CICERO fusion output you can load the final fusion output file at https://proteinpaint.stjude.org/FusionEditor/.

## Dependencies <a name="dependencies"></a>

- [GNU parallel](https://www.gnu.org/software/parallel/)
- [Samtools 1.3.1](http://www.htslib.org/doc/samtools-1.3.1.html)
- [Cap3](https://www.ncbi.nlm.nih.gov/pubmed/10508846)
- [Blat](https://genome.ucsc.edu/goldenpath/help/blatSpec.html)
- Java 1.8.0
- Perl 5.10.1 with libraries:
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

CICERO can be run with Docker. Pre-built Docker images are provided for each release in GitHub Packages.

Invoke the CICERO wrapper using the Docker image available in GitHub Packages. You will likely need to add an additional bind mount for the output and input (BAM + junctions) files. Note the following command pulls the `latest` tag for the Docker image. For reproducible results, it is advisable to specify the exact version to run.

```bash
docker run -v <path to reference directory>:/reference ghcr.io/stjude/cicero:latest [-n cores] -b <bam file path> -g <genome> -r /reference -o <output directory> [-j junctions file] [-p] [-s int] [-t int] [-c int]
```

See [Running CICERO](#running) for details of the parameters.

## Running with St. Jude Cloud <a name="cloud"></a>

CICERO is integrated in the St. Jude Cloud Rapid RNA-Seq workflow. To run CICERO in St. Jude Cloud, access the tool through the [platform page](https://platform.stjude.cloud/workflows/rapid_rna-seq). Documentation
for running and interpreting results is available in the [user guide](https://university.stjude.cloud/docs/genomics-platform/workflow-guides/rapid-rnaseq/).

## Generate junctions file with RNApeg <a name="junctions"></a>

RNApeg is required to generate a junctions file for use by CICERO. You can get RNApeg from both Docker and Singularity. Once RNApeg is complete, the `*.junctions.tab.shifted.tab` file can be provided to CICERO using the `-j` argument. 

### Running RNApeg via Docker:
RNApeg is authored by Michael N. Edmonson ([@mnedmonson](https://github.com/mnedmonson)).

### RNApeg overview

This software analyzes nextgen RNA sequencing data which has been mapped to whole-genome coordinates, identifying evidence of both known and novel splicing events from the resulting alignments. The raw junction sites in the mapped BAMs undergo postprocessing to correct various issues related to mapping ambiguity. The result is a more compact and consistent set of junction calls, simplifying downstream quantification, analysis, and comparison.

### RNApeg key features

#### Raw junction extraction

First, the BAM read mappings are analyzed to identify putative junction sites. This produces a list of junction sites along with counts of supporting reads and several associated quality metrics. While reflective of the BAM data, this output typically requires refinement by the following steps.

#### Correction vs. reference junctions

Novel junctions are compared with reference exon junction boundaries and evaluated for mapping ambiguity which can justify adjusting the sites to match. Even small ambiguities such as the presence of the same nucleotide on either side of a junction can be enough to nudge a prediction that would otherwise perfectly match a reference isoform out of place.

#### Self-correction for novel junctions

Mapping ambiguity is next evaluated within the novel junctions themselves. Ambiguous junctions are combined where possible, merging their counts of supporting reads and related annotations. This reduces the population of novel junctions while simultaneously improving the evidence for those remaining. Combining evidence for poorly-covered sites also improves the chances of these sites passing the default minimum level of 3 supporting reads required for reporting junctions in the final output.

#### Correction vs. novel skips of known exons

Additional correction of novel junctions is also performed to identify previously unknown skips of an exon (or exons) within known reference isoforms. Special handling is required in these cases because while the corrected boundaries are known, the events themselves are novel.

#### Edge correction of novel junctions vs. reference exons

Junctions may also be shifted in cases of ambiguity involving a single edge (i.e. junction start or end). While not doubly-anchorable as with known reference junctions or novel skips of known exons, this adjustment can standardize evidence e.g. for novel exons.

#### Junction calling

Novel junctions are subjected to additional scrutiny before being reported:

- must be supported by a minimum of 3 reads
- at least one read must pass minimum flanking sequence requirements, to avoid false positives near read ends due to insufficient anchoring
- the junction must be either observed bidirectionally, or be supported by very clean alignments (either perfect or with very few high-quality mismatches, insertions, deletions, or soft clips)

While these requirements are minimal, they substantially reduce background noise.

#### Cross-sample correction

This step pools results for a set of samples and does additional standardization of novel exons based on the combined set. Mostly this has the effect of standardizing ambiguous novel junction sites across samples, but it can occasionally result in combinations of sites as well.

#### Output

The primary output files are tab-delimited text.

Output is also written in UCSC .bed format, which can be used to visualize the junctions and supporting read counts within the UCSC genome browser.

### Running RNApeg via Docker

```bash
docker run -v <outdir>:/results ghcr.io/stjude/rnapeg:latest -b bamfile -f fasta -r refflat
```

- `fasta` reference genome; i.e. "Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa" or "Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa" from [Reference Files](#reference).
- `refflat` i.e. "Homo_sapiens/GRCh38_no_alt/mRNA/RefSeq/refFlat.txt" or "Homo_sapiens/GRCh37-lite/mRNA/Combined/all_refFlats.txt" from [Reference Files](#reference).

### Running RNApeg via Singularity:

```bash
singularity run --containall --bind <outdir>:/results docker://ghcr.io/stjude/rnapeg:latest -b bamfile -f fasta -r refflat -O /results
```

You will also need to add `--bind` arguments to mount the file paths for `bamfile`, `fasta`, and `refflat` into the container. 

## Downloading reference files <a name="reference"></a>

Reference files are required to run CICERO. They can be found at the following location:

- GRCh37-lite: https://doi.org/10.5281/zenodo.3817656
- GRCh38_no_alt: https://doi.org/10.5281/zenodo.3894739

## Supported Genome Version <a name="genomes"></a>

CICERO currently supports `GRCh37-lite` and `GRCh38_no_alt`. 

## Demo <a name="demo"></a>

A demo of CICERO can be found at the following location:
* https://www.stjuderesearch.org/site/lab/zhang/cicero

## Output Fields <a name="outputfields"></a>
| Field    | Description |
|----------|-------------|
| sample   | Sample ID   |
| geneA / geneB | gene at breakpoint A / B |
| chrA / chrB | chromosome at breakpoint A / B |
| posA / posB | coordinate at breakpoint A / B |
| ortA / ortB | Mapping strand of assembled contig at breakpoint A / B |
| featureA / featureB | 5utr / 3utr / coding / intron / intergenic at breakpoint A / B |
| sv_ort | Whether the mapping orientation of assembled contig has confident biological meaning; if confident, then '>', else '?' (e.g. the contig mapping is from sense strand of gene A to antisense strand of gene B). |
| readsA / readsB | number of junction reads that support the fusion at breakpoint A / B |
| matchA / matchB | contig matched length at breakpoint A / B region |
| repeatA / repeatB | repeat score (0~1) at breakpoint A / B region, the higher the more repetitive |
| coverageA / coverageB | coverage of junction reads that support the fusion at breakpoint A / B (add the sequence length that can be mapped to the assembled contig for each junction read) |
| ratioA / ratioB | MAF of soft-clipped reads at breakpoint A / B (calculate the MAF for plus mapped reads and minus mapped reads, respectively; use the maximum MAF). |
| qposA / qposB | breakpoint position in the contig that belongs to A / B part |
| total_readsA / total_readsB | total reads number at the breakpoint at breakpoint A / B |
| contig | Assembled contig sequence that support the fusion |
| type | CTX (interchromosomal translocation) / Internal_dup / ITX (inversion) / DEL (deletion) / INS (insertion) / read_through |
| score | Fusion score, the higher the better |
| rating | HQ (known fusions) / RT (read_through) / LQ (others) |
| medal | Estimated pathogenicity assessment using St. Jude Medal Ceremony. Value: 0/1/2/3/4, the bigger the better |
| functional effect | ITD (Internal_dup) / Fusion / upTSS / NLoss / CLoss / other |
| frame | 0 (event is not in frame) / 1 (event is in-frame) / 2 (geneB portion contains canonical coding start site (i.e. the entire CDS for geneB)) / 3 (possible 5' UTR fusion in geneB) |

## Citation <a name="citation"></a>

Tian, L., Li, Y., Edmonson, M.N. et al. CICERO: a versatile method for detecting complex and diverse driver fusions using cancer RNA sequencing data. Genome Biol 21, 126 (2020). https://doi.org/10.1186/s13059-020-02043-x

## License <a name="license"></a>
Copyright 2020 St. Jude Children's Research Hospital

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
