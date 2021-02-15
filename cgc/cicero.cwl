class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: an.thrasher/cicero-benchmark/cicero/8
baseCommand: []
inputs:
  - id: disable_optimize
    type: boolean?
    inputBinding:
      position: 100
      prefix: '-no-optimize'
    label: Disable optimization
    doc: >-
      Disables optimizations. This can increase sensitivity, but increases the
      computational requirements.
  - id: num_soft_clip
    type: int?
    inputBinding:
      position: 103
      prefix: '-s'
    label: Soft Clip
    doc: minimum number of soft clip support required
  - id: threshold_soft_clip
    type: int?
    inputBinding:
      position: 102
      prefix: '-t'
    label: Threshold Soft Clip
    doc: threshold for enabling increased soft clip cutoff
  - id: cluster_distance
    type: int?
    inputBinding:
      position: 104
      prefix: '-c'
    label: Cluster distance
    doc: clustering distance for grouping similar sites
  - id: junction_file
    type: File?
    inputBinding:
      position: 100
      prefix: '-j'
    label: Junction File
    doc: junctions file from RNApeg
  - id: num_cores
    type: int?
    inputBinding:
      position: 1
      prefix: '-n'
    label: ncores
    doc: is the number of cores to be run on (with GNU parallel)
  - id: bam_file
    type: File
    inputBinding:
      position: 2
      prefix: '-b'
      shellQuote: false
    label: BAM File
    doc: input bamfile mapped to human genome builds GRCh37-lite or GRCh38_no_alt.
    'sbg:fileTypes': BAM
    secondaryFiles:
      - .bai
  - id: genome
    type:
      type: enum
      symbols:
        - GRCh37-lite
        - GRCh38_no_alt
      name: genome
    inputBinding:
      position: 3
      prefix: '-g'
  - id: refdir
    type: Directory
    inputBinding:
      position: 4
      prefix: '-r'
  - id: disable_excluded
    type: boolean?
    inputBinding:
      position: 104
      prefix: '-d'
    doc: disable excluded regions file use
outputs:
  - id: '#output_file'
    type: File?
    outputBinding:
      glob: >-
        $(inputs.bam_file.nameroot)/CICERO_DATADIR/$(inputs.bam_file.nameroot)/final_fusions.txt
label: cicero
arguments:
  - position: 101
    prefix: '-o'
    valueFrom: (inputs.bam_file.nameroot)
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'ghcr.io/stjude/cicero:latest'
