{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "baseCommand": [],
    "inputs": [
        {
            "id": "disable_optimize",
            "type": "boolean?",
            "inputBinding": {
                "prefix": "-no-optimize",
                "shellQuote": false,
                "position": 100
            },
            "label": "Disable optimization",
            "doc": "Disables optimizations. This can increase sensitivity, but increases the computational requirements."
        },
        {
            "id": "num_soft_clip",
            "type": "int?",
            "inputBinding": {
                "prefix": "-s",
                "shellQuote": false,
                "position": 103
            },
            "label": "Soft Clip",
            "doc": "minimum number of soft clip support required"
        },
        {
            "id": "threshold_soft_clip",
            "type": "int?",
            "inputBinding": {
                "prefix": "-t",
                "shellQuote": false,
                "position": 102
            },
            "label": "Threshold Soft Clip",
            "doc": "threshold for enabling increased soft clip cutoff"
        },
        {
            "id": "cluster_distance",
            "type": "int?",
            "inputBinding": {
                "prefix": "-c",
                "shellQuote": false,
                "position": 104
            },
            "label": "Cluster distance",
            "doc": "clustering distance for grouping similar sites"
        },
        {
            "id": "junction_file",
            "type": "File?",
            "inputBinding": {
                "prefix": "-j",
                "shellQuote": false,
                "position": 100
            },
            "label": "Junction File",
            "doc": "junctions file from RNApeg"
        },
        {
            "id": "num_cores",
            "type": "int?",
            "inputBinding": {
                "prefix": "-n",
                "shellQuote": false,
                "position": 1
            },
            "label": "ncores",
            "doc": "is the number of cores to be run on (with GNU parallel)"
        },
        {
            "id": "bam_file",
            "type": "File",
            "inputBinding": {
                "prefix": "-b",
                "shellQuote": false,
                "position": 2
            },
            "label": "BAM File",
            "doc": "input bamfile mapped to human genome builds GRCh37-lite or GRCh38_no_alt.",
            "sbg:fileTypes": "BAM",
            "secondaryFiles": [
                {
                    "pattern": ".bai?"
                },
                {
                    "pattern": "$(self.nameroot).bai?"
                }
            ]
        },
        {
            "id": "genome",
            "type": {
                "type": "enum",
                "symbols": [
                    "GRCh37-lite",
                    "GRCh38_no_alt"
                ],
                "name": "genome"
            },
            "inputBinding": {
                "prefix": "-g",
                "shellQuote": false,
                "position": 3
            },
            "doc": "Reference to which the input BAM is aligned. Must specify the matching reference data."
        },
        {
            "loadListing": "deep_listing",
            "id": "refdir",
            "type": "Directory",
            "inputBinding": {
                "prefix": "-r",
                "shellQuote": false,
                "position": 4
            },
            "doc": "Directory of reference files. Download and extract the reference archives available from https://doi.org/10.5281/zenodo.3817656 (GRCh37-lite) or https://doi.org/10.5281/zenodo.3894739 (GRCh38_no_alt)."
        },
        {
            "id": "disable_excluded",
            "type": "boolean?",
            "inputBinding": {
                "prefix": "-d",
                "shellQuote": false,
                "position": 104
            },
            "doc": "disable excluded regions file use"
        }
    ],
    "outputs": [
        {
            "id": "#output_file",
            "type": "File?",
            "outputBinding": {
                "glob": "$(inputs.bam_file.nameroot)/CICERO_DATADIR/$(inputs.bam_file.nameroot)/$(inputs.bam_file.nameroot).final_fusions.txt"
            }
        }
    ],
    "doc": "# Description\n\nTo discover driver fusions beyond canonical exon-to-exon chimeric transcripts, we developed CICERO, a local assembly-based algorithm that integrates RNA-seq read support with extensive annotation for candidate ranking. CICERO outperforms commonly used methods, achieving a 95% detection rate for 184 independently validated driver fusions including internal tandem duplications and other non-canonical events in 170 pediatric cancer transcriptomes.\n\nThe first step is to search the BAM file for soft clipped reads. The soft clipped reads are then filtered, assembled, and annotated to find candidate structural variations. The candidate structural variations are then filtered to produce a final call set of structural variants.\n\n\n## Inputs:\n* **Tumor RNA-Seq BAM** - Indexed tumor sample BAM\n* **Reference directory** - Directory of reference files for the genome of interest. These can be acquired from: https://doi.org/10.5281/zenodo.3817656 (GRCh37-lite) or https://doi.org/10.5281/zenodo.3894739 (GRCh38_no_alt).\n* **Junction file** - Junction annotations from RNApeg\n* **Genome** - Human genome to use: [GRCh37-lite, GRCh38_no_alt]\n\n\n## Outputs:\n* **Final Fusions** - CICERO called fusion events",
    "label": "cicero",
    "arguments": [
        {
            "prefix": "-o",
            "shellQuote": false,
            "position": 101,
            "valueFrom": "$(inputs.bam_file.nameroot)"
        },
        {
            "prefix": "-f",
            "shellQuote": false,
            "position": 102,
            "valueFrom": "$(inputs.bam_file.nameroot).final_fusions.txt"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": 24000,
            "coresMin": 2
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "hints": [
        {
            "class": "DockerRequirement",
            "dockerPull": "cgc-images.sbgenomics.com/stjude/cicero:latest"
        }
    ],
    "sbg:categories": [
        "RNA",
        "Variant Calling"
    ],
    "sbg:links": [
        {
            "id": "https://github.com/stjude/CICERO",
            "label": "Source Code"
        },
        {
            "id": "https://doi.org/10.1186/s13059-020-02043-x",
            "label": "Publication"
        }
    ]
}
