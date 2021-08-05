{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
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
            }
        },
        {
            "loadListing": "deep_listing",
            "id": "refdir",
            "type": "Directory",
            "inputBinding": {
                "prefix": "-r",
                "shellQuote": false,
                "position": 4
            }
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
            "dockerPull": "ghcr.io/stjude/cicero:latest"
        }
    ]
}
