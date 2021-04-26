{
    "class": "Workflow",
    "cwlVersion": "v1.1",
    "label": "RNApeg + CICERO",
    "inputs": [
        {
            "id": "in_bam",
            "sbg:fileTypes": "BAM",
            "type": "File",
            "label": "BAM File",
            "source": "in_bam",
        },
        {
            "id": "in_fasta",
            "sbg:fileTypes": "FA",
            "type": "File",
            "label": "FASTA File",
            "source": "in_fasta",
        },
        {
            "id": "in_refflat",
            "sbg:fileTypes": "TXT",
            "type": "File",
            "label": "RefFlaf Text File",
            "source": "in_refflat",
        },
        {
            "id": "in_ref",
            "type": "Directory",
            "label": "Reference Directory",
            "source": "in_refflat",
            "loadListing": "deep_listing"
        },
        {
            "id": "in_genome",
            "type": {
                "type": "enum",
                "symbols": [
                    "GRCh37-lite",
                    "GRCh38_no_alt"
                ],
                "name": "in_genome"
            },
            "label": "Genome",
            "default": "GRCh38_no_alt",
        },
        {
            "id": "num_cores",
            "type": "int?",
            "label": "ncores",
            "doc": "is the number of cores to be run on (with GNU parallel)",
            "sbg:exposed": true
        }
    ],
    "outputs": [
        {
            "id": "output_file",
            "outputSource": [
                "cicero/output_file"
            ],
            "type": "File",
        }
    ],
    "steps": [
        {
            "id": "cicero",
            "in": [
                {
                    "id": "disable_optimize",
                    "default": false
                },
                {
                    "id": "junction_file",
                    "source": "rnapeg/output_file"
                },
                {
                    "id": "num_cores",
                    "source": "num_cores"
                },
                {
                    "id": "bam_file",
                    "source": "in_bam"
                },
                {
                    "id": "genome",
                    "source": "in_genome"
                },
                {
                    "id": "refdir",
                    "source": "in_ref"
                }
            ],
            "out": [
                {
                    "id": "output_file"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.1",
                "id": "cicero",
                "baseCommand": [
                    "Cicero.sh"
                ],
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
                        "id": "output_file",
                        "type": "File",
                        "outputBinding": {
                            "glob": "$(inputs.bam_file.nameroot)/CICERO_DATADIR/$(inputs.bam_file.nameroot)/final_fusions.txt"
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
                ],
            },
            "label": "cicero",
        },
        {
            "id": "rnapeg",
            "in": [
                {
                    "id": "fasta",
                    "source": "in_fasta"
                },
                {
                    "id": "bam_file",
                    "source": "in_bam"
                },
                {
                    "id": "refFlat",
                    "source": "in_refflat"
                }
            ],
            "out": [
                {
                    "id": "output_file"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.1",
                "id": "rnapeg",
                "baseCommand": [],
                "inputs": [
                    {
                        "id": "fasta",
                        "type": "File",
                        "inputBinding": {
                            "prefix": "-f",
                            "shellQuote": false,
                            "position": 32
                        },
                        "label": "FASTA File",
                        "doc": "FASTA file",
                        "sbg:fileTypes": "FA",
                        "secondaryFiles": [
                            {
                                "pattern": ".fai"
                            }
                        ]
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
                        "id": "refFlat",
                        "type": "File",
                        "inputBinding": {
                            "prefix": "-r",
                            "shellQuote": false,
                            "position": 50
                        },
                        "label": "refFlat Text",
                        "doc": "Uncompressed refFlat.txt file from the UCSC genome annotation database",
                        "sbg:fileTypes": "TXT"
                    }
                ],
                "outputs": [
                    {
                        "id": "output_file",
                        "type": "File",
                        "outputBinding": {
                            "glob": "$(inputs.bam_file.nameroot).bam.junctions.tab.shifted.tab"
                        }
                    }
                ],
                "label": "rnapeg",
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 4000,
                        "coresMin": 0
                    },
                    {
                        "class": "InitialWorkDirRequirement",
                        "listing": [
                            "$(inputs.bam_file)"
                        ]
                    },
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "hints": [
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "ghcr.io/stjude/rnapeg:latest"
                    }
                ]
            },
            "label": "rnapeg"
        }
    ],
    "requirements": [
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ]
}
