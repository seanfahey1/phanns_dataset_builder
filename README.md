# phanns_dataset_builder

- [phanns\_dataset\_builder](#phanns_dataset_builder)
  - [Installation](#installation)
  - [How to Use](#how-to-use)
    - [Setting up a configuration file](#setting-up-a-configuration-file)
    - [Running the program](#running-the-program)
  - [How to generate an Entrez API key](#how-to-generate-an-entrez-api-key)
  - [Output](#output)



## Installation
`pip install phanns-dataset-builder`

## How to Use
### Setting up a configuration file
Create a config.toml file with the following format:
```
[positive_labels]
LABEL_1 = ["term1", "term2"]
LABEL_2 = ["term1", "term2", ... "termN"]
...
LABEL_N = ["term1"]

[query]
additional_query = [
    "AND refseq[filter]",
    "AND phage[Title]",
    "NOT hypothetical[Title]",
    "NOT putative[Title]",
    "NOT putitive[Title]",
    "NOT probable[Title]",
    "NOT possible[Title]",
    "NOT unknown[Title]",
    "AND 50:1000000[SLEN]"
    ]
```

The LABEL_1 key would build the following Entrez query:

```
(term1[Title] OR term2[Title]) AND refseq[filter] AND phage[Title]
NOT hypothetical[Title] NOT putative[Title] NOT putitive[Title] NOT probable[Title] NOT
possible[Title] NOT unknown[Title] AND 50:1000000[SLEN]
```

### Running the program
`dataset-builder -c [./relative/path/to/config.toml] -d [data directory] -l [logs directory] -e [your email here] -a [api key]`

## How to generate an Entrez API key
See the [NCBI walkthrough](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us)
for instructions on obtaining your own Entrez api key.

## Output
__Data Files__

The program will create one directory for each label in the configuration file inside the user specified data folder. Each record retrieved that matches the query will be stored as a `.fasta` formatted file using the accession number as the name.

__Log Files__

A single log file will be written in the user specified logging directory with the name `Entrez_info_[timestamp of start time].log`. The log file will contain the query executed, the number of expected records, and the status of each record retrieved.
