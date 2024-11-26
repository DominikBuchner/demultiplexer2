# Demultiplexer2

![logo](https://github.com/user-attachments/assets/e9c034d1-be0f-4e06-a78d-95fcaf03e926)

## Introduction 

**Demultiplexer2** is a Python package designed to efficiently demultiplex paired-end Illumina sequencing reads by identifying and sorting inline tags, enabling streamlined downstream analysis.

## Installation

You can install **demultiplexer2** using pip:

```
pip install demultiplexer2
```

When updates are released, users will be notified and can upgrade to the latest version with:

```
pip install --upgrade demultiplexer2
```

##  Usage

A *primerset* is a configuration file that stores details about the primers and tags used for a specific dataset. It is automatically saved in the `/demultiplexer2/data` directory for future use.

### Step 1: Create a Primer Set

To create a new primerset, use the following command:

```
demultiplexer2 create_primerset --name NameOfPrimerset --n_primers NumberOfPrimers
```
* --name: Specifies the name of the primerset (e.g., fwh2F2-fwhR2n).
* --n_primers: Indicates the number of primers in your dataset.

The *primerset* is an Excel file that contains critical information organized into three sheets:

1. General Information: Stores details about the primers used for amplification.
2. Forward Tags: Contains the names and sequences of the tags associated with the forward primers.
3. Reverse Tags: Contains the names and sequences of the tags associated with the reverse primers.

Fill out the primerset file before continuing with the next step.

### Step 2: Create a Tagging scheme

To create a new tagging scheme, use the following command:

```
demultiplexer2 create_tagging_scheme --name NameOfTaggingScheme --data_dir InputDirectory --primerset_path PathToPrimerset
```
* --name: Specifies the name of the tagging scheme (e.g., MyFirstStudy).
* --data_dir: Specifies the directory with all files you want to demultiplex.
* --primerset_path: Specifies the path to the primerset you want to use to demultiplex this dataset.

The *tagging scheme* is an Excel file that links your input files to sample names after demultiplexing. It will be save to the current working directory.
The sample names have to be added in the tagging scheme. Fill out the tagging scheme before continuing with the next step.

### Step 3: Demultiplex

To run the demultiplexing algorithm use this command:

```
demultiplexer2 demultiplex --primerset_path PathToPrimerset --tagging_scheme_path PathToTaggingScheme --output_dir OutputDirectory
```
* --primerset_path: Specifies the path to the primerset you want to use to demultiplex this dataset.
* --tagging_scheme_path: Specifies the path to the tagging scheme you want to use to demultiplex this dataset.
* --output_dir: Specifies the output directory to write to.

Given this information, demultiplexer2 will demultiplex your input to your output directory and give some statistics about how many reads could be assigned to tags. Unmatched reads will be directly discarded. 

