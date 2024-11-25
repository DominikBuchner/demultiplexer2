# Demultiplexer2

![logo](https://github.com/user-attachments/assets/e9c034d1-be0f-4e06-a78d-95fcaf03e926)

## Introduction 

**Demultiplexer2** is a Python package designed to efficiently demultiplex paired-end Illumina sequencing reads by identifying and sorting inline tags, enabling streamlined downstream analysis.

## Installation

You can install **demultiplexer2** using pip:

```
pip install demultiplexer2
```
When 

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
