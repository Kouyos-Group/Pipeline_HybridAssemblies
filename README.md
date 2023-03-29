# Genome Analysis pipeline using long and short reads (hybrid assemblies)

by Judith Bergadà Pijuan

This pipeline is aimed to perform the hybrid assemblies of bacterial genomes
using long reads polished with short paired-end reads.
Later, it performs the genomes annotation, compares the genome content of the
given DNA sequences, and carries out a gene variant calling analysis
in order to detect the SNPs across sequences. Finally, it identifies the
differentially located insertion sequences among strains.
Given multiple long reads and short paired-end sequencing reads (FASTQ files),
it provides a table file showing the genome content comparison, and (multiple)
tables showing the SNPs detected across strains. It also provides (multiple)
tables with the identification of the differentially located insertion
sequences among strains.
Outputs have the same format as given by software Roary, Snippy and ISCompare.
The pipeline also provides the hybrid assembly of the sequencing reads
and their annotation with prokka.

## Installation

To use this pipeline, you need to install the following dependencies:
- Filtlong
- Flye
- bwa
- Polypolish
- Prokka
- Roary
- Snippy
- ISCompare

You also need to have the following Python modules:
- Numpy
- Pandas
- Biopython
- Mechanize
- DNA_features_viewer
- lxml

Later, you need to download the tool:
```bash
cd $HOME
git clone https://github.com/judithbergada/Pipeline_HybridAssemblies
```

## Usage

The pipeline expects you to have the following folders:
- FASTQ folder with long reads: this is a folder containing only
your long sequencing reads (FASTQ files).
You must have all your long-reads FASTQ files here, and nothing else.
- FASTQ folder with short reads: this is a folder containing only
your short paired-end sequencing reads (FASTQ files).
You must have all your short-reads FASTQ files here, and nothing else.
It is important that the order of the long and the short reads belonging to
the same strain is the same within each folder. For this reason, we highly
recommend to always use the same prefix for each strain.


To get information about the usage, please try:

```bash
./hybridassembly.sh -h
```

The MLST tool can be used with these parameters:

```
Usage: hybridassembly.sh  [-h or --help]
                          [-l or --long_fastqfolder]
                          [-s or --short_fastqfolder]
                          [-o or --outname]
                          [-t or --threads]

Optional arguments:
    -h, --help:
                Show this help message and exit.
    -o, --outname:
                Name of your analysis.
                It will be used to name the output files.
                Default: mymlst.
    -t, --threads:
                Number of threads that will be used.
                It must be an integer.
                Default: 8.
Required arguments:
    -l, --long_fastqfolder:
                Path to the folder that contains ALL your Long-Reads.
                Only FASTQ files should be placed in it.
                You need long reads.
    -s, --short_fastqfolder:
                Path to the folder that contains ALL your Short-Reads.
                Only FASTQ files should be placed in it.
                You need forward and reverse paired-end reads.
                The order must be the same as for the long reads.
```

Enjoy using the tool!
