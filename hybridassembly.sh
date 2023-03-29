#!/bin/bash

set -e

###################
## Assign inputs ##
###################

# Define usage of the script
function print_usage {
  printf """
Usage: hybridassembly.sh  [-h or --help]
                          [-l or --long_fastqfolder]
                          [-s or --short_fastqfolder]
                          [-o or --outname]
                          [-t or --threads]
"""
}
# Describe usage of the tool and provide help
function print_help {
  print_usage
  printf """
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
"""
}

# Define inputs
for ARGS in "$@"; do
  shift
        case "$ARGS" in
                "--long_fastqfolder") set -- "$@" "-l" ;;
                "--short_fastqfolder") set -- "$@" "-s" ;;
                "--outname") set -- "$@" "-o" ;;
                "--threads") set -- "$@" "-t" ;;
                "--help") set -- "$@" "-h" ;;
                *) set - "$@" "$ARGS"
        esac
done

# Define defaults
outn="myLongShortAnalysis"; threads=8

# Define all parameters
while getopts 'l:s:o::t::h' flag; do
        case "${flag}" in
                l) long_fastqfolder=${OPTARG} ;;
                s) short_fastqfolder=${OPTARG} ;;
                o) outn=${OPTARG} ;;
                t) threads=${OPTARG} ;;
                h) print_help
                   exit 1;;
                *) print_usage
                    exit 1;;
        esac
done


##############################
## Identify Software Errors ##
##############################

printf "\nChecking if required software is installed...\n"
# Check installation of the required software.
# If something is missing, show how to install it.
if ! [ -x "$(command -v filtlong)" ]; then
  echo "Missing: Filtlong not found"
  echo "Information on the installation:"
  echo "https://github.com/rrwick/Filtlong"
  exit 127
fi
if ! [ -x "$(command -v flye)" ]; then
  echo "Missing: Flye not found"
  echo "Information on the installation:"
  echo "https://github.com/fenderglass/Flye"
  exit 127
fi
if ! [ -x "$(command -v bwa)" ]; then
  echo "Missing: bwa not found"
  echo "Information on the installation:"
  echo "https://github.com/lh3/bwa"
  exit 127
fi
if ! [ -x "$(command -v polypolish)" ]; then
  echo "Missing: Polypolish not found"
  echo "Information on the installation:"
  echo "https://github.com/rrwick/Polypolish/wiki"
  exit 127
fi
if ! [ -x "$(command -v prokka)" ]; then
  echo "Missing: Prokka not found"
  echo "Information on the installation:"
  echo "https://github.com/tseemann/prokka"
  exit 127
fi
if ! [ -x "$(command -v roary)" ]; then
  echo "Missing: Roary not found"
  echo "Information on the installation:"
  echo "https://github.com/sanger-pathogens/Roary"
  exit 127
fi
if ! [ -x "$(command -v snippy)" ]; then
  echo "Missing: Snippy not found"
  echo "Information on the installation:"
  echo "https://github.com/tseemann/snippy"
  exit 127
fi
if ! [ -x "$(command -v ISCompare.py)" ]; then
  echo "Missing: ISCompare not found"
  echo "Information on the installation:"
  echo "https://github.com/maurijlozano/ISCompare"
  exit 127
fi
if ! python -c 'import pkgutil; exit(not pkgutil.find_loader("Bio"))'; then
  echo "Missing: biopython package not found"
  echo "Try to install it this way: pip install biopython"
exit 127
fi
if ! python -c 'import pkgutil; exit(not pkgutil.find_loader("mechanize"))';
  then
  echo "Missing: mechanize package not found"
  echo "Try to install it this way: pip install mechanize"
exit 127
fi
if ! python -c 'import pkgutil; exit(not pkgutil.find_loader("pandas"))'; then
  echo "Missing: pandas package not found"
  echo "Try to install it this way: pip install pandas"
exit 127
fi
if ! python -c 'import pkgutil; exit(not pkgutil.find_loader("numpy"))'; then
  echo "Missing: numpy package not found"
  echo "Try to install it this way: pip install numpy"
exit 127
fi
if ! python -c 'import pkgutil; exit(not pkgutil.find_loader("lxml"))'; then
  echo "Missing: lxml package not found"
  echo "Try to install it this way: pip install lxml"
exit 127
fi
if ! python -c 'import pkgutil; \
  exit(not pkgutil.find_loader("dna_features_viewer"))'; then
  echo "Missing: dna_features_viewer package not found"
  echo "Try to install it this way: pip install dna_features_viewer"
exit 127
fi

echo "Required software is properly installed."


########################################
## Identify Errors in Inputs Required ##
########################################

printf "\nChecking if required inputs are correct...\n"
# Check if directory containing FASTQ files exist
if [ ! -d ${long_fastqfolder} ]; then
  echo "Error: --long_fastqfolder doesn't exist."
  echo "Solution: check if the path to this directory is correct."
  exit 1
fi
if [ ! -d ${short_fastqfolder} ]; then
  echo "Error: --short_fastqfolder doesn't exist."
  echo "Solution: check if the path to this directory is correct."
  exit 1
fi

# Check if only FASTQ files are provided in the two fastq folders
for seqfile in ${long_fastqfolder}/*; do
  # Get the extension of the file, if file is compressed
  if file --mime-type "${seqfile}" | grep -q gzip; then
    filename=${seqfile%.*}
    extension=${filename##*.}
  else # Get the extension of the file, if file is NOT compressed
    extension=${seqfile##*.}
  fi
  # Check if extension is fastq or fq
  extension=$(tr "[:upper:]" "[:lower:]" <<< ${extension})
  if [[ ${extension} != "fq" && ${extension} != "fastq" ]]; then
    echo "Error: --long_fastqfolder should only contain FASTQ files."
    echo "Solution: remove any other file from this directory."
    exit 1
  fi
done
for seqfile in ${short_fastqfolder}/*; do
  # Get the extension of the file, if file is compressed
  if file --mime-type "${seqfile}" | grep -q gzip; then
    filename=${seqfile%.*}
    extension=${filename##*.}
  else # Get the extension of the file, if file is NOT compressed
    extension=${seqfile##*.}
  fi
  # Check if extension is fastq or fq
  extension=$(tr "[:upper:]" "[:lower:]" <<< ${extension})
  if [[ ${extension} != "fq" && ${extension} != "fastq" ]]; then
    echo "Error: --short_fastqfolder should only contain FASTQ files."
    echo "Solution: remove any other file from this directory."
    exit 1
  fi
done
echo "Required inputs seem correct."


########################################
## Identify Errors in Optional Inputs ##
########################################

printf "\nChecking if optional inputs are correct...\n"
# Check if the number of threads is an integer
if ! [[ ${threads} =~ ^[0-9]+$ ]]; then
  echo "Error: --threads is not an integer."
  echo "Solution: remove this optional parameter or use an integer."
  exit 1
fi
echo "Optional inputs seem correct."

# Create folder of outputs
rm -rf ${outn} && \
  mkdir ${outn} && \
  chmod +xwr ${outn}

############################
## Perfrom Reads Assembly ##
############################

printf "\nRunning long-read-first hybrid assembly...\n"
echo "--QC filtering: remove the worst reads until only 500 Mbp remain--"

# Create folder with filtered reads
rm -rf ${outn}/reads_filtered && \
  mkdir ${outn}/reads_filtered && \
  chmod +xwr ${outn}/reads_filtered

# Filter reads with low quality
for fq in ${long_fastqfolder}/* ; do
  fq_name=$(echo ${fq} | sed 's/\.fastq.gz//' | sed 's/.*\///')
  echo ${fq}
  filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 \
    ${fq} | gzip > ${outn}/reads_filtered/${fq_name}_filtered.fastq.gz
done

echo "--Generating assemblies with flye--"

# Create folder with assembly
rm -rf ${outn}/assemblies && \
  mkdir ${outn}/assemblies && \
  chmod +xwr ${outn}/assemblies

# Perform long read assemblies
for fq in ${outn}/reads_filtered/* ; do
  fq_name=$(echo ${fq} | sed 's/\_filtered.fastq.gz//' | sed 's/.*\///')
  echo ${fq}
  flye --nano-hq ${fq} --threads ${threads} \
    --out-dir ${outn}/assemblies/${fq_name}_flye
done

echo "Polishing long-read assemblies with short reads"

# Create folder with polished assemblies
rm -rf ${outn}/polished_assemblies && \
  mkdir ${outn}/polished_assemblies && \
  chmod +xwr ${outn}/polished_assemblies

# Polish assemblies with short reads
idx=0
short_reads=(${short_fastqfolder}/*)
for assembly_folder in ${outn}/assemblies/* ; do
  fq=$(echo $assembly_folder | sed 's/\_flye//' | sed 's/.*\///')
  echo ${fq}
  # Get the reference strain for later analysis (the first strain)
  bwa index ${assembly_folder}/assembly.fasta
  echo ${short_reads[${idx}]}
  echo ${assembly_folder}
  echo ${short_reads[${idx}+1]}
  bwa mem -t 16 -a ${assembly_folder}/assembly.fasta \
    ${short_reads[${idx}]} > ${assembly_folder}/alignments_1.sam
  bwa mem -t 16 -a ${assembly_folder}/assembly.fasta \
    ${short_reads[${idx}+1]} > ${assembly_folder}/alignments_2.sam
  polypolish_insert_filter.py --in1 ${assembly_folder}/alignments_1.sam \
    --in2 ${assembly_folder}/alignments_2.sam \
    --out1 ${assembly_folder}/filtered_1.sam \
    --out2 ${assembly_folder}/filtered_2.sam
  polypolish ${assembly_folder}/assembly.fasta \
    ${assembly_folder}/filtered_1.sam ${assembly_folder}/filtered_2.sam > \
    ${outn}/polished_assemblies/${fq}_polished_assembly.fasta
  idx=${idx}+2
done

echo "Hybrid assemblies finished successfully."


###############################
## Compute Prokka annotation ##
###############################

printf "\nAnnotating polished assemblies...\n"
for fq in ${outn}/polished_assemblies/* ; do
  fq_name=$(echo ${fq} | sed 's/\_polished_.*//' | sed 's/.*\///')
  echo ${fq}
  prokka --outdir ${outn}/annotation --prefix ${fq_name}_prokka \
  --cpus ${threads} --force ${fq}
done
echo "Gene annotations finished successfully."


############################
## Compare genome content ##
############################

printf "\nComparing genome content...\n"
roary -e --mafft -p ${threads} -f ${outn}/GenomeContent \
  -r ${outn}/annotation/*.gff
echo "Genome content analysis finished sucessfully."


###########################
# Compute SNPs detection ##
###########################

printf "\nPerforming variant calling analysis...\n"

rm -rf ${outn}/SNPs && \
  mkdir ${outn}/SNPs && \
  chmod +xwr ${outn}/SNPs

# Compute variant calling with Snippy and identify SNPs
idx=0
for polished_file in ${outn}/polished_assemblies/* ; do
  fq=$(echo ${polished_file} | sed 's/\_polished_.*//' | sed 's/.*\///')
  # Use as a reference genome the first strain that is provided
  if [ ${idx} -eq 0 ]; then
    reference_genome="${outn}/annotation/${fq}_prokka.ffn"
  else
    snippy --report --cpus ${threads} \
    --outdir ${outn}/SNPs/${fq} --ref ${reference_genome} \
    --ctgs ${polished_file}
  fi
  idx=${idx}+1
done

# If there are more than 2 strains to compare, identify core SNPs.
end=$(ls -1q ${long_fastqfolder} | wc -l)
if [ ${end} -gt 2 ]; then
  snippy-core --ref ${reference_genome} \
    --prefix ${outn}/SNPsCore ${outn}/SNPs/*
fi
echo "SNPs were detected sucessfully."

##################################
## Identify insertion sequences ##
##################################

printf "\nIdentifying insertion sequences...\n"

rm -rf ${outn}/IS_found && \
  mkdir ${outn}/IS_found && \
  chmod +xwr ${outn}/IS_found

idx=1
for gbk_file in ${outn}/annotation/*.gbk ; do
  # Use as a reference genome the first strain that is provided
  if [ "${idx}" == "1" ]; then
    reference_genome="${gbk_file}"
  else
    ISCompare.py -q ${gbk_file} -r ${reference_genome} \
      -o ${outn}/IS_found/comparison_${idx} -c -p -I <<< "Y"
  fi
  idx=${idx}+1
done

echo "Done!"
echo "All analyses finished successfully. Good luck with the results!"
