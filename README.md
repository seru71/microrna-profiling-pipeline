
# Micro RNA profiling pipeline 

Micro RNA profiling pipeline for NGS reads with UMI. 

Developed at [BTM](http://biostat.umed.pl), [Medical University of Lodz](http://umed.pl).


### Pipeline

Conceptualy data processing is following:

 1) read trimming and extraction of UMI
 2) mapping to miRBase (must be downloaded and indexed separately)

    1) deduping based on UMI
 	
    2) separate counting of unique and multimapping reads
    
 3) mapping of miRBase-unmapped to entire genome 
 
 	 1. (TBD!) deduping 
   
 	 2. counting features
    
 4) QC stats generation, both on read and alignment levels


### Dependencies

The default path consists of following tools, stitched together with [Ruffus](http://www.ruffus.org.uk):

 - (optional) bcl2fastq
 - FastQC
 - umi_tools for extracting UMI bacodes and deduplication
 - bowtie2
 - SAMTools in numerous places
 - featureCount
 - custom utilities for SAM processing

Supplementary/legacy tasks are implemented for:

 - trimming with Trimmomatic and UMI extraction with a custom script (legacy)
 - mapping with BWA and STAR


### Setup

#### Dockerized mode
Not used. Not tested.

#### Native execution mode

1. Install Python package [Ruffus](http://www.ruffus.org.uk/) and all dependencies. 
Running jobs on a cluster (PBS, Slurm, etc) requires additional `drmaa` package, and is currently disabled.

2. Download [miRBase](ftp://mirbase.org/pub/mirbase/) reference, reference genome and annotation.

3. Clone the pipeline repository and cd to it:
```
git clone https://github.com/seru71/mirna-profiling-pipeline.git <PIPELINE_HOME>
cd <PIPELINE_HOME>
```

5. Edit pipeline_settings.cfg specifing:
 
   - input data paths - either runfolder path or regex pointing to FASTQ files
   - reference-root - path to where indexed reference genome is expected
   - scratch-root - work dir
   - mirbase-fasta - fasta file with miRBase reference sequences, tested on mature human miRNAs
   - genome-fasta - reference for mapping non-miRNA reads
   - annotation-gtf - a file with genome feature annotations
   - adapters-fasta - fasta file with adapters for Trimmomatic (not required in default analysis path)
   - paths to executables



### Usage

The pipeline is run using `pipeline.py` script:

#### Running the script

You can run the script using `python <PIPELINE_HOME>/pipeline.py`.
A list of possible options will be presented. 
The only required argument is `-s SETTINGS_FILE`, which specifies the location of the settings file.

The settings file contains paths to input files, resources (e.g. reference FASTA, adapters file), docker settings, and docker containers to use. 
See an exemplary file for all required options in <PIPELINE_HOME>/pipeline_settings.cfg.

Running specific tasks can be done with `-t` / `--target` argument, e.g:

 - `complete_run` task runs entire pipeline with full-QC
 - `count_mirs` runs data processing and miR counting without QC
 - `qc_mapping` generates mapping stat tables, no miR counting
 - `qc_reads` runs FastQC on unprocessed and UMI-trimmed reads

Dry-run (useful for testing and troubleshooting) is enabled with `-n`. 

Verbosity switch allows to print pipeline flow (in dry-run) and follow progress at different detail-level, e.g (`-vv` - tasks, `-vvv` - jobs in tasks).

In order to run several jobs in parallel use the `-j N`.



#### Outputs

The pipeline script creates following directory structure in SCRATCH-ROOT directory (given in settings), or in 
SCRATCH-ROOT/RUN_ID (if the pipeline is run with bcl2fastq step):

 - `SAMPLE_ID/` - one dir per sample, named after samples found in the input data
 - `fastqs/`    - fastq files (if bcl2fastq conversion is run)
 - `drmaa/`     - slurm scripts created automatically (for debugging purposes)
 - `qc/`        - QC output


After finishing, the sample directories will contain:

 - trimmed FASTQ files
 - BAM files
 - sample's mapping stats


The main results are in SCRATCH-ROOT directory which contains three miRNA count tables in tab-delimited format:

 1) mirna_unique  - contains counts of reads mapping uniqly
 2) mirna_allbest - additionally to (1) contains counts of reads mapping with equal alignment score in more than one location
 3) mirna_single  - same as (2) but only first alignment for each read is kept/counted

Tables with mapping stats for all samples at different processing stages are in the qc/ directory.



#### Typical usage

For running the complete analysis using 12 concurrent threads:

```
pipeline.py -s my_settings.cfg \
        -t complete_run \
        -vvv -j 12
```



