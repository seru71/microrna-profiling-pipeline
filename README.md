
# miRNA profiling pipeline 



## Pipeline

The pipeline consists of 

  bcl2fastq (optional),
  Trimmomatic, 
  umi trimming script (umi_trimmer.py)
  BWA & SAMtools,
  umitools,
  FastQC,



## Setup

### Dockerized executables mode

TBD

### Native execution mode

1. Install Python package `Ruffus` (http://www.ruffus.org.uk/). 
Running jobs on a cluster (PBS, Slurm, etc) requires additional `drmaa` package. 

2. Clone the pipeline repository:
`git clone https://github.com/seru71/mirna-profiling-pipeline.git <PIPELINE_HOME>`

3. Change directory to newly created pipeline dir
```
cd <PIPELINE_HOME>
```
4. Specify in pipeline_settings.cfg:
 
  input data paths - either runfolder path or regex pointing to FASTQ files
  reference-root - path to where BWA indexed reference genome is expected
  reference-fasta - fasta file with reference genome
  scratch-root - work dir
  adapters-fasta - fasta file with adapters for Trimmomatic
  paths to executables



## Usage

The pipeline is run using `pipeline.py` script:

* Running the script

    You can run the script using `python <PIPELINE_HOME>/pipeline.py`.
    A list of possible options will be presented. The only required option is `--pipeline_settings FILE`, 
    which specifies the location of the settings file.
    
    The settings file contains paths to input files, resources 
    (e.g. reference FASTA, adapters file), docker settings, 
    and docker containers to use. See an exemplary file for all required options 
    in <PIPELINE_HOME>/pipeline_settings.cfg.
  
    If you want to follow progress of the script, use the verbose option (`-vvvvv`).
    In order to use multithreading, use the `-j` option (`-j 12`).

* Outputs

    The pipeline script creates following directory structure in scratch-root directory (given in settings), or in 
    scratch-root/RUN_ID (if the pipeline is run with bcl2fastq step):
    	SAMPLE_ID/ - one dir per sample, named after samples found in the input data
    	fastqs/    - fastq files (if bcl2fastq conversion is run)
    	drmaa/     - slurm scripts created automatically (for debugging purposes)
    	qc/        - qc output

    After finishing, the sample directories will contain:
    	- trimmed FASTQ files
    	- a BAM file

* Typical usage

    For running the profiling using 12 concurrent threads:

	pipeline.py --settings /my_dir/pipeline_settings/my_settings.cfg \
				--target complete_run \
				-vvv -j 12 \
				--log_file /my_dir/pipeline_run_XXX.log





