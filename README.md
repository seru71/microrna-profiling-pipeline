
# miRNA profiling pipeline 



## Pipeline

The pipeline consists of (optional) bcl2fastq, Trimmomatic, FastQC, BWA, SAMtools, Picard-tools, and ... 



## Setup

### Pipeline

1. Install Python if you don't have it from before, and a cool python package - `Ruffus` (http://www.ruffus.org.uk/). 
Running jobs on a cluster (PBS, Slurm, etc) requires `drmaa` package. 
You might also need following packages: optparse, logging, shutil

2. Clone the pipeline repository:
`git clone https://github.com/seru71/mirna-profiling-pipeline.git <PIPELINE_HOME>`

3. Change directory to newly created pipeline dir, and select the desired version
```
cd <PIPELINE_HOME>
git checkout v0.1
```


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





