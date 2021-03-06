#!/usr/bin/env python
"""

	    pipeline.py
                        --run_folder PATH
                        [--settings PATH] (by default RUN_FOLDER/settings.cfg)
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]  (by default the last task in the pipeline)
                        [--jobs N]        (by default 1)
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]
                        [--run_on_bcl_tile TILE_REGEX]

"""
import sys
import os
import glob


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --settings pipeline_settings.cfg [--target_task TASK] [more_options]")
    
                                
    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", 
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")


    #
    #   pipeline
    #
    parser.add_option("-s", "--settings", dest="pipeline_settings",
                        metavar="FILE",
                        type="string",
                        help="File containing all the settings for the analysis. By default settings.cfg in the run_folder.")                  
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", 
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")

    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true",
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")
    parser.add_option("--rebuild_mode", dest="rebuild_mode",
                        action="store_false", 
                        help="gnu_make_maximal_rebuild_mode")
    parser.add_option("--run_on_bcl_tile", dest="run_on_bcl_tile",
                        type="string",                        
                        help="Use only this tile when doing bcl2fastq conversion. For testing purposes.")
    
    parser.set_defaults(pipeline_settings=None, 
                        jobs=1, verbose=0, 
                        target_tasks=list(), forced_tasks=list(), 
                        just_print=False, key_legend_in_graph=False,
                        rebuild_mode=True, run_on_bcl_tile=None)
    

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = ['pipeline_settings']

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have been defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)
    
        if not len(missing_options):
            return
    
        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
        
    check_mandatory_options(options, mandatory_options, helpstr)
            
    



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888



if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    module_name = "exome"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})



    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Get pipeline settings from a config file  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    if not os.path.exists(options.pipeline_settings): 
	raise Exception('Provided settings file [%s] does not exist or cannot be read.' % options.pipeline_settings)

    import ConfigParser
    config = ConfigParser.ConfigParser()
    config.read(options.pipeline_settings)



    # Should dockerized execution be used?
    dockerize = True
    try:
        docker_bin = config.get('Docker','docker-binary')
        logger.info('Found docker-binary setting. Using dockerized execution mode.')
    except ConfigParser.NoOptionError:
        logger.info('Docker-binary setting is missing. Using regular execution mode.')
        dockerize = False
 

    # Get the pipeline input
    run_folder = None
    input_fastqs = None
    input_bams = None
    try:
	run_folder = config.get('Inputs','run-directory') 
	logger.info('Found run-directory setting. Starting from bcl2fastq conversion.')
        # check presence of the run folder, and sample sheet file
        if not os.path.exists(run_folder) or not os.path.exists(os.path.join(run_folder,'SampleSheet.csv')):
            raise Exception("Missing sample sheet file: %s.\n" % os.path.join(run_folder,'SampleSheet.csv'))
    except ConfigParser.NoOptionError:
        try:
            input_fastqs = os.path.join(runs_scratch_dir if dockerize else '', config.get('Inputs','input-fastqs'))
            input_fastqs_resolved = glob.glob(input_fastqs)
            if len(input_fastqs_resolved) < 2:
                raise Exception("Missing input FASTQs. Found %s FASTQ files in [%s].\n" % (len(input_fastqs_resolved), config.get('Inputs','input-fastqs')))
            logger.info('Found %s input FASTQs. Starting from read trimming.' % len(input_fastqs_resolved))
        except ConfigParser.NoOptionError:
            try:
                input_bams = os.path.join(runs_scratch_dir if dockerize else '', config.get('Inputs','input-bams'))
                input_bams_resolved = glob.glob(input_bams)
                if len(input_bams_resolved) < 1:
                    raise Exception("Missing input BAMs. Found %s BAM files in [%s].\n" % (len(input_bams_resolved), config.get('Inputs','input-bams')))
                logger.info('Found %s input BAMs. Starting from indexing BAMs.' % len(input_bams_resolved))
            except ConfigParser.NoOptionError:
	        raise Exception('Found no valid input setting in [%s]. Please provide one of [run_directory|input-fastqs|input-bams] in the pipeline settings file.' % options.pipeline_settings)


 


    # Root paths
    
    reference_root = config.get('Paths','reference-root')
    scratch_root = os.getcwd()
    try:
        scratch_root   = config.get('Paths','scratch-root')
    except ConfigParser.NoOptionError:
        logger.info('Scratch-root setting is missing. Using current directory: %s' % scratch_root)
    
    run_id = os.path.basename(run_folder) if run_folder != None else os.path.basename(scratch_root)
    runs_scratch_dir = os.path.join(scratch-root, run_id) if run_folder != None else scratch_root
    logger.info('Run\'s scratch directory: %s' % runs_scratch_dir)
      
   
    # optional /tmp dir
    tmp_dir = None
    try:
        tmp_dir = config.get('Paths','tmp-dir')
    except ConfigParser.NoOptionError:
        logger.info('No tmp-dir provided. %s\'s /tmp will be used.' % ('Container' if dockerize else 'Execution host'))
    

    if dockerize:
        # Docker args
        docker_args = config.get('Docker', 'docker-args')
        docker_args += " -v " + ":".join([run_folder, run_folder,"ro"])
        docker_args += " -v " + ":".join([reference_root,reference_root,"ro"])
    
        # Tmp, if should be different than the default  
        if tmp_dir != None: 
            docker_args += " -v " + ":".join([tmp_dir,tmp_dir,"rw"])
            
        docker_args += " -v " + ":".join([runs_scratch_dir,runs_scratch_dir,"rw"])
        docker_args += " -w " + runs_scratch_dir
        docker = " ".join([docker_bin, docker_args]) 
    
    # set the default value if the tmp-dir was unset
    tmp_dir = "/tmp" if tmp_dir==None else tmp_dir
      
    
    # reference files
    mirbase_reference = os.path.join(reference_root, config.get('Resources','mirbase-fasta'))
    genome_reference  = os.path.join(reference_root, config.get('Resources','genome-fasta'))
    annotation_gtf    = os.path.join(reference_root, config.get('Resources','annotation-gtf'))
    
    try:
        adapters = os.path.join(reference_root, config.get('Resources', 'adapters-fasta'))
    except ConfigParser.NoOptionError:
        if run_folder != None or input_fastqs != None: raise Exception('Found no adapters-fasta setting, which is required for read trimming.')
    
    # tools
    bcl2fastq = config.get('Tools','bcl2fastq')
    fastqc = config.get('Tools','fastqc')
    trimmomatic = config.get('Tools', 'trimmomatic') 
    bwa = config.get('Tools','bwa')
    bowtie2 = config.get('Tools','bowtie2')
    star = config.get('Tools','star')
    samtools = config.get('Tools','samtools')
    umitools = config.get('Tools','umitools')
    picard = config.get('Tools','picard-tools')
    featurecounts = config.get('Tools','featurecounts')
        


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#import drmaa
#drmaa_session = drmaa.Session()
#drmaa_session.initialize()
drmaa_session=None

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

   
"""
cmd is given in a form:
    
    command {args}
    interpreter {interpreter_args} command {atgs}

The strings {args} and {interpreter_args} are replaced with args and interpreter_args values.
Examples of correct commands:
    cmd = "samtools {args}"
    cmd = "samtools -p param {args}"
    cmd = "samtools -p2 param {args} -p2 param2"
    cmd = "java {interpreter_args} -jar myjarfile.jar {args} -extras extra_arg
    cmd = "java -XmX4G {interpreter_args} -jar myjarfile.jar {args} -extras extra_arg
"""
def run_cmd(cmd, args, dockerize, interpreter_args=None, run_locally=True,
            cpus=1, mem_per_cpu=1024, walltime='24:00:00', 
            retain_job_scripts = True, job_script_dir = os.path.join(runs_scratch_dir, "drmaa")):
    
    full_cmd = ("{docker} "+cmd).format(docker = docker if dockerize else "",
                                        args=args, 
                                        interpreter_args = interpreter_args if interpreter_args!=None else "")
#    print full_cmd
    stdout, stderr = "", ""
    slurm_account = "default"

    job_options = "--account={account} \
		   --ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(account=slurm_account, cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
                   
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=logger, working_directory=os.getcwd(),
                                 drmaa_session = drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", cmd, err, stdout, stderr])))


""" 
Currently not available in dockerized mode. 
Only default job scheduling params of run_command available when executing via SLURM.
"""
def run_piped_command(*args):
    run_locally=True
    retain_job_scripts = True
    job_script_dir = os.path.join(runs_scratch_dir, "drmaa")	
    cpus=1
    mem_per_cpu=1024
    walltime="24:00:00"
 
    stdout, stderr = "", ""
    job_options = "--ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
	
    full_cmd = expand_piped_command(*args)
#    print full_cmd	
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=logger, working_directory=os.getcwd(),
                                 drmaa_session = drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", full_cmd, err, stdout, stderr])))
	
def expand_piped_command(cmd, cmd_args, interpreter_args=None, *args):
	expanded_cmd = cmd.format(args=cmd_args, interpreter_args = interpreter_args if interpreter_args!=None else "")
	expanded_cmd += (" | "+expand_piped_command(*args)) if len(args) > 0 else ""
	return expanded_cmd



def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd(samtools, "index %s" % bam, dockerize=dockerize)


def produce_fastqc_report(fastq_file, output_dir=None):
    args = fastq_file
    args += (' -o '+output_dir) if output_dir != None else ''
    run_cmd(fastqc, args, dockerize=dockerize)


def count_fastq_reads(fastqs, output):
    run_piped_command("zcat {args}", " ".join(fastqs), None,
                      "grep {args}", "-c ^@ > %s" % output, None)

def transpose_tsv_table(inputf, outputf, columns):
    if os.path.exists(outputf):
        os.remove(outputf)
    for c in columns:
        run_piped_command("cut {args}", "-f %s %s" % (str(c), inputf), None,
                          "xargs {args}", "echo", None,
                          "sed {args}", "'s/ /\t/g' >> %s" % outputf, None)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *


#
#
# Prepare FASTQ
# 

@active_if(run_folder != None)
@follows(mkdir(runs_scratch_dir), mkdir(os.path.join(runs_scratch_dir,'fastqs')))
@files(run_folder, os.path.join(runs_scratch_dir,'fastqs','completed'))
@posttask(touch_file(os.path.join(runs_scratch_dir,'fastqs','completed')))
def bcl2fastq_conversion(run_directory, completed_flag):
    """ Run bcl2fastq conversion and create fastq files in the run directory"""
    out_dir = os.path.join(runs_scratch_dir,'fastqs')
    interop_dir = os.path.join(out_dir,'InterOp')

    # r, w, d, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual) 
    args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -d2 -p4 \
            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
    if options.run_on_bcl_tile != None:
        args += " --tiles %s" % options.run_on_bcl_tile
        
    run_cmd(bcl2fastq, args, dockerize=dockerize, cpus=8, mem_per_cpu=2048)
    


#
# Prepare directory for every sample and link the input fastq files
# Expected format:
#    /path/to/file/[SAMPLE_ID]_S[1-9]\d?_L\d\d\d_R1_001.fastq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@active_if(run_folder != None or input_fastqs != None)
@follows(mkdir(runs_scratch_dir))
@jobs_limit(1)    # to avoid problems with simultanous creation of the same sample dir
@subdivide(os.path.join(runs_scratch_dir,'fastqs','*.fastq.gz') if run_folder != None else input_fastqs,
           formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_L\d\d\d_R1_001\.fastq\.gz$'), 
           os.path.join(runs_scratch_dir,'{SAMPLE_ID[0]}/{basename[0]}{ext[0]}'))
def link_fastqs(fastq_in, fastq_out):
    """Make working directory for every sample and link fastq files in"""
    if not os.path.exists(os.path.dirname(fastq_out)):
        os.mkdir(os.path.dirname(fastq_out))
    if not os.path.exists(fastq_out):
        os.symlink(fastq_in, fastq_out) 

   

    #88888888888888888888888888888888888888888888888888
    #
    #   R e a d   t r i m m i n g
    #
    #88888888888888888888888888888888888888888888888888


PRE_UMI_ADAPTER='AACTGTAGGCACCATCAAT'

# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[S_NUM]_[LANE_ID]_R1_001.fastq.gz
# The output will be written to:
#    [SAMPLE_ID]_[LANE_ID].fq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"

@active_if(run_folder != None or input_fastqs != None)
@transform(link_fastqs, regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R1_001\.fastq\.gz$'),  r'\1/\2_\3.fq.gz')
def extract_umi(input_fq, output_fq):
    """ Trim reads to the beginning of UMI adaper and extract UMI sequence to read identifier """ 
    
    pattern = ".{15}.*(?P<discard_1>%s){s<=2}(?P<umi_1>.{12})(?P<discard_2>.*)" % PRE_UMI_ADAPTER
    logfile = output_fq+'.log'
    
    args = "extract -I {infq} -S {outfq} -L {log} \
                    --extract-method='regex' \
                    --bc-pattern '{pattern}' \
                    ".format(infq=input_fq, 
                             outfq=output_fq, 
                             log=logfile, 
                             pattern=pattern)
    
    run_cmd(umitools, args, dockerize=dockerize)
 
    
#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[S_NUM]_[LANE_ID]_R1_001.fastq.gz
# The output will be written to:
#    [SAMPLE_ID]_[LANE_ID].fq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@active_if(run_folder != None or input_fastqs != None)
@transform(link_fastqs, regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R1_001\.fastq\.gz$'),  r'\1/\2_\3.fq.gz')
def trim_reads_and_extract_umi(input_fq, output_fq):

    intermediate_fq = output_fq+".intermediate.gz"
    
    # trim adapters
    args = "SE -phred33 -threads 1 \
            {in_fq} {out_fq} \
            ILLUMINACLIP:{adapter}:2:30:8 \
            SLIDINGWINDOW:4:15".format(in_fq=input_fq, out_fq=intermediate_fq, adapter=adapters)

#    # if umis are used, set minimum length of reat after adapter trimming to 12nt + lem(pre_umi_adapter) + umi_length
#    if PRE_UMI_ADAPTER:
#        minumilen=10
#        minlen = 12 + len(PRE_UMI_ADAPTER) + minumilen
#        args += " MINLEN:%s" % minlen

    max_mem = 2048
    run_cmd(trimmomatic, args, interpreter_args="-Xmx"+str(max_mem)+"m", 
            dockerize=dockerize, cpus=1, mem_per_cpu=max_mem)

    if PRE_UMI_ADAPTER:
        # trim UMIs 
        umi_stats_file = output_fq+'.umistats'
        args = " ".join([intermediate_fq, output_fq, PRE_UMI_ADAPTER, ">", umi_stats_file])
        run_cmd("python umi_trimmer.py {args}", args, dockerize=dockerize)
    else:
        os.symlink(intermediate_fq, output_fq)





    #88888888888888888888888888888888888888888888888888
    #
    # QC the FASTQ files
    #
    #88888888888888888888888888888888888888888888888888


@follows(mkdir(os.path.join(runs_scratch_dir,'qc')), mkdir(os.path.join(runs_scratch_dir,'qc','read_qc')))
@transform(link_fastqs, formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fastq\.gz$'), 
           os.path.join(runs_scratch_dir,'qc','read_qc/')+'{SAMPLE_ID[0]}_fastqc.html')
def qc_raw_reads(input_fastq, report):
    """ Generate FastQC report for raw FASTQs """
    produce_fastqc_report(input_fastq, os.path.dirname(report))


@follows(mkdir(os.path.join(runs_scratch_dir,'qc')), mkdir(os.path.join(runs_scratch_dir,'qc','read_qc')))
@transform(extract_umi, formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$'), 
	  os.path.join(runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html')
def qc_trimmed_reads(input_fastq, report):
    """ Generate FastQC report for trimmed FASTQs """
    produce_fastqc_report(input_fastq, os.path.dirname(report))


@follows(mkdir(os.path.join(runs_scratch_dir,'qc')), mkdir(os.path.join(runs_scratch_dir,'qc','read_qc')))
@collate(link_fastqs, 
            formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_L\d\d\d\_R1_001.fastq\.gz$"), 
            "{subpath[0][1]}/qc/read_qc/{subdir[0][0]}.raw_readcount")
def count_raw_reads_per_sample(fastqs, count_file):
    count_fastq_reads(fastqs, count_file)


@follows(mkdir(os.path.join(runs_scratch_dir,'qc')), mkdir(os.path.join(runs_scratch_dir,'qc','read_qc')))
@collate(extract_umi, 
        formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_L\d\d\d\.fq\.gz$"), 
        "{subpath[0][1]}/qc/read_qc/{subdir[0][0]}.trimmed_readcount")
def count_trimmed_reads_per_sample(fastqs, count_file):
    count_fastq_reads(fastqs, count_file)

@merge([count_raw_reads_per_sample, count_trimmed_reads_per_sample], 
        os.path.join(runs_scratch_dir,'qc', 'read_qc', 'all_samples.readcounts'))
def join_fastq_readcounts(count_files, output):
    raw = count_files[0:len(count_files)/2]
    trimmed = count_files[len(count_files)/2:]
    with open(output+".tmp", 'w') as out:
        out.write("sample\traw\ttrimmed\n")
        for i in range(0,len(raw)):
            sample = os.path.basename(raw[i])[:-len(".raw_readcount")]
            with open(raw[i]) as r, open(trimmed[i]) as t:
                out.write('\t'.join([sample, r.readline().strip(), t.readline()]))

    transpose_tsv_table(output+".tmp", output, [1,2,3])
    os.remove(output+".tmp")

@follows(qc_raw_reads, qc_trimmed_reads, join_fastq_readcounts)
def qc_reads():
    pass


    #8888888888888888888888888888888888888888888888888888
    #
    #                   M a p p i n g 
    #
    #8888888888888888888888888888888888888888888888888888



def bwa_map_and_sort_pe(output_bam, ref_genome, fq1, fq2, read_group=None, threads=1):
    raise Exception('Not implemented')
#
#   NOT TESTED!
#
#	bwa_aln_args = "aln -l10 -k1 -t {threads} {ref} {fq1} \
#			".format(threads=threads, 
#	                         ref=ref_genome, fq1=fq1)
#
#	bwa_sampe_args = "sampe {rg} {ref} - {fq1} {fq2} \
#	            ".format(rg="-r '%s'" % read_group if read_group!=None else "", 
#                        ref=ref_genome, fq1=fq1, fq2=fq2)
#
#	samtools_args = "sort -o {out} -".format(out=output_bam)
#
#	run_piped_command(bwa, bwa_aln_args, None,
#                          bwa, bwa_sampe_args, None,
#                          samtools, samtools_args, None)



def bwa_map_and_sort_se(output_bam, ref_genome, fq1, read_group=None, threads=1):
	
	#bwa_aln_args = "aln -l10 -k1 -t {threads} {ref} {fq1} \
    
    bwa_aln_args = "aln -n1 -t {threads} {ref} {fq1} \
                   ".format(threads=threads, 
                            ref=ref_genome, fq1=fq1)

    bwa_samse_args = "samse {rg} {ref} - {fq} \
                     ".format(rg="-r '%s'" % read_group if read_group!=None else "", 
                        ref=ref_genome, fq=fq1)
	
    samtools_args = "sort -o {out} -".format(out=output_bam)

    run_piped_command(bwa, bwa_aln_args, None,
                          bwa, bwa_samse_args, None,
                          samtools, samtools_args, None)
	


def bowtie_map_and_sort_se(output_bam, ref_genome, fq1, read_group=None, threads=1):
    
    rgs = read_group.split('\\t')
    rg_id = rgs[1]
    rg_fields = ' '.join(['--rg '+s for s in rgs[2:]])
    # --score-min G,28,0
    bowtie2_args = '-N1 -L9 --norc -k10 --local \
                    --score-min L,4,1.3 --mp 4 \
                    --rg-id {rg_id} {rg} \
                    -x {ref} -U {fq}'.format(ref=ref_genome, fq=fq1,
                                            rg_id=rg_id, rg=rg_fields)
    convert_args = "-t convert-secondary"
    samtools_args = "view -h -F256"
    sort_args = "sort -o {out} -".format(out=output_bam)
    
    run_piped_command(bowtie2, bowtie2_args, None,
                      os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mirbase_prep.py {args}'), convert_args, None,
                      samtools, samtools_args, None,
                      samtools, sort_args, None)
    
    
def star_map_and_sort_se(output_bam, ref_genome, fq, read_group=None, threads=1, 
                         all_best_primary=False, keep_unmapped=True):
                             
    ref_dir = ref_genome+"_STAR"
    star_args = "--genomeDir {ref} --readFilesIn {fq} --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate \
                --outFileNamePrefix {bam} --outFilterMismatchNmax {nm} --outFilterScoreMin 14 \
                ".format(ref=ref_dir, fq=fq, bam=output_bam, nm=1) 
                                                   
    if all_best_primary: star_args+=" --outSAMprimaryFlag AllBestScore"
    if keep_unmapped: star_args += " --outSAMunmapped Within"
    
    samtools_args = 'view -hb -F16 -F256 > %s' % output_bam
    
    run_piped_command(star, star_args, None,
                        samtools, samtools_args, None)
    

def merge_bams(out_bam, *in_bams):
	threads = 1
	mem = 4096
	
	args = "merge %s" % out_bam
	for bam in in_bams:
		args += (" "+bam)
		
	run_cmd(samtools, args, dockerize=dockerize, cpus=threads, mem_per_cpu=int(mem/threads))
	
	
def map_reads(fastq_list, ref_genome, output_bam, read_groups=None, mapper='bowtie', **args):

    # If no read groups is provided, we could make up default ones based on fq filenames.
    if read_groups==None:
        s_ids = [ os.path.basename(s[0][:-len('_R1.fq.gz')] if isinstance(s, tuple) else s[:-len('fq.gz')]) for s in fastq_list ]
        read_groups = [ '@RG\\tID:{sid}\\tSM:{sid}\\tLB:{sid}'.format(sid=s) for s in s_ids]
    
    tmp_bams = [ output_bam+str(i) for i in range(0, len(fastq_list)) ]
    for i in range(0, len(fastq_list)):
        if isinstance(fastq_list[i], tuple):
            if mapper == 'bwa':
                bwa_map_and_sort_pe(tmp_bams[i], ref_genome, fastq_list[i][0], fastq_list[i][1], read_groups[i])
            else:
                raise Exception('Bowtie2/STAR PE mapping is not implemented')
        else:
            if mapper == 'bwa':
                bwa_map_and_sort_se(tmp_bams[i], ref_genome, fastq_list[i], read_group=read_groups[i])   
            elif mapper == 'bowtie':
                bowtie_map_and_sort_se(tmp_bams[i], ref_genome, fastq_list[i], read_group=read_groups[i])   
            else:
               star_map_and_sort_se(tmp_bams[i], ref_genome, fastq_list[i], read_group=read_groups[i], **args)
                
    merge_bams(output_bam, *tmp_bams)
    
    for f in tmp_bams:
        os.remove(f)


@collate(extract_umi, 
            formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_L\d\d\d\.fq\.gz$"), 
            "{subpath[0][0]}/{subdir[0][0]}.bam",
            "{SAMPLE_ID[0]}")
def map_to_mirbase(fastqs, bam_file, sample_id):
    """ Maps trimmed reads from all lanes """
    read_groups = ['@RG\\tID:{rgid}\\tSM:{lb}\\tLB:{lb}'\
			.format(rgid=sample_id+"_"+lane_id, lb=sample_id) for lane_id in ['L001', 'L002', 'L003', 'L004']]
    map_reads(fastqs, mirbase_reference, bam_file, read_groups, mapper='bowtie')    

    
@transform(map_to_mirbase, suffix('.bam'), '.unmapped.fq.gz')
def extract_unmapped(bam, fastq):
    samtools_args = "view -h -f4 {}".format(bam)
    picard_args = "SamToFastq I=/dev/stdin F={} VALIDATION_STRINGENCY=LENIENT".format(fastq)
    
    run_piped_command(samtools, samtools_args, None,
                        picard, picard_args, None)


@transform(extract_unmapped, suffix('.fq.gz'), '.bam')
def map_unmapped(fastq, bam):
    map_reads([fastq], genome_reference, bam, mapper='bwa')



    #888888888888888888888888888888888888888888
    #
    #        U M I   d e d u p i n g
    #
    #88888888888888888888888888888888888888888888



@transform(map_to_mirbase, suffix('.bam'), '.bam.bai')
def index_mirbase_bam(bam, _):
    """ Index mirbase bam """
    index_bam(bam)


@follows(index_mirbase_bam)
@transform(map_to_mirbase, suffix('.bam'), '.dedup.bam', r'\1.dedup.log')
def dedup_by_umi(input_bam, dedupped_bam, logfile):
    """ Dedup reads with the same mapping start and UMI """
    args = "dedup -I {inbam} -S {outbam} -L {log} \
            --method unique \
	       ".format(inbam=input_bam, outbam=dedupped_bam, log=logfile)
    run_cmd(umitools, args, dockerize=dockerize)



    #888888888888888888888888888888888888888888
    #
    #        C o u n t i n g
    #
    #88888888888888888888888888888888888888888888
    
'''
    
    BWA aln: 
     - does not allow single-strand mapping
     - keeps alternative hits, also with equal score in XA tag, not as separate records
     - does not soft-clip outside of the reference (reference sequences need to be padded)
      
     Cannot be used because
     - due to (1) cannot map directional reads
     - due to (2) cannot filter on rc-strand flag because that can lead to removing equally scoring fwd-strand hit in XA-tag
     - filtering of SAM would be necessary to extract the correctly mapping reads, but XA tag does not contain alignment score and selecting best scoring alignments would be probelematic
    
    Bowtie2:
     + allows for one-strand mapping (--norc flag)
     + keeps alternative hits in separate records
     - only one hit gets primary flag, all alternatives (incl. best scoring) get supplementary hit flag
     + soft-clipping in --local mode
     
     Uniqly mapping reads can be counted without problems:
      - get multimapping output and filter MQ>5
     Best mapping multimappers require SAM parsing and modifying mapping flag prior to filtering and counting:
      - get multimapping output and swap 256 flag for the lines where alignment score (AS) is equal to the primary alignment AS
      - filter on supplementary alignment flag (256)
    
    STAR:
     - does not allow for one-strand mapping
     + keeps alternative hits in separate records, allowing to compensate for (1) with flag-filtering reverse mappers
     + all best hits can get primary flag, lower scoring hits get supplementary alignment flag
     + allows soft-clipping, but does not allow tuning of the soft-clipping extent
     
     Uniqly mapping fwd-reads cannot be easily filtered becasue MQ and NH tag (number of hits) are calculated based on all mapped reads (including reverse strand) 
     Best mapping multimappers can be counted without problems.
    
    
    
'''


def count_ref_hits(bam, counts_file):
    run_piped_command(samtools, "idxstats %s" % bam, None,
                      "cut {args}", "-f1,3 > %s" % counts_file, None)

def merge_count_tables(count_files, table_file, suffix=".dedup.filt.bam.mirbase_counts.txt"):
    sample_ids = [os.path.basename(f)[:-len(suffix)] for f in count_files]
    header = "mirna\t" + '\t'.join(sample_ids) + '\n'
    with open(table_file, 'w') as f:
        f.write(header)
    
    inputs = " ".join(count_files)
    run_cmd("paste {inputs} | \
             cut -f1,$(echo `seq 2 2 {num}` | sed 's/ /,/g') \
             >> {table}".format(inputs=inputs, num=2*len(inputs), table=table_file), "", None)


def filter_alignments(inbam, outbam, filters):
    samtools_view_args = "view {filt} -hb {bam} > {out}\
                         ".format(filt=filters, bam=inbam, out=outbam)
    run_cmd(samtools, samtools_view_args, None)


#
# count uniqly mapping (works well only with bowtie2)

@transform(dedup_by_umi, suffix('.bam'), '.uniq.bam')
def filter_unique(inbam, outbam):
    ''' drop multimapping alignments '''
    # bowtie2 --norc presets
    filter_alignments(inbam, outbam, "-q10")
    
    # STAR presets - reverse-complement and suboptimal mappings are removed in the mapping step
    # WARNING: -q4 can remove reads mapping uniqly to forward-strand, because STAR maps to both 
    #          strands and sets low MQ if a read maps in more than 1 location. There are several 
    #          reverse complement miRs which will be affected by this
    #filter_alignments(inbam, outbam, "-q4")

@transform(filter_unique, suffix('.bam'), '.bam.bai')
def index_deduped_uniq_bam(bam, _):
    """ Index deduped bam """
    index_bam(bam)

@follows(index_deduped_uniq_bam)
@transform(filter_unique, suffix('.bam'), '.bam.mirbase_counts.txt')
def count_unique_mirbase_reads(bam, counts_file):
    """ Count uniqly mapping mirbase hits """
    count_ref_hits(bam, counts_file)

@merge(count_unique_mirbase_reads, os.path.join(runs_scratch_dir, "mirna_unique_count_table.txt"))
def produce_mirna_unique_counts_table(count_files, table_file):
    """ Join per sample count tables for unique counts"""
    merge_count_tables(count_files, table_file, ".dedup.uniq.bam.mirbase_counts.txt")


#
# count multimappers once
# (removes seconds and subsequent best scoring alignments for each read)
#

@transform(dedup_by_umi, suffix('.bam'), '.single.bam')
def filter_single(inbam, outbam):
    ''' drop second and subsequent alignments of each read'''
    
    mirbase_prep=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'mirbase_prep.py {args}')
    run_piped_command(samtools, "view -h %s" % inbam, None,
                      mirbase_prep, "-t drop-multi", None,  
                      samtools, "view -hSb > %s" % outbam, None)

@transform(filter_single, suffix('.bam'), '.bam.bai')
def index_single_bam(bam, _):
    """ Index deduped bam """
    index_bam(bam)

@follows(index_single_bam)
@transform(filter_single, suffix('.bam'), '.bam.mirbase_counts.txt')
def count_single_mirbase_reads(bam, counts_file):
    """ Count uniqly mapping mirbase hits """
    count_ref_hits(bam, counts_file)

@merge(count_single_mirbase_reads, os.path.join(runs_scratch_dir, "mirna_single_count_table.txt"))
def produce_mirna_single_counts_table(count_files, table_file):
    """ Join per sample count tables for unique counts"""
    merge_count_tables(count_files, table_file, ".dedup.single.bam.mirbase_counts.txt")



#
# count all best scoring hits (includes best scoring multimappers)
# STAR could be used if mapping was done with flag editing (secondary&best_score -> primary)

@transform(dedup_by_umi, suffix('.bam'), '.bam.bai')
def index_deduped_bam(bam, _):
    """ Index deduped bam """
    index_bam(bam)

@follows(index_deduped_bam)
@transform(dedup_by_umi, suffix('.bam'), '.bam.mirbase_counts.txt')
def count_allbest_mirbase_reads(bam, counts_file):
    """ Count mirbase hits including equally good multimappers"""
    count_ref_hits(bam, counts_file)

@merge(count_allbest_mirbase_reads, os.path.join(runs_scratch_dir, "mirna_allbest_count_table.txt"))
def produce_mirna_allbest_counts_table(count_files, table_file):
    """ Join per sample count tables for alignments including equally good multimappers"""
    merge_count_tables(count_files, table_file, ".dedup.bam.mirbase_counts.txt")


@follows(produce_mirna_unique_counts_table, produce_mirna_single_counts_table, produce_mirna_allbest_counts_table)
def count_mirs():
    pass


#
# unmapped

@transform(map_unmapped, suffix('.bam'), '.bam.bai')
def index_genome_bam(bam, _):
    """ Index genome bam """
    index_bam(bam)


@follows(index_genome_bam)
@transform(map_unmapped, suffix('.bam'), '.bam.counts')
def count_unmapped_reads_by_category(bam, counts_file):
    """ Count genome features hit by mirbase-unmapped reads """
    args = "-t gene -g gene_biotype -O -M --fraction \
            -a {gtf} -o {out} {bam} \
            ".format(gtf=annotation_gtf, out=counts_file, bam=bam)

    run_cmd(featurecounts, args, dockerize=dockerize)
    
    # extract counts only
    tmp_name = counts_file+".list"
    run_cmd("cut -f1,7 {} > {} \
            ".format(counts_file,tmp_name) , "", dockerize=dockerize)
    os.rename(tmp_name, counts_file)


    #888888888888888888888888888888888888888888
    #
    #        Q C   t a b l e s 
    #
    #88888888888888888888888888888888888888888888



def run_flagstat(bam):
    run_cmd(samtools, 
            "flagstat {} > {} ".format(bam, bam+'.flagstat'), 
            dockerize=dockerize)   
    
@transform(map_to_mirbase, suffix(".bam"), ".bam.flagstat")
def flagstat_mirbase_bam(bam, _):
    run_flagstat(bam)
    
@transform(dedup_by_umi, suffix(".bam"), ".bam.flagstat")
def flagstat_deduped_allbest_bam(bam, _):
    run_flagstat(bam)

@transform(filter_unique, suffix(".bam"), ".bam.flagstat")
def flagstat_deduped_unique_bam(bam, _):
    run_flagstat(bam)
    
@transform(map_unmapped, suffix(".bam"), ".bam.flagstat")
def flagstat_unmapped_bam(bam, _):
    run_flagstat(bam)    
    
def list_to_tsv_line(l):
    return('\t'.join(l) + '\n')

def get_total_and_mapped_from_flagstats(flagstat_files, file_suffix='.bam.flagstat', header_list=['sample','total','mapped']):
    out = list_to_tsv_line(header_list)
    for fname in flagstat_files:
        sample_id = os.path.basename(fname)[:-len(file_suffix)]
        with open(fname) as f:
            top_lines = f.readlines()[0:5]
            lines = [top_lines[i] for i in [0,4]]
            total, mapped = [l.split(" ")[0] for l in lines]
            out += list_to_tsv_line([sample_id, total, mapped])
    return out
    
'''
@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(flagstat_mirbase_bam, os.path.join(runs_scratch_dir, "qc", "mirbase_mapping_stats.tsv"))
def aggregate_mirbase_mapping_stats(flagstats, out_table):
    with open(out_table,"wt") as out:
        out.write(get_total_and_mapped_from_flagstats(flagstats, 
                                                      header_list=['sample','mirbase_total','mirbase_mapped']))

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(flagstat_deduped_allbest_bam, os.path.join(runs_scratch_dir, "qc", "deduped_allbest_mapping_stats.tsv"))
def aggregate_deduped_allbest_mapping_stats(flagstats, out_table):
     with open(out_table,"wt") as out:
        out.write(get_total_and_mapped_from_flagstats(flagstats, file_suffix='.dedup.bam.flagstat',
                                                      header_list=['sample','allbest_total','allbest_mapped']))

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(flagstat_deduped_unique_bam, os.path.join(runs_scratch_dir, "qc", "deduped_uniq_mapping_stats.tsv"))
def aggregate_deduped_unique_mapping_stats(flagstats, out_table):
     with open(out_table,"wt") as out:
        out.write(get_total_and_mapped_from_flagstats(flagstats, file_suffix='.dedup.uniq.bam.flagstat',
                                                      header_list=['sample','unique_total','unique_mapped']))

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(flagstat_unmapped_bam, os.path.join(runs_scratch_dir, "qc", "unmapped_mapping_stats.tsv"))
def aggregate_unmapped_mapping_stats(flagstats, out_table):
     with open(out_table,"wt") as out:
        out.write(get_total_and_mapped_from_flagstats(flagstats, file_suffix='.unmapped.bam.flagstat',
                                                      header_list=['sample','unmapped_total','unmapped_mapped']))
'''
    
#
# mapstats
# Flagstat counts not reads but alignments, so if a read has two alignments it is counted twice. Mapstat counts unique reads in a BAM
#

    
def run_mapstat(bam):
    """ Mapstat counts unique reads in a BAM, all and mapped """
    run_cmd(samtools, "view {} | cut -f1 | sort | uniq | wc -l > {} ".format(bam, bam+'.mapstat'), dockerize=dockerize)
    run_cmd(samtools, "view -F4 {} | cut -f1 | sort | uniq | wc -l >> {} ".format(bam, bam+'.mapstat'), dockerize=dockerize)

@transform(map_to_mirbase, suffix(".bam"), ".bam.mapstat")
def mapstat_mirbase_bam(bam, _):
    run_mapstat(bam)
    
@transform(dedup_by_umi, suffix(".bam"), ".bam.mapstat")
def mapstat_deduped_allbest_bam(bam, _):
    run_mapstat(bam)

@transform(filter_unique, suffix(".bam"), ".bam.mapstat")
def mapstat_deduped_unique_bam(bam, _):
    run_mapstat(bam)
    
@transform(filter_single, suffix(".bam"), ".bam.mapstat")
def mapstat_deduped_single_bam(bam, _):
    run_mapstat(bam)

@transform(map_unmapped, suffix(".bam"), ".bam.mapstat")
def mapstat_unmapped_bam(bam, _):
    run_mapstat(bam)    
    

def get_total_and_mapped_from_mapstats(mapstat_files, output_table, file_suffix='.bam.mapstat', header_list=['sample','total','mapped']):
    tmp_file = "/tmp/smpls"
    with open(tmp_file, 'w') as f:
        f.write("\n".join(header_list[1:]))

    with open(output_table, 'w') as f:
        f.write("\t".join([header_list[0]]+[os.path.basename(fname)[:-len(file_suffix)] for fname in mapstat_files]))
        f.write("\n")

    run_cmd("paste {args}","%s >> %s" % (" ".join([tmp_file]+mapstat_files), output_table), dockerize=dockerize)


@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(mapstat_mirbase_bam, os.path.join(runs_scratch_dir, "qc", "mirbase_mapping_stats.tsv"))
def aggregate_mirbase_mapping_stats(mapstats, out_table):
    get_total_and_mapped_from_mapstats(mapstats, out_table, header_list=['sample','mirbase_total','mirbase_mapped'])

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(mapstat_deduped_allbest_bam, os.path.join(runs_scratch_dir, "qc", "deduped_allbest_mapping_stats.tsv"))
def aggregate_deduped_allbest_mapping_stats(mapstats, out_table):
    get_total_and_mapped_from_mapstats(mapstats, out_table, file_suffix='.dedup.bam.mapstat',
                                                      header_list=['sample','allbest_total','allbest_mapped'])

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(mapstat_deduped_unique_bam, os.path.join(runs_scratch_dir, "qc", "deduped_uniq_mapping_stats.tsv"))
def aggregate_deduped_unique_mapping_stats(mapstats, out_table):
    get_total_and_mapped_from_mapstats(mapstats, out_table, file_suffix='.dedup.uniq.bam.mapstat',
                                                      header_list=['sample','unique_total','unique_mapped'])

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(mapstat_deduped_single_bam, os.path.join(runs_scratch_dir, "qc", "deduped_single_mapping_stats.tsv"))
def aggregate_deduped_single_mapping_stats(mapstats, out_table):
    get_total_and_mapped_from_mapstats(mapstats, out_table, file_suffix='.dedup.single.bam.mapstat',
                                                      header_list=['sample','single_total','single_mapped'])

@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(mapstat_unmapped_bam, os.path.join(runs_scratch_dir, "qc", "unmapped_mapping_stats.tsv"))
def aggregate_unmapped_mapping_stats(mapstats, out_table):
    get_total_and_mapped_from_mapstats(mapstats, out_table, file_suffix='.unmapped.bam.mapstat',
                                                      header_list=['sample','unmapped_total','unmapped_mapped'])


# unmapped categories



@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(count_unmapped_reads_by_category, os.path.join(runs_scratch_dir, "qc", "unmapped_biofeature_stats.tsv"))
def aggregate_unmapped_biofeatures(input_files, out_table):
    header = "sample"
    content = ""
    for i, fname in enumerate(input_files):
        sample_id = os.path.basename(fname)[:-len(".unmapped.bam.counts")]
        content += sample_id
        with open(fname) as f:
            for l in f.xreadlines():
                if l[0]=='#' or l.find("Geneid") == 0: continue
                cat, count = l.strip().split()
                content += '\t'+count
                if i == 0: header += '\t'+cat
        content += '\n'
    
    with open(out_table, 'wt') as f:
        f.write(header + '\n' + content)


@transform(join_fastq_readcounts, formatter(), 
          add_inputs(aggregate_mirbase_mapping_stats, aggregate_deduped_allbest_mapping_stats, \
                     aggregate_deduped_unique_mapping_stats, aggregate_deduped_single_mapping_stats, \
                     aggregate_unmapped_mapping_stats, aggregate_unmapped_biofeatures), 
          os.path.join(runs_scratch_dir, "qc", "mapping_stats.tsv"))
def join_mapping_stats(input_stats, out_stats):
    
    
    fq_stats, mirbase_stats, allbest_stats, unique_stats, single_stats, unmapped_stats, bf_stats = \
        input_stats[0], input_stats[1], input_stats[2], input_stats[3], input_stats[4], input_stats[5], input_stats[6]
    
    # cat mapstats first
    run_piped_command("cat {args}", "{} {} {} {} {} ".format(mirbase_stats, allbest_stats, unique_stats, single_stats, unmapped_stats), None,
                       "grep {args}", "-v sample", None,
                       "cat {args}", "{} - > {}".format(fq_stats, out_stats), None)
    
    
    # drop sample id from all except mirbase_stats
    # drop dedup_mapped because it is equal to dedup_total
    #args = '{} {} {} {} {} {} | cut -f 1-3,5,8,11,14-15,17- > {}\
    #       '.format(mirbase_stats, allbest_stats, unique_stats, single_stats, unmapped_stats, bf_stats, out_stats)
    #run_cmd('paste {args}', args, dockerize=dockerize)


@follows(join_mapping_stats)
def qc_mapping():
    pass



#@posttask(cleanup_files)
@follows(count_mirs, qc_reads, qc_mapping)
def complete_run():
    pass





#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            verbose=options.verbose, #verbose_abbreviated_path=0,
                            checksum_level = 0)

    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                        gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:        
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multithread     = options.jobs,
#                            logger          = stderr_logger,
                            logger          = logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)
    
        
    #drmaa_session.exit()
    
