import groovy.json.JsonOutput
import java.time.format.DateTimeFormatter

def fastp_version = '1.3.1'
def bowtie2_version = '2.5.5'
def samtools_version = '1.23.1'


params.skip_fastp = false
params.fastq_ext  = ".fastq.gz"
params.fwd_suffix = "_1"
params.rev_suffix = "_2"
params.n_base_limit = 0
params.length_required = 60
def available_adapters = [
    'Illumina_TruSeq' : [
        fwd: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
        rev: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        fastp_opts: '--trim_poly_g'
    ],
    'Illumina_Nextera' : [
        fwd: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
        rev: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
        fastp_opts: '--trim_poly_g'
    ],
    'DNBseq' : [
        fwd: 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA',
        rev: 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG',
        fastp_opts: null
    ]
]
params.adapters = null
params.fastp_opts = ""

params.skip_host_removal = false
def available_hosts = [
    'Homo_sapiens',
    'Sus_scrofa',
    'Mus_musculus',
    'Rattus_norvegicus',
    'Oryctolagus_cuniculus',
    'Gallus_gallus',
    'Bos_taurus',
    'Canis_lupus_familiaris',
    'Oncorhynchus_mykiss',
    'Felis_catus',
    'Equus_caballus'
]
params.host = null
params.host_index_parent_dir = '/biodatabank/meteor-data/reference.ng'
params.save_host_reads = false
params.host_reads_suffix = '_Host'
params.host_subdir= 'host'

params.sample_list = null
params.project_name = null
params.help = false

if (params.help || !workflow.commandLine.contains('--')) {
    log.info """
USAGE:
    nextflow run ${workflow.scriptName} --input_dir INPUT_DIR --output_dir OUTPUT_DIR [...]

DESCRIPTION:
    This workflow performs quality control on paired-end or single-end FASTQ files.
    It includes:
      • Read trimming and adapter removal using fastp
      • Removal of host-derived reads by aligning against the host genome

MANDATORY PARAMETERS:
    --input_dir INPUT_DIR             Path to the folder containing sample directories with raw FASTQ files.

    --output_dir OUTPUT_DIR           Path to the folder containing sample directories after quality control.

FASTQ FILES PARAMETERS:
    --fastq_ext FASTQ_EXT             FASTQ files extension.
                                      Default: ${params.fastq_ext}

    --fwd_suffix FWD_SUFFIX           Suffix of forward read files.
                                      Default: ${params.fwd_suffix}

    --rev_suffix REV_SUFFIX           Suffix of reverse read files.
                                      Default: ${params.rev_suffix}

FASTP PARAMETERS:
    --skip_fastp                      Skip read trimming and adapter removal
                                      Default: ${params.skip_fastp}

    --adapters SEQUENCING_ADAPTERS    Name of sequencing adapters to be removed from raw reads
                                      Supported adapters: ${available_adapters.keySet().join(', ')}
                                      Default: ${params.adapters}

    --n_base_limit N_BASE_LIMIT       Maximum number of undetermined bases (Ns) in a read
                                      Default: ${params.n_base_limit}

    --length_required LENGTH          Reads shorter than LENGTH after fastp will be discarded
                                      Default: ${params.length_required}

    --fastp_opts PARAMS               Additional command-line options to pass directly to fastp

HOST REMOVAL PARAMETERS:
    --skip_host_removal               Skip host-derived reads removal
                                      Default: ${params.skip_host_removal}

    --host HOST                       Reads mapped to the HOST genome will be discarded
                                      Supported hosts: ${available_hosts.join(', ')}
                                      Default: ${params.host}

    --save_host_reads                 Save host-derived reads to separate FASTQ files
                                      Default: ${params.save_host_reads}

    --host_reads_suffix HOST_SUFFIX   Suffix inserted into FASTQ filenames with host-derived reads.
                                      Default: ${params.host_reads_suffix}

    --host_subdir DIRNAME             Name of the subdirectory (inside each sample folder) where
                                      FASTQ files with host-derived reads will be saved.
                                      Default: ${params.host_subdir}

OTHER PARAMETERS:
    --sample_list SAMPLE_LIST         Path to a text file listing the samples to process (one sample name per line).
                                      Default: ${params.sample_list}

    --project_name PROJECT_NAME       Name of the sequencing project
                                      Default: Name of the parent directory containing the input_dir

WORKING DIRECTORY:
    By default, Nextflow will create its working directory under ${workDir}
    You can override the default location using the Nextflow -w option.
    IMPORTANT: The working directory should be a NFS-mounted storage accessible from all nodes

EXAMPLE:
    nextflow run ${workflow.scriptName} \\
        -w                  ${System.getenv('HOME')} \\
        --input_dir         /projects/my_project/sample_raw \\
        --output_dir        /projects/my_project/sample_fastp_noHost \\
        --fastq_ext         .fq.gz \\
        --fwd_suffix        _1.clean.rmHost \\
        --rev_suffix        _2.clean.rmHost \\
        --adapters          DNBseq \\
        --n_base_limit      1 \\
        --length_required   60 \\
        --fastp_opts='--trim_poly_x'  \\
        --host              Homo_sapiens \\
        --save_host_reads \\
        --host_reads_suffix '_fromHost' \\
        --host_subdir host_reads
        --project_name      my_project
"""
    exit 0
}


// Mandatory parameters
if (!params.input_dir)            error "Missing required parameter: --input_dir"
if (!params.output_dir)           error "Missing required parameter: --output_dir"

def project_name = !params.project_name ? file(params.input_dir).toAbsolutePath().parent.name : params.project_name 

def selected_adapters = null
def fastp_opts = params.fastp_opts ?: ''
if (!params.skip_fastp && params.adapters) {
    if( !available_adapters.containsKey(params.adapters) ) {
        error "Unknown sequencing adapters: ${params.adapters}\n" +
              "Supported adapters: ${available_adapters.keySet().join(', ')}"
    }
    selected_adapters = available_adapters[params.adapters]
    fastp_opts += " ${selected_adapters.fastp_opts}"
}

def host_genome = null
def host_index = null
if (!params.skip_host_removal) {
    if (!params.host)
    {
        error "Missing required parameter: --host\n" +
              "If you do not want to perform host removal, use --skip_host_removal"
    }

    if (!available_hosts.contains(params.host) ) {
        error "Unknown host: ${params.host}\n" +
            "Supported hosts: ${available_hosts.join(', ')}"
    }
    host_index_dir = new File("${params.host_index_parent_dir}/${params.host}_default").getCanonicalPath()
    host_genome = new File(host_index_dir).getName()
    host_index = "${host_index_dir}/fasta/${host_genome}_dnaspace_index"
}

process fastp_pe {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_name), val(run_name), path(fwd_fastq_file), path(rev_fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fwd_fastq_file}"),
    path("${sample_name}/${rev_fastq_file}"), emit: fastq_files

    tuple path("${sample_name}/${run_name}.fastp_report.json"),
    path("${sample_name}/${run_name}.fastp_report.html"), emit: fastp_reports

    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*.fastp_report.{html,json}"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*${params.fastq_ext}", enabled: params.skip_host_removal

    beforeScript 'source /usr/share/modules/init/bash'
    module "fastp/${fastp_version}"

    script:
    def fwd_adapter = selected_adapters ? " --adapter_sequence    ${selected_adapters.fwd}" : ""
    def rev_adapter = selected_adapters ? " --adapter_sequence_r2 ${selected_adapters.rev}" : ""
    """
    mkdir $sample_name
    fastp \\
    --in1 ${fwd_fastq_file} \\
    --in2 ${rev_fastq_file} \\
    --cut_front --cut_tail \\
    ${fwd_adapter} \\
    ${rev_adapter} \\
    ${fastp_opts} \\
    --n_base_limit ${params.n_base_limit} \\
    --length_required ${params.length_required} \\
    --out1 ${sample_name}/${fwd_fastq_file} \\
    --out2 ${sample_name}/${rev_fastq_file} \\
    --json ${sample_name}/${run_name}.fastp_report.json \\
    --html ${sample_name}/${run_name}.fastp_report.html \\
    --thread ${task.cpus}
    """
}

process fastp_se {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_name), val(run_name), path(fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fastq_file}"), emit: fastq_file

    tuple path("${sample_name}/${run_name}.fastp_report.json"),
    path("${sample_name}/${run_name}.fastp_report.html"), emit: fastp_reports

    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*.fastp_report.{html,json}"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*${params.fastq_ext}", enabled: params.skip_host_removal

    beforeScript 'source /usr/share/modules/init/bash'
    module "fastp/${fastp_version}"

    script:
    """
    mkdir $sample_name
    fastp \\
    -i ${fastq_file} \\
    --cut_front --cut_tail \\
    ${fastp_opts} \\
    --n_base_limit ${params.n_base_limit} \\
    --length_required ${params.length_required} \\
    -o ${sample_name}/${fastq_file} \\
    --json ${sample_name}/${run_name}.fastp_report.json \\
    --html ${sample_name}/${run_name}.fastp_report.html \\
    --thread ${task.cpus}
    """
}


process host_removal_pe {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory 16.GB 

    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_name), val(run_name), path(fwd_fastq_file), path(rev_fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fwd_fastq_file}"),
    path("${sample_name}/${rev_fastq_file}"), emit: fastq_files

    tuple path("${sample_name}/${params.host_subdir}/*${params.host_reads_suffix}${params.fwd_suffix}${params.fastq_ext}"),
    path("${sample_name}/${params.host_subdir}/*${params.host_reads_suffix}${params.rev_suffix}${params.fastq_ext}"),
    optional: !params.save_host_reads, emit: host_reads_fastq_file

    publishDir "${params.output_dir}", mode: 'copy'

    beforeScript 'source /usr/share/modules/init/bash'
    module "bowtie2/${bowtie2_version}"
    module "samtools/${samtools_version}"

    script:
    def host_reads_fwd_fastq_file =
    fwd_fastq_file
        .toString()
        .replaceAll(
            "${params.fwd_suffix}${params.fastq_ext}\$",
            "${params.host_reads_suffix}${params.fwd_suffix}${params.fastq_ext}"
        )

    def host_reads_rev_fastq_file =
    rev_fastq_file
        .toString()
        .replaceAll(
            "${params.rev_suffix}${params.fastq_ext}\$",
            "${params.host_reads_suffix}${params.rev_suffix}${params.fastq_ext}"
        )

    def host_subdir = "${sample_name}/${params.host_subdir}"
    def save_host_reads_cmd = params.save_host_reads ?
    """
    mkdir ${host_subdir}

    samtools \\
    fastq \\
    -G 12 \\
    -1 ${host_subdir}/${host_reads_fwd_fastq_file} \\
    -2 ${host_subdir}/${host_reads_rev_fastq_file} \\
    -N -c 6 \\
    -@ ${task.cpus} \\
    ${run_name}.bam
    """ : ""
    """
    bowtie2 \\
    -1 ${fwd_fastq_file} \\
    -2 ${rev_fastq_file} \\
    -x ${host_index} \\
    -p ${task.cpus}  \\
    --mm | \\
    samtools view \\
    -o ${run_name}.bam \\
    -b -

    mkdir ${sample_name}
    samtools \\
    fastq \\
    -f 12 -F 256 \\
    -1 ${sample_name}/${fwd_fastq_file} \\
    -2 ${sample_name}/${rev_fastq_file} \\
    -N -c 6 \\
    -@ ${task.cpus} \\
    ${run_name}.bam

    ${save_host_reads_cmd}

    # Cleanup
    rm ${run_name}.bam
    """
}

process host_removal_se {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory 16.GB 

    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_name), val(run_name), path(fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fastq_file}"), emit: fastq_file

    path("${sample_name}/${params.host_subdir}/*${params.host_reads_suffix}${params.fastq_ext}"),
    optional: !params.save_host_reads, emit: host_reads_fastq_file

    publishDir "${params.output_dir}", mode: 'copy'

    beforeScript 'source /usr/share/modules/init/bash'
    module "bowtie2/${bowtie2_version}"
    module "samtools/${samtools_version}"

    script:
    def host_reads_fastq_file =
    fastq_file
        .toString()
        .replaceAll(
            "${params.fastq_ext}\$",
            "${params.host_reads_suffix}${params.fastq_ext}"
        )

    def host_subdir = "${sample_name}/${params.host_subdir}"
    def save_host_reads_cmd = params.save_host_reads ?
    """
    mkdir ${host_subdir}

    samtools \\
    fastq \\
    -G 4 \\
    -0 ${host_subdir}/${host_reads_fastq_file} \\
    -N -c 6 \\
    -@ ${task.cpus} \\
    ${run_name}.bam
    """ : ""
    """
    bowtie2 \\
    -U ${fastq_file} \\
    -x ${host_index} \\
    -p ${task.cpus}  \\
    --mm | \\
    samtools view \\
    -o ${run_name}.bam \\
    -b -

    mkdir ${sample_name}
    samtools \\
    fastq \\
    -f 4 -F 256 \\
    -0 ${sample_name}/${fastq_file} \\
    -N -c 6 \\
    -@ ${task.cpus} \\
    ${run_name}.bam

    ${save_host_reads_cmd}

    # Cleanup
    rm ${run_name}.bam
    """
}

process create_ini_json_files {
    tag "sample: $sample_name run: $run_name tag: $tag"

    cpus 4
    memory 8.GB

    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_name), val(run_name), path(fastq_file), val(tag),
    val(raw_read_count), val(after_fastp_read_count), val(after_host_removal_read_count)

    output:
    path("${sample_name}/*_census_stage_0.ini")
    path("${sample_name}/*_census_stage_0.json")

    publishDir "${params.output_dir}", mode: 'copy'

    script:
    def full_sample_name = fastq_file.getName().replaceAll("${params.fastq_ext}\$", "")
    def census_stage_0_ini_file = full_sample_name + "_census_stage_0.ini"
    def census_stage_0_json_file = full_sample_name + "_census_stage_0.json"
    def sequencing_device =  params.adapters ? params.adapters.split('_')[0] : "NA"
    def low_quality_read_count = raw_read_count - after_fastp_read_count
    def host_read_count = after_fastp_read_count -after_host_removal_read_count
    """
    mkdir ${sample_name}
    sed \\
    -e s/__PROJECT_NAME/${project_name}/ \\
    -e s/__SAMPLE_NAME/${sample_name}/ \\
    -e s/__FULL_SAMPLE_NAME/${full_sample_name}/ \\
    -e s/__FASTQ_FILE/${fastq_file}/ \\
    -e s/__TAG/${tag}/ \\
    -e s/__SEQUENCING_DEVICE/${sequencing_device}/ \\
    -e s/__RAW_READ_COUNT/${raw_read_count}/ \\
    -e s/__AFTER_HOST_REMOVAL_READ_COUNT/${after_host_removal_read_count}/ \\
    -e s/__HOST_READ_COUNT/${host_read_count}/ \\
    -e s/__LOW_QUALITY_READ_COUNT/${low_quality_read_count}/ \\
    ${baseDir}/templates/template_census_stage_0.ini > ${sample_name}/${census_stage_0_ini_file}

     sed \\
    -e s/__PROJECT_NAME/${project_name}/ \\
    -e s/__SAMPLE_NAME/${sample_name}/ \\
    -e s/__FULL_SAMPLE_NAME/${full_sample_name}/ \\
    -e s/__FASTQ_FILE/${fastq_file}/ \\
    -e s/__TAG/${tag}/ \\
    -e s/__SEQUENCING_DEVICE/${sequencing_device}/ \\
    -e s/__RAW_READ_COUNT/${raw_read_count}/ \\
    -e s/__AFTER_HOST_REMOVAL_READ_COUNT/${after_host_removal_read_count}/ \\
    -e s/__HOST_READ_COUNT/${host_read_count}/ \\
    -e s/__LOW_QUALITY_READ_COUNT/${low_quality_read_count}/ \\
    ${baseDir}/templates/template_census_stage_0.json > ${sample_name}/${census_stage_0_json_file}
    """
}


include { validate_pe_runs as validate_raw_pe_runs } from './include/validate_runs.nf'
include { validate_se_runs as validate_raw_se_runs } from './include/validate_runs.nf'
include { validate_pe_runs as validate_after_fastp_pe_runs } from './include/validate_runs.nf'
include { validate_se_runs as validate_after_fastp_se_runs } from '././include/validate_runs.nf'
include { validate_pe_runs as validate_after_host_removal_pe_runs } from './include/validate_runs.nf'
include { validate_se_runs as validate_after_host_removal_se_runs } from './include/validate_runs.nf'

workflow {
    def sample_dirs_ch = params.sample_list \
        ? Channel // If sample_list is provided, grab samples from list
            .fromPath(params.sample_list)
            .splitText()
            .map { line ->
                def sample_name = line.trim()
                def sample_dir = file("${params.input_dir}/${sample_name}")

                if (!sample_dir.isDirectory()) {
                      error "Sample '${sample_name}' (from ${params.sample_list}) " +
                      "was not found in input directory ${params.input_dir}"
                }

                return sample_dir
            }
        : Channel // If sample_list is not provided, use all subdirectories
            .fromPath("${params.input_dir}/*", type: 'dir')

    def sample_count_ch = sample_dirs_ch.count()
    sample_count_ch
    .subscribe { sample_count ->
        log.info "Found ${sample_count} sample${sample_count > 1 ? 's' : ''}"
    }

    raw_runs_ch = sample_dirs_ch
        .flatMap { sample_dir ->

            def sample_name = sample_dir.getName()

            def all_fastq_files = sample_dir.listFiles().findAll {
                it.name.endsWith(params.fastq_ext)
            }.sort()

            if (all_fastq_files.isEmpty()) {
                error "No FASTQ files found in directory ${sample_dir}\n" +
                "Check that your FASTQ file extension matches the value of --fastq_ext (currently '${params.fastq_ext}')"
            }

            def fwd_fastq_files = all_fastq_files.findAll { it.name.endsWith("${params.fwd_suffix}${params.fastq_ext}") }
            def rev_fastq_files = all_fastq_files.findAll { it.name.endsWith("${params.rev_suffix}${params.fastq_ext}") }
            def se_fastq_files  = all_fastq_files.findAll {
                !it.name.endsWith("${params.fwd_suffix}${params.fastq_ext}") &&
                !it.name.endsWith("${params.rev_suffix}${params.fastq_ext}")
            }

            // Extract run name from filename
            def extractRun = { fastq_file_name ->
                fastq_file_name
                    .replaceAll("${params.fwd_suffix}${params.fastq_ext}\$", "")
                    .replaceAll("${params.rev_suffix}${params.fastq_ext}\$", "")
                    .replaceAll("${params.fastq_ext}\$", "")
            }

            // Group by run name
            def runs = [:].withDefault { [fwd: null, rev: null, se: null] }

            fwd_fastq_files.each { fastq_file ->
                def run_name = extractRun(fastq_file.name)
                runs[run_name].fwd = fastq_file
            }

            rev_fastq_files.each { fastq_file ->
                def run_name = extractRun(fastq_file.name)
                runs[run_name].rev= fastq_file
            }

            se_fastq_files.each { fastq_file ->
                def run_name = extractRun(fastq_file.name)
                runs[run_name].se = fastq_file
            }

            // ----------------------------------------------------------
            // Produce one entry per run
            // ----------------------------------------------------------
            return runs.collect { run_name, run_fastq_files ->
                boolean has_pe = (run_fastq_files.fwd != null || run_fastq_files.rev != null)
                boolean has_se = run_fastq_files.se != null

                // PE run
                if (has_pe && !has_se) {
                    if (run_fastq_files.rev == null) {
                        error "Paired-end mismatch in sample '${sample_name}', run '${run_name}'\n" +
                        " Found forward file ${run_fastq_files.fwd.name} but not reverse file."
                    } else if (run_fastq_files.fwd == null) {
                        error "Paired-end mismatch in sample '${sample_name}', run '${run_name}'\n" +
                        " Found reverse file ${run_fastq_files.rev.name} but not forward file."
                    }

                    return [
                        sample_name,
                        run_name,
                        'paired_end',
                        run_fastq_files.fwd,
                        run_fastq_files.rev
                    ]
                }
                // SE run
                else if (has_se && !has_pe) {
                    return [
                        sample_name,
                        run_name,
                        'single_end',
                        run_fastq_files.se,
                        null
                    ]
                }
                // Mixed run (error)
                else if (has_pe && has_se) {
                    throw new RuntimeException(
                        "Run '${run_name}' of sample '${sample_name}' contains both PE and SE reads."
                    )
                }
            }
        }
    
    // ------------------------------------------------------------
    // Split channels
    // ------------------------------------------------------------
    raw_runs_ch = raw_runs_ch
        .branch { sample_name, run_name, library_layout, fwd_fastq_file, rev_fastq_file ->
            paired_end: library_layout == 'paired_end'
                return [sample_name, run_name, fwd_fastq_file, rev_fastq_file]

            single_end: library_layout == 'single_end'
                return [sample_name, run_name, fwd_fastq_file]
    }

    def pe_run_count_ch = raw_runs_ch.paired_end.count()
    def se_run_count_ch = raw_runs_ch.single_end.count()
    pe_run_count_ch
    .combine(se_run_count_ch)
    .subscribe { pe_run_count, se_run_count ->
        log.info "Found ${pe_run_count} paired-end run${pe_run_count > 1 ? 's' : ''}"
        log.info "Found ${se_run_count} single-end run${se_run_count > 1 ? 's' : ''}"

        if (pe_run_count == 0)
        {
            log.warn "No paired-end runs found."
            log.warn "If you expected paired-end data, please check the --fwd_suffix and --rev_suffix parameters"
        }
    }

    def extract_read_count = { read_count_file ->
        return read_count_file.readLines()[0].trim().toLong()
    }
    def extract_read_counts_pe = { pe_runs_ch ->
        pe_runs_ch.map { items ->
            def (sample_name, run_name, fwd_fastq_file, rev_fastq_file, read_count_file) = items
            def read_count = extract_read_count(read_count_file)
            [sample_name, run_name, fwd_fastq_file, rev_fastq_file, read_count]
        }
    }
    def extract_read_counts_se = { se_runs_ch ->
        se_runs_ch.map { items ->
            def (sample_name, run_name, fastq_file, read_count_file) = items
            def read_count = extract_read_count(read_count_file)
            [sample_name, run_name, fastq_file, read_count]
        }
    }

    // Check integrity of raw paired FASTQ files
    validated_raw_pe_runs_ch = validate_raw_pe_runs(raw_runs_ch.paired_end)
    raw_read_counts_pe_ch = extract_read_counts_pe(validated_raw_pe_runs_ch)
    validated_raw_pe_runs_ch = validated_raw_pe_runs_ch.map { items -> items[0..-2] }
    
    // Check integrity of raw single FASTQ files
    validated_raw_se_runs_ch = validate_raw_se_runs(raw_runs_ch.single_end)
    raw_read_counts_se_ch = extract_read_counts_se(validated_raw_se_runs_ch)
    validated_raw_se_runs_ch = validated_raw_se_runs_ch.map { items -> items[0..-2] }
   
    if (!params.skip_fastp) {
        // Run fastp on raw paired FASTQ files
        after_fastp_pe_runs_ch = fastp_pe(validated_raw_pe_runs_ch).fastq_files
        // Run fastp on raw single FASTQ files
        after_fastp_se_runs_ch = fastp_se(validated_raw_se_runs_ch).fastq_file

        // Check integrity of paired FASTQ files after fastp
        after_fastp_validated_pe_runs_ch = validate_after_fastp_pe_runs(after_fastp_pe_runs_ch)
        after_fastp_read_counts_pe_ch = extract_read_counts_pe(after_fastp_validated_pe_runs_ch)
        after_fastp_validated_pe_runs_ch = after_fastp_validated_pe_runs_ch.map { items -> items[0..-2] }

        // Check integrity of single FASTQ files after fastp
        after_fastp_validated_se_runs_ch = validate_after_fastp_se_runs(after_fastp_se_runs_ch)
        after_fastp_read_counts_se_ch = extract_read_counts_se(after_fastp_validated_se_runs_ch)
        after_fastp_validated_se_runs_ch = after_fastp_validated_se_runs_ch.map { items -> items[0..-2] }
    } else {
        after_fastp_validated_pe_runs_ch = validated_raw_pe_runs_ch
        after_fastp_read_counts_pe_ch = raw_read_counts_pe_ch

        after_fastp_validated_se_runs_ch = validated_raw_se_runs_ch
        after_fastp_read_counts_se_ch = raw_read_counts_se_ch
    }


    if (!params.skip_host_removal) {
        log.info "Selected reference genome for host ${params.host}: ${host_genome}"
        // Discard reads that map to the host genome in paired FASTQ files
        after_host_removal_pe_runs_ch = host_removal_pe(after_fastp_validated_pe_runs_ch)
        // Discard reads that map to the host genome in single FASTQ files
        after_host_removal_se_runs_ch = host_removal_se(after_fastp_validated_se_runs_ch)

        // Check integrity of paired FASTQ files after host removal
        after_host_removal_validated_pe_runs_ch =  validate_after_host_removal_pe_runs(after_host_removal_pe_runs_ch.fastq_files)
        after_host_removal_read_counts_pe_ch = extract_read_counts_pe(after_host_removal_validated_pe_runs_ch)
        after_host_removal_validated_pe_runs_ch = after_host_removal_validated_pe_runs_ch.map { items -> items[0..-2] }

        // Check integrity of single FASTQ files after host removal
        after_host_removal_validated_se_runs_ch =  validate_after_host_removal_se_runs(after_host_removal_se_runs_ch.fastq_file)
        after_host_removal_read_counts_se_ch = extract_read_counts_se(after_host_removal_validated_se_runs_ch)
        after_host_removal_validated_se_runs_ch = after_host_removal_validated_se_runs_ch.map { items -> items[0..-2] }
    } else {
        after_host_removal_read_counts_pe_ch = after_fastp_read_counts_pe_ch
        after_host_removal_read_counts_se_ch = after_fastp_read_counts_se_ch
    }
    // Create a channel that indicates for each FASTQ file
    // 1. the number of raw reads
    // 2. the number of reads after fastp
    // 3. the number of reads after fastp and host removal
    combined_read_counts_pe_ch = raw_read_counts_pe_ch
        .join(after_fastp_read_counts_pe_ch, by: [0,1])
        .join(after_host_removal_read_counts_pe_ch, by: [0,1])
        .flatMap { items -> 
            def sample_name = items[0]
            def run_name = items[1]
            def fwd_fastq_file = items[8]
            def rev_fastq_file = items[9]
            def raw_read_count = items[4]
            def after_fastp_read_count = items[7]
            def after_host_removal_read_count = items[10]
            
            return [
                [sample_name, run_name, fwd_fastq_file, "1",
                raw_read_count, after_fastp_read_count, after_host_removal_read_count],
                [sample_name, run_name, rev_fastq_file, "2",
                raw_read_count, after_fastp_read_count, after_host_removal_read_count]
            ]
        }
    
    combined_read_counts_se_ch = raw_read_counts_se_ch
        .join(after_fastp_read_counts_se_ch, by: [0,1])
        .join(after_host_removal_read_counts_se_ch, by: [0,1])
        .map { items ->
            def sample_name = items[0]
            def run_name = items[1]
            def fastq_file = items[6]
            def raw_read_count = items[3]
            def after_fastp_read_count = items[5]
            def after_host_removal_read_count = items[7]

            [sample_name, run_name, fastq_file, "single",
            raw_read_count, after_fastp_read_count, after_host_removal_read_count]
    }
     
    combined_read_counts_ch = combined_read_counts_pe_ch.mix(combined_read_counts_se_ch)

    create_ini_json_files(combined_read_counts_ch)
}

workflow.onComplete {
    def datetime_formatter = DateTimeFormatter.ISO_OFFSET_DATE_TIME

    def global_params = [
        project_name: project_name,
        input_dir: file(params.input_dir).toAbsolutePath().toString(),
        output_dir: file(params.output_dir).toAbsolutePath().toString()
    ]
    global_params += params.subMap(['sample_list', 'fastq_ext', 'fwd_suffix', 'rev_suffix'])

    def fastp_params = [ skip_fastp: params.skip_fastp ]
    if (!params.skip_fastp) {
        fastp_params['fastp_version'] = fastp_version
        fastp_params += params.subMap(['adapters', 'n_base_limit', 'length_required'])
        if (params.adapters)
        {
            fastp_params['adapters_fwd'] = selected_adapters.fwd
            fastp_params['adapters_rev'] = selected_adapters.rev
        }
        fastp_params['fastp_opts'] = fastp_opts
    }

    def host_removal_params = [ skip_host_removal: params.skip_host_removal ]
    if (!params.skip_host_removal) {
        host_removal_params += [
            bowtie2_version: bowtie2_version,
            samtools_version: samtools_version,
            host: params.host,
            host_index_dir: host_index_dir,
            save_host_reads: params.save_host_reads
        ]

        if (params.save_host_reads) {
            host_removal_params += params.subMap(['host_reads_suffix', 'host_subdir'])
        }
    }

    def metadata = [
        workflow: [
            name: workflow.scriptName,
            version: "git describe --tags".execute().text.trim()
        ],
        execution_info: [
            username: System.getProperty("user.name"),
            start_time: workflow.start.format(datetime_formatter),
            end_time: workflow.complete.format(datetime_formatter),
            duration: workflow.duration.toString()
        ],
        workflow_parameters: [
            global: global_params,
            fastp: fastp_params,
            host_removal: host_removal_params
        ]
    ]

    // Définition du chemin du fichier de sortie
    def reportFile = file("${params.output_dir}/fastq_qc_workflow_summary.json")

    // Création du répertoire parent si nécessaire et écriture du fichier
    reportFile.parent.mkdirs()
    reportFile.text = JsonOutput.prettyPrint(JsonOutput.toJson(metadata))

    log.info "Workflow metadata saved to ${reportFile}"
}
