params.fastp_version = '1.3.1'
params.bowtie2_version = '2.5.5'
params.samtools_version = '1.23.1'


params.skip_fastp = false
params.fastq_ext  = ".fastq.gz"
params.fwd_suffix = "_1"
params.rev_suffix = "_2"
params.n_base_limit = 0
params.length_required = 60
params.adapters = null
params.fastp_opts = ""

params.skip_host_removal = false
params.host = null
params.host_index_parent_dir = '/biodatabank/meteor-data/reference.ng'
params.save_host_reads = false
params.host_reads_suffix = '_Host'
params.host_subdir= 'host'

params.sample_list = null
params.project_name = null
params.help = false

def available_adapters() {
    [
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
}


def available_hosts() {
    [
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
}


def print_help() {
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
                                      Supported adapters: ${available_adapters().keySet().join(', ')}
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
                                      Supported hosts: ${available_hosts().join(', ')}
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
}


def project_name() {
    params.project_name ?: file(params.input_dir).toAbsolutePath().parent.name
}


def selected_adapters() {
    if (params.skip_fastp || !params.adapters) {
        return null
    }

    if (!available_adapters().containsKey(params.adapters)) {
        error "Unknown sequencing adapters: ${params.adapters}\n" +
              "Supported adapters: ${available_adapters().keySet().join(', ')}"
    }

    available_adapters()[params.adapters]
}


def fastp_opts() {
    def opts = params.fastp_opts ?: ''
    def adapters = selected_adapters()

    if (adapters?.fastp_opts) {
        opts += " ${adapters.fastp_opts}"
    }

    opts
}


def host_index_dir() {
    new File("${params.host_index_parent_dir}/${params.host}_default").getCanonicalPath()
}


def host_genome() {
    new File(host_index_dir()).getName()
}


def host_index() {
    "${host_index_dir()}/fasta/${host_genome()}_dnaspace_index"
}


// Check the mandatory parameters before any task is launched.
def validate_params() {
    if (!params.input_dir) {
        error "Missing required parameter: --input_dir"
    }

    if (!params.output_dir) {
        error "Missing required parameter: --output_dir"
    }

    // Rejects an unknown adapters name.
    selected_adapters()

    if (!params.skip_host_removal) {
        if (!params.host) {
            error "Missing required parameter: --host\n" +
                  "If you do not want to perform host removal, use --skip_host_removal"
        }

        if (!available_hosts().contains(params.host)) {
            error "Unknown host: ${params.host}\n" +
                  "Supported hosts: ${available_hosts().join(', ')}"
        }
    }
}


process FASTP_PE {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*.fastp_report.{html,json}"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*${params.fastq_ext}", enabled: params.skip_host_removal

    beforeScript 'source /usr/share/modules/init/bash'
    module "fastp/${params.fastp_version}"

    input:
    tuple val(sample_name), val(run_name), path(fwd_fastq_file), path(rev_fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fwd_fastq_file}"),
    path("${sample_name}/${rev_fastq_file}"), emit: fastq_files

    tuple path("${sample_name}/${run_name}.fastp_report.json"),
    path("${sample_name}/${run_name}.fastp_report.html"), emit: fastp_reports

    script:
    def adapters = selected_adapters()
    def fwd_adapter = adapters ? " --adapter_sequence    ${adapters.fwd}" : ""
    def rev_adapter = adapters ? " --adapter_sequence_r2 ${adapters.rev}" : ""
    """
    mkdir $sample_name
    fastp \\
    --in1 ${fwd_fastq_file} \\
    --in2 ${rev_fastq_file} \\
    --cut_front --cut_tail \\
    ${fwd_adapter} \\
    ${rev_adapter} \\
    ${fastp_opts()} \\
    --n_base_limit ${params.n_base_limit} \\
    --length_required ${params.length_required} \\
    --out1 ${sample_name}/${fwd_fastq_file} \\
    --out2 ${sample_name}/${rev_fastq_file} \\
    --json ${sample_name}/${run_name}.fastp_report.json \\
    --html ${sample_name}/${run_name}.fastp_report.html \\
    --thread ${task.cpus}
    """
}

process FASTP_SE {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*.fastp_report.{html,json}"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "${sample_name}/*${params.fastq_ext}", enabled: params.skip_host_removal

    beforeScript 'source /usr/share/modules/init/bash'
    module "fastp/${params.fastp_version}"

    input:
    tuple val(sample_name), val(run_name), path(fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fastq_file}"), emit: fastq_file

    tuple path("${sample_name}/${run_name}.fastp_report.json"),
    path("${sample_name}/${run_name}.fastp_report.html"), emit: fastp_reports

    script:
    """
    mkdir $sample_name
    fastp \\
    -i ${fastq_file} \\
    --cut_front --cut_tail \\
    ${fastp_opts()} \\
    --n_base_limit ${params.n_base_limit} \\
    --length_required ${params.length_required} \\
    -o ${sample_name}/${fastq_file} \\
    --json ${sample_name}/${run_name}.fastp_report.json \\
    --html ${sample_name}/${run_name}.fastp_report.html \\
    --thread ${task.cpus}
    """
}


process HOST_REMOVAL_PE {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory 16.GB 

    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.output_dir}", mode: 'copy'

    beforeScript 'source /usr/share/modules/init/bash'
    module "bowtie2/${params.bowtie2_version}"
    module "samtools/${params.samtools_version}"

    input:
    tuple val(sample_name), val(run_name), path(fwd_fastq_file), path(rev_fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fwd_fastq_file}"),
    path("${sample_name}/${rev_fastq_file}"), emit: fastq_files

    tuple path("${sample_name}/${params.host_subdir}/*${params.host_reads_suffix}${params.fwd_suffix}${params.fastq_ext}"),
    path("${sample_name}/${params.host_subdir}/*${params.host_reads_suffix}${params.rev_suffix}${params.fastq_ext}"),
    optional: !params.save_host_reads, emit: host_reads_fastq_file

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
    set -o pipefail

    bowtie2 \\
    -1 ${fwd_fastq_file} \\
    -2 ${rev_fastq_file} \\
    -x ${host_index()} \\
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

process HOST_REMOVAL_SE {
    tag "sample: $sample_name run: $run_name"

    cpus 8
    memory 16.GB 

    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.output_dir}", mode: 'copy'

    beforeScript 'source /usr/share/modules/init/bash'
    module "bowtie2/${params.bowtie2_version}"
    module "samtools/${params.samtools_version}"

    input:
    tuple val(sample_name), val(run_name), path(fastq_file)

    output:
    tuple val(sample_name), val(run_name),
    path("${sample_name}/${fastq_file}"), emit: fastq_file

    path("${sample_name}/${params.host_subdir}/*${params.host_reads_suffix}${params.fastq_ext}"),
    optional: !params.save_host_reads, emit: host_reads_fastq_file

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
    set -o pipefail

    bowtie2 \\
    -U ${fastq_file} \\
    -x ${host_index()} \\
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

process CREATE_INI_JSON_FILES {
    tag "sample: $sample_name run: $run_name tag: $tag"

    cpus 4
    memory 8.GB

    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(sample_name), val(run_name), path(fastq_file), val(tag),
    val(raw_read_count), val(after_fastp_read_count), val(after_host_removal_read_count)

    output:
    path("${sample_name}/*_census_stage_0.ini")
    path("${sample_name}/*_census_stage_0.json")

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
    -e s/__PROJECT_NAME/${project_name()}/ \\
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
    -e s/__PROJECT_NAME/${project_name()}/ \\
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


include { VALIDATE_PE_RUNS as VALIDATE_RAW_PE_RUNS } from './include/validate_runs.nf'
include { VALIDATE_SE_RUNS as VALIDATE_RAW_SE_RUNS } from './include/validate_runs.nf'
include { VALIDATE_PE_RUNS as VALIDATE_AFTER_FASTP_PE_RUNS } from './include/validate_runs.nf'
include { VALIDATE_SE_RUNS as VALIDATE_AFTER_FASTP_SE_RUNS } from '././include/validate_runs.nf'
include { VALIDATE_PE_RUNS as VALIDATE_AFTER_HOST_REMOVAL_PE_RUNS } from './include/validate_runs.nf'
include { VALIDATE_SE_RUNS as VALIDATE_AFTER_HOST_REMOVAL_SE_RUNS } from './include/validate_runs.nf'

// Remove `suffix` from the end of `name`. Only ever called on a name already
// known to end with the suffix, so a plain substring is exact.
def strip_suffix(String name, String suffix) {
    name.substring(0, name.length() - suffix.length())
}


// Return the sample directories to process, either read from the sample list
// file or discovered as the subdirectories of the input directory.
def resolve_sample_dirs() {
    def input_dir = file(params.input_dir)

    if (!input_dir.isDirectory()) {
        error "Input directory does not exist, or is not a directory: ${input_dir}"
    }

    if (params.sample_list) {
        def sample_names = file(params.sample_list).readLines()
            .collect { line -> line.trim() }
            .findAll { line -> line && !line.startsWith('#') }

        if (!sample_names) {
            error "No sample name found in ${params.sample_list}"
        }

        return sample_names.collect { sample_name -> file("${input_dir}/${sample_name}") }
    }

    // listFiles() does not descend into a symlink to a directory, it returns the
    // link itself, so always list the target the link resolves to.
    def sample_dirs = file(input_dir.toRealPath()).listFiles()
        .findAll { entry -> entry.isDirectory() && !entry.name.startsWith('.') }
        .sort()

    if (!sample_dirs) {
        error "No sample directory found in ${input_dir}"
    }

    return sample_dirs
}


// Scan every sample directory and group its FASTQ files by run. Every problem
// found is collected so that a broken input is reported in full, in one go,
// before any task is launched.
def parse_runs() {
    def fwd_suffix = "${params.fwd_suffix}${params.fastq_ext}"
    def rev_suffix = "${params.rev_suffix}${params.fastq_ext}"

    def pe_runs  = []
    def se_runs  = []
    def problems = []

    resolve_sample_dirs()
        .each { sample_dir ->
            def sample_name = sample_dir.getName()

            // Collect this as a problem rather than letting listFiles() throw,
            // so that a stale path does not hide the other samples' problems.
            if (!sample_dir.isDirectory()) {
                problems << (params.sample_list
                    ? "Sample '${sample_name}' (from ${params.sample_list}) " +
                      "was not found in input directory ${params.input_dir}"
                    : "Sample directory does not exist, or is not a directory:\n" +
                      "  - ${sample_dir}")
                return
            }

            // Index the FASTQ files by run name, i.e. the file name with its mate
            // suffix removed, so that mates are matched by name rather than by
            // position in two independently sorted lists.
            def fwd_by_run = [:]
            def rev_by_run = [:]
            def se_by_run  = [:]

            // A symlinked sample directory has to be listed through its target, but
            // the sample keeps the name it was given here rather than the target's.
            file(sample_dir.toRealPath()).listFiles().each { fastq_file ->
                if (fastq_file.name.endsWith(fwd_suffix)) {
                    fwd_by_run[strip_suffix(fastq_file.name, fwd_suffix)] = fastq_file
                }
                else if (fastq_file.name.endsWith(rev_suffix)) {
                    rev_by_run[strip_suffix(fastq_file.name, rev_suffix)] = fastq_file
                }
                else if (fastq_file.name.endsWith(params.fastq_ext)) {
                    se_by_run[strip_suffix(fastq_file.name, params.fastq_ext)] = fastq_file
                }
            }

            def run_names = (fwd_by_run.keySet() + rev_by_run.keySet() + se_by_run.keySet()).sort()

            if (!run_names) {
                problems << "No FASTQ files found in directory ${sample_dir}\n" +
                    "Check that your FASTQ file extension matches the value of --fastq_ext " +
                    "(currently '${params.fastq_ext}')"
                return
            }

            // A run is paired-end as soon as one mate suffix matched, and it then
            // owes both mates. The two layouts cannot be mixed within a run.
            def mixed = run_names.findAll { run_name ->
                (fwd_by_run[run_name] || rev_by_run[run_name]) && se_by_run[run_name]
            }

            def orphans = run_names.findAll { run_name ->
                !se_by_run[run_name] && (!fwd_by_run[run_name] || !rev_by_run[run_name])
            }

            def sample_problems = []

            if (mixed) {
                sample_problems << "Paired-end and single-end FASTQ files for the same run " +
                    "of sample '${sample_name}':\n" +
                    mixed.collect { run_name ->
                        def fastq_files = [fwd_by_run[run_name], rev_by_run[run_name], se_by_run[run_name]]
                        "  - run '${run_name}': " +
                            fastq_files.findAll().collect { fastq_file -> fastq_file.name }.join(', ')
                    }.join("\n") + "\n" +
                    "If you expected only paired-end data, please check the --fwd_suffix and --rev_suffix parameters."
            }

            if (orphans) {
                sample_problems << "Unpaired FASTQ files for sample '${sample_name}':\n" +
                    orphans.collect { run_name ->
                        def present = fwd_by_run[run_name] ?: rev_by_run[run_name]
                        def missing = fwd_by_run[run_name] ? rev_suffix : fwd_suffix
                        "  - ${present.name} has no mate (expected ${run_name}${missing})"
                    }.join("\n")
            }

            if (sample_problems) {
                problems << sample_problems.join("\n\n")
                return
            }

            run_names.each { run_name ->
                if (se_by_run[run_name]) {
                    se_runs << [sample_name, run_name, se_by_run[run_name]]
                }
                else {
                    pe_runs << [sample_name, run_name, fwd_by_run[run_name], rev_by_run[run_name]]
                }
            }
        }

    if (problems) {
        error("Found ${problems.size()} invalid sample(s) in '${params.sample_list ?: params.input_dir}':\n\n" +
            problems.join("\n\n"))
    }

    return [paired_end: pe_runs, single_end: se_runs]
}


// Read the number of reads counted by the validation processes.
def extract_read_count(read_count_file) {
    read_count_file.readLines()[0].trim().toLong()
}


// Replace the read count file of each run by the read count it holds.
def extract_read_counts_pe(pe_runs_ch) {
    pe_runs_ch.map { items ->
        def (sample_name, run_name, fwd_fastq_file, rev_fastq_file, read_count_file) = items
        def read_count = extract_read_count(read_count_file)
        [sample_name, run_name, fwd_fastq_file, rev_fastq_file, read_count]
    }
}


def extract_read_counts_se(se_runs_ch) {
    se_runs_ch.map { items ->
        def (sample_name, run_name, fastq_file, read_count_file) = items
        def read_count = extract_read_count(read_count_file)
        [sample_name, run_name, fastq_file, read_count]
    }
}


// Version of the workflow, as `git describe` reports it. The -C option makes git
// describe the repository holding the workflow, not the one the run happens to be
// launched from. Anything that stops git from answering, such as a copy that is not
// a checkout or a checkout without any tag, is reported as an unknown version.
def workflow_version() {
    try {
        def git = ['git', '-C', "${projectDir}", 'describe', '--tags'].execute()
        def version = git.text.trim()
        git.waitFor()

        return git.exitValue() == 0 && version ? version : 'unknown'
    }
    catch (IOException _e) {
        return 'unknown'
    }
}


// Write a JSON report describing the workflow execution.
def write_workflow_summary() {
    def datetime_formatter = java.time.format.DateTimeFormatter.ISO_OFFSET_DATE_TIME

    def global_params = [
        project_name: project_name(),
        input_dir: file(params.input_dir).toAbsolutePath().toString(),
        output_dir: file(params.output_dir).toAbsolutePath().toString()
    ]
    global_params += params.subMap(['sample_list', 'fastq_ext', 'fwd_suffix', 'rev_suffix'])

    def fastp_params = [ skip_fastp: params.skip_fastp ]
    if (!params.skip_fastp) {
        fastp_params['fastp_version'] = params.fastp_version
        fastp_params += params.subMap(['adapters', 'n_base_limit', 'length_required'])
        if (params.adapters)
        {
            def adapters = selected_adapters()
            fastp_params['adapters_fwd'] = adapters.fwd
            fastp_params['adapters_rev'] = adapters.rev
        }
        fastp_params['fastp_opts'] = fastp_opts()
    }

    def host_removal_params = [ skip_host_removal: params.skip_host_removal ]
    if (!params.skip_host_removal) {
        host_removal_params += [
            bowtie2_version: params.bowtie2_version,
            samtools_version: params.samtools_version,
            host: params.host,
            host_index_dir: host_index_dir(),
            save_host_reads: params.save_host_reads
        ]

        if (params.save_host_reads) {
            host_removal_params += params.subMap(['host_reads_suffix', 'host_subdir'])
        }
    }

    def metadata = [
        workflow: [
            name: workflow.scriptName,
            version: workflow_version()
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
    reportFile.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(metadata))

    log.info "Workflow metadata saved to ${reportFile}"
}


workflow {
    if (params.help || !workflow.commandLine.contains('--')) {
        print_help()
        return
    }

    validate_params()

    def runs = parse_runs()

    def sample_names = (runs.paired_end + runs.single_end).collect { run -> run[0] }.unique()
    def pe_run_count = runs.paired_end.size()
    def se_run_count = runs.single_end.size()

    log.info "Found ${sample_names.size()} sample${sample_names.size() == 1 ? '' : 's'}"
    log.info "Found ${pe_run_count} paired-end run${pe_run_count == 1 ? '' : 's'}"
    log.info "Found ${se_run_count} single-end run${se_run_count == 1 ? '' : 's'}"

    if (pe_run_count == 0) {
        log.warn "No paired-end runs found."
        log.warn "If you expected paired-end data, please check the --fwd_suffix and --rev_suffix parameters"
    }

    def raw_pe_runs_ch = channel.fromList(runs.paired_end)
    def raw_se_runs_ch = channel.fromList(runs.single_end)


    // Check integrity of raw paired FASTQ files
    validated_raw_pe_runs_ch = VALIDATE_RAW_PE_RUNS(raw_pe_runs_ch)
    raw_read_counts_pe_ch = extract_read_counts_pe(validated_raw_pe_runs_ch)
    validated_raw_pe_runs_ch = validated_raw_pe_runs_ch.map { items -> items[0..-2] }
    
    // Check integrity of raw single FASTQ files
    validated_raw_se_runs_ch = VALIDATE_RAW_SE_RUNS(raw_se_runs_ch)
    raw_read_counts_se_ch = extract_read_counts_se(validated_raw_se_runs_ch)
    validated_raw_se_runs_ch = validated_raw_se_runs_ch.map { items -> items[0..-2] }
   
    if (!params.skip_fastp) {
        // Run fastp on raw paired FASTQ files
        after_fastp_pe_runs_ch = FASTP_PE(validated_raw_pe_runs_ch).fastq_files
        // Run fastp on raw single FASTQ files
        after_fastp_se_runs_ch = FASTP_SE(validated_raw_se_runs_ch).fastq_file

        // Check integrity of paired FASTQ files after fastp
        after_fastp_validated_pe_runs_ch = VALIDATE_AFTER_FASTP_PE_RUNS(after_fastp_pe_runs_ch)
        after_fastp_read_counts_pe_ch = extract_read_counts_pe(after_fastp_validated_pe_runs_ch)
        after_fastp_validated_pe_runs_ch = after_fastp_validated_pe_runs_ch.map { items -> items[0..-2] }

        // Check integrity of single FASTQ files after fastp
        after_fastp_validated_se_runs_ch = VALIDATE_AFTER_FASTP_SE_RUNS(after_fastp_se_runs_ch)
        after_fastp_read_counts_se_ch = extract_read_counts_se(after_fastp_validated_se_runs_ch)
        after_fastp_validated_se_runs_ch = after_fastp_validated_se_runs_ch.map { items -> items[0..-2] }
    } else {
        after_fastp_validated_pe_runs_ch = validated_raw_pe_runs_ch
        after_fastp_read_counts_pe_ch = raw_read_counts_pe_ch

        after_fastp_validated_se_runs_ch = validated_raw_se_runs_ch
        after_fastp_read_counts_se_ch = raw_read_counts_se_ch
    }


    if (!params.skip_host_removal) {
        log.info "Selected reference genome for host ${params.host}: ${host_genome()}"
        // Discard reads that map to the host genome in paired FASTQ files
        after_host_removal_pe_runs_ch = HOST_REMOVAL_PE(after_fastp_validated_pe_runs_ch)
        // Discard reads that map to the host genome in single FASTQ files
        after_host_removal_se_runs_ch = HOST_REMOVAL_SE(after_fastp_validated_se_runs_ch)

        // Check integrity of paired FASTQ files after host removal
        after_host_removal_validated_pe_runs_ch =  VALIDATE_AFTER_HOST_REMOVAL_PE_RUNS(after_host_removal_pe_runs_ch.fastq_files)
        after_host_removal_read_counts_pe_ch = extract_read_counts_pe(after_host_removal_validated_pe_runs_ch)
        after_host_removal_validated_pe_runs_ch = after_host_removal_validated_pe_runs_ch.map { items -> items[0..-2] }

        // Check integrity of single FASTQ files after host removal
        after_host_removal_validated_se_runs_ch =  VALIDATE_AFTER_HOST_REMOVAL_SE_RUNS(after_host_removal_se_runs_ch.fastq_file)
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

    CREATE_INI_JSON_FILES(combined_read_counts_ch)

    workflow.onComplete {
        write_workflow_summary()
    }
}
