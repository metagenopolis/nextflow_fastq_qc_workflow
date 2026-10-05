process VALIDATE_PE_RUNS {
    tag "sample: $sample_name run: $run_name"

    cpus 4
    memory 8.GB

    errorStrategy 'retry'
    maxRetries 1

    // Publish only when nothing runs after this validation, i.e. when both fastp
    // and host removal are skipped. Otherwise the last stage that runs publishes.
    publishDir "${params.output_dir}", mode: 'copy', pattern: "*${params.fastq_ext}", enabled: params.skip_host_removal && params.skip_fastp, saveAs: { file_name -> "${sample_name}/${file_name}" }

    beforeScript 'source /usr/share/modules/init/bash'
    module 'seqtk/1.4'

    input:
    tuple val(sample_name), val(run_name), path(fwd_fastq_file), path(rev_fastq_file)

    output:
    tuple val(sample_name), val(run_name), path(fwd_fastq_file), path(rev_fastq_file),
    path("${run_name}.read_count.tsv")

    script:
    """
    # Check files integrity
    gzip -tvf ${fwd_fastq_file}
    gzip -tvf ${rev_fastq_file}

    FWD_READ_COUNT=\$(seqtk size ${fwd_fastq_file} | cut -f1)
    REV_READ_COUNT=\$(seqtk size ${rev_fastq_file} | cut -f1)

    # Check if FASTQ files are empty
    if [ "\${FWD_READ_COUNT}" -eq 0 ] ; then
        echo "Error: Forward FASTQ file contains zero reads." >&2
        echo "sample: ${sample_name} run: ${run_name} file: ${fwd_fastq_file}" >&2
        exit 1
    fi
    if [ "\${REV_READ_COUNT}" -eq 0 ] ; then
        echo "Error: Reverse FASTQ file contains zero reads." >&2
        echo "sample: ${sample_name} run: ${run_name} file: ${rev_fastq_file}" >&2
        exit 1
    fi

    # Check if FASTQ files have the same number of reads
    if [ "\${FWD_READ_COUNT}" -ne "\${REV_READ_COUNT}"  ] ; then
        echo "ERROR: Paired-end read count mismatch." >&2
        echo "sample: ${sample_name} run: ${run_name}" >&2
        echo "Forward read count: \$FWD_READ_COUNT  (file: ${fwd_fastq_file})" >&2
        echo "Reverse read count: \$REV_READ_COUNT  (file: ${rev_fastq_file})" >&2 
        exit 1
    fi

    echo \${FWD_READ_COUNT} > ${run_name}.read_count.tsv
    """
}

process VALIDATE_SE_RUNS {
    tag "sample: $sample_name run: $run_name"

    cpus 4
    memory 8.GB

    errorStrategy 'retry'
    maxRetries 1

    // Publish only when nothing runs after this validation, i.e. when both fastp
    // and host removal are skipped. Otherwise the last stage that runs publishes.
    publishDir "${params.output_dir}", mode: 'copy', pattern: "*${params.fastq_ext}", enabled: params.skip_host_removal && params.skip_fastp, saveAs: { file_name -> "${sample_name}/${file_name}" }

    beforeScript 'source /usr/share/modules/init/bash'
    module 'seqtk/1.4'

    input:
    tuple val(sample_name), val(run_name), path(fastq_file)

    output:
    tuple val(sample_name), val(run_name), path(fastq_file),
    path("${run_name}.read_count.tsv")

    script:
    """
    # Check file integrity
    gzip -tvf ${fastq_file}

    READ_COUNT=\$(seqtk size ${fastq_file} | cut -f1)

    # Check if FASTQ file is empty
    if [ "\${READ_COUNT}" -eq 0 ] ; then
        echo "Error: FASTQ file contains zero reads." >&2
        echo "sample: ${sample_name} run: ${run_name} file: ${fastq_file}" >&2
        exit 1
    fi

    echo \${READ_COUNT} > ${run_name}.read_count.tsv
    """
}