process SEQGAIN {
    tag "${meta.id}"
    publishDir "${params.outdir}/seqgain", mode: 'copy'
    conda 'bioconda::seqgain=0.3.1'

    input:
    tuple val(meta), path(reads), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    def args = task.ext.args ?: ''
    def hasBam = bam && bam.toString() != '[]'
    def readArgs = meta.mapped_only ? '' : "--reads1 ${reads[0]}" + (meta.single_end ? '' : " --reads2 ${reads[1]}")
    def bamArg = hasBam ? "--bam ${bam}" : ''
    def sourceArg = hasBam && !meta.mapped_only ? '--kmer-source fastq' : ''
    """
    seqgain \\
        ${readArgs} \\
        ${bamArg} \\
        ${sourceArg} \\
        --sample '${meta.id}' \\
        --outdir ${meta.id} \\
        ${args}
    """
}
