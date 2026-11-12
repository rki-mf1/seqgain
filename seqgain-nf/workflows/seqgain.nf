process SAMTOOLS_MAPPED_ONLY {
    tag "${meta.id}"
    conda 'bioconda::samtools=1.21'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.mapped.bam"), path("${meta.id}.mapped.bam.bai"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -b -F 4 -o ${meta.id}.mapped.bam ${bam}
    samtools index -@ ${task.cpus} ${meta.id}.mapped.bam
    """
}

include { FASTP } from '../modules/nf-core/fastp/main'
include { BOWTIE2_BUILD } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SEQGAIN_ANALYSIS } from '../subworkflows/local/seqgain_analysis'

workflow SEQGAIN_WORKFLOW {
    take:
    ch_samples

    main:
    ch_reads = ch_samples.map { meta, reads, ref -> tuple(meta, reads, []) }
    FASTP(ch_reads, false, false, false)

    ch_ref_by_id = ch_samples.filter { meta, reads, ref -> ref != null }
        .map { meta, reads, ref -> tuple(meta.id, meta, ref) }
    ch_trimmed_by_id = FASTP.out.reads.map { meta, reads -> tuple(meta.id, meta, reads) }

    ch_with_ref = ch_trimmed_by_id.join(ch_ref_by_id, by: 0)
    ch_without_ref = ch_trimmed_by_id.join(
        ch_samples.filter { meta, reads, ref -> ref == null }
            .map { meta, reads, ref -> tuple(meta.id) }, by: 0
    ).map { id, meta, reads -> tuple(meta, reads, [], []) }

    ch_reference = ch_with_ref.map { id, meta, reads, ref -> tuple(meta, ref) }
    BOWTIE2_BUILD(ch_reference)

    ch_index_keyed = BOWTIE2_BUILD.out.index.map { meta, index -> tuple(meta.id, meta, index) }
    ch_align = ch_with_ref.join(ch_index_keyed, by: 0)
        .map { id, meta, reads, ref, idxmeta, index -> tuple(meta, reads, index, ref) }
        .multiMap { meta, reads, index, ref ->
            reads: tuple(meta, reads)
            index: tuple(meta, index)
            fasta: tuple(meta, ref)
        }
    BOWTIE2_ALIGN(ch_align.reads, ch_align.index, ch_align.fasta, false, true)
    SAMTOOLS_INDEX(BOWTIE2_ALIGN.out.bam)

    ch_mapped_bam = BOWTIE2_ALIGN.out.bam
        .map { meta, bam -> tuple(meta.id, meta, bam) }
        .join(SAMTOOLS_INDEX.out.index.map { meta, bai -> tuple(meta.id, bai) }, by: 0)
    ch_mapped = ch_trimmed_by_id.join(ch_mapped_bam, by: 0)
        .map { id, meta, reads, bammeta, bam, bai -> tuple(meta, reads, bam, bai) }

    if (params.mapped_only) {
        SAMTOOLS_MAPPED_ONLY(ch_mapped.map { meta, reads, bam, bai -> tuple(meta, bam, bai) })
        ch_seqgain_input = SAMTOOLS_MAPPED_ONLY.out.bam.map { meta, bam, bai -> tuple(meta, [], bam, bai) }
    } else {
        ch_seqgain_input = ch_mapped.mix(ch_without_ref)
    }
    SEQGAIN_ANALYSIS(ch_seqgain_input)
}
