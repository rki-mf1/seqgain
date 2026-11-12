/* Pair trimmed reads with the corresponding mapped, sorted BAM and index. */
include { SEQGAIN } from '../../../modules/local/seqgain/main'

workflow SEQGAIN_ANALYSIS {
    take:
    ch_reads_trimmed
    ch_mapped_sorted_bam_bai

    main:
    ch_seqgain_input = ch_reads_trimmed
        .map { meta, reads -> tuple(meta.id, meta, reads) }
        .join(
            ch_mapped_sorted_bam_bai.map { meta, bam, bai ->
                tuple(meta.id, meta, bam, bai)
            },
            by: 0
        )
        .map { sample_id, reads_meta, reads, bam_meta, bam, bai ->
            tuple(reads_meta, reads, bam, bai)
        }

    SEQGAIN(ch_seqgain_input)
}
