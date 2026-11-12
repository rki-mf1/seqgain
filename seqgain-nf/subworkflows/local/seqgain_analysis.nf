include { SEQGAIN } from '../../modules/local/seqgain/main'

workflow SEQGAIN_ANALYSIS {
    take:
    ch_input

    main:
    SEQGAIN(ch_input)
}
