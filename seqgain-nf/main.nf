#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SEQGAIN_WORKFLOW } from './workflows/seqgain'

params.input = params.input ?: null
params.outdir = params.outdir ?: 'results'

workflow {
    if (!params.input) error 'Provide --input samplesheet.csv'
    if (!(params.mapped_only in [true, false])) error '--mapped_only must be true or false'

    ch_samples = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def id = row.sampleID?.toString()?.trim()
            def r1 = row.read1?.toString()?.trim()
            def r2 = row.read2?.toString()?.trim()
            def ref = row.reference?.toString()?.trim()
            if (!id || !r1) error "Sample ${id ?: '<missing>'}: sampleID and read1 are required"
            if (params.mapped_only && !ref) error "Sample ${id}: --mapped_only requires a reference"
            if (ref && !file(ref).exists()) error "Sample ${id}: reference not found: ${ref}"
            if (!file(r1).exists()) error "Sample ${id}: read1 not found: ${r1}"
            if (r2 && !file(r2).exists()) error "Sample ${id}: read2 not found: ${r2}"
            def meta = [id: id, single_end: !r2, mapped_only: params.mapped_only]
            tuple(meta, r2 ? [file(r1), file(r2)] : [file(r1)], ref ? file(ref) : null)
        }

    SEQGAIN_WORKFLOW(ch_samples)
}
