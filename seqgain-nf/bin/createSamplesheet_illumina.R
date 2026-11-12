#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
options <- list(
  make_option(c('-f', '--folder'), type='character', help='Directory containing FASTQ files'),
  make_option(c('-o', '--outfile'), type='character', help='Output CSV'),
  make_option(c('-r', '--reffile'), type='character', default='', help='Optional reference FASTA for all samples')
)
opt <- parse_args(OptionParser(option_list=options))
if (is.null(opt$folder) || is.null(opt$outfile)) stop('Provide --folder and --outfile')
if (!dir.exists(opt$folder)) stop('FASTQ folder does not exist')
if (nzchar(opt$reffile) && !file.exists(opt$reffile)) stop('Reference FASTA does not exist')
fastqs <- sort(list.files(opt$folder, pattern='\\.(fastq|fq)(\\.gz)?$', full.names=TRUE))
if (!length(fastqs)) stop('No FASTQ files found')
files <- normalizePath(fastqs, mustWork=TRUE)
names(files) <- basename(files)
used <- character()
rows <- list()
for (f in files) {
  base <- basename(f)
  if (base %in% used) next
  if (grepl('_R2(_001)?\\.(fastq|fq)(\\.gz)?$', base)) {
    mate <- sub('_R2(_001)?(\\.(fastq|fq)(\\.gz)?)$', '_R1\\1\\2', base)
    if (mate %in% names(files)) next
    stop(sprintf('R2 without R1: %s', base))
  }
  r2name <- sub('_R1(_001)?(\\.(fastq|fq)(\\.gz)?)$', '_R2\\1\\2', base)
  is_r1 <- grepl('_R1(_001)?\\.(fastq|fq)(\\.gz)?$', base)
  r2 <- if (is_r1 && r2name %in% names(files)) files[[r2name]] else ''
  if (is_r1 && !nzchar(r2)) stop(sprintf('R1 without R2: %s (rename single-end file without _R1)', base))
  sample <- sub('(_R1(_001)?)?\\.(fastq|fq)(\\.gz)?$', '', base)
  rows[[length(rows)+1]] <- data.frame(sampleID=sample, reference=opt$reffile, read1=unname(f), read2=unname(r2))
  used <- c(used, base, if (nzchar(r2)) r2name)
}
df <- do.call(rbind, rows)
if (anyDuplicated(df$sampleID)) stop('Duplicate sample IDs detected')
write.csv(df, opt$outfile, row.names=FALSE, quote=TRUE, na='')
cat(sprintf('Wrote %d samples to %s\n', nrow(df), opt$outfile))
