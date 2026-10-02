/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TWOBIT_TO_FA } from '../../modules/custom/ucsc/twobittofa/main'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/custom/gunzip/main'
include { CHROMSIZE } from '../../modules/custom/chromsize/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 
workflow GENOME {
    take:
      genome  // file: /path/to/genome.{2bit/fasta}

    main:
      // INFO: accepts .2bit, .gz (gzipped FASTA) or plain FASTA
      // INFO: genome and chrom_sizes are value channels (bare path, no meta), reusable by any number of tasks
      ch_versions = Channel.empty()
      ch_fasta = Channel.empty()

      def genome_file = file(genome, checkIfExists: true)
      def genome_path = genome_file.toString()

      // INFO: chromsize reads the input file as given (before any conversion)
      ch_chrom_sizes = CHROMSIZE([[:], genome_file]).chromsize.map { it[1] }

      // INFO: if fasta is .2bit or .gz, convert or uncompress it
      if (genome_path.endsWith(".2bit")) {
          ch_fasta = TWOBIT_TO_FA([[:], genome_file]).fasta.map { it[1] }
          ch_versions = ch_versions.mix(TWOBIT_TO_FA.out.versions)
      } else if (genome_path.endsWith(".gz")) {
          ch_fasta = GUNZIP_FASTA([[:], genome_file]).gunzip.map { it[1] }
          ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
      } else {
          ch_fasta = Channel.value(genome_file)
      }

      ch_versions = ch_versions.mix(CHROMSIZE.out.versions)

    emit:
      genome      = ch_fasta // file: /path/to/genome.fasta
      chrom_sizes = ch_chrom_sizes
      versions    = ch_versions
}

