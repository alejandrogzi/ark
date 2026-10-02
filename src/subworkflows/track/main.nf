/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RSYNC_SSH } from '../../modules/custom/ssh/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOAD_TRACK {
    take:
      bigbed                 // channel: [ [ meta ], [ bigbed ] ]
      user                   // string
      server                 // string
      target_dir             // string
      web                    // string
      species                // string
      ch_versions            // [ meta, versions.yml ]

    main:
      // INFO: one call per bigBed class (pass, trash, fusions, nmd, ...), only with params.load_track
      // INFO: per bigBed: rsync to <user>@<server>:<target_dir>/<species>/isopipe/<name>.bb,
      // INFO: then symlink it into <web>/<species>/ on the server so the browser can serve it
      // WARN: needs passwordless ssh to <server>; nothing is emitted except versions
      RSYNC_SSH(
        bigbed,
        user,
        server,
        target_dir,
        web,
        species,
      )

      ch_versions = ch_versions.mix(RSYNC_SSH.out.versions)

    emit:
      versions              = ch_versions
}
