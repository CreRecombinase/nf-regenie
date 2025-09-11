/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { REGENIE_STEP1; FIX_PRED; REGENIE_STEP2_SNV } from '../modules/local/regenie/processes.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REGENIE {

    main:
    // Parameter summary and versions
    ch_versions = Channel.empty()
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_regenie_software_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // Input checks
    if( !params.phenoFile || !params.covarFile ) {
      exit 1, 'ERROR: --phenoFile and --covarFile are required.'
    }
    def mode = params.input_mode
    if( !mode ) {
      if (params.bgen)      mode = 'bgen'
      else if (params.bed)  mode = 'bed'
      else if (params.pgen) mode = 'pgen'
    }
    if( !mode || !(mode in ['bgen','bed','pgen']) ) {
      exit 1, "ERROR: Provide one of --bgen/--bed/--pgen or set --input_mode to bgen|bed|pgen"
    }

    // Chromosomes
    Channel
      .from( (params.chr_list ?: '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X').split(',').collect{ it.trim() }.findAll{ it } )
      .set { CHR_CH }

    def bgen_has_at = (mode=='bgen' && params.bgen && params.bgen.contains('@'))
    def bed_has_at  = (mode=='bed'  && params.bed  && params.bed.contains('@'))
    def pgen_has_at = (mode=='pgen' && params.pgen && params.pgen.contains('@'))

    Channel.value([ mode,
                    mode=='bgen' ? params.bgen : (mode=='bed' ? params.bed : params.pgen),
                    params.bgi,
                    params.sample ]).set { GENO_STEP1 }

    def shard_geno = CHR_CH.map{ chr ->
      if( mode=='bgen' ) {
        def g = params.bgen
        def i = params.bgi
        if( bgen_has_at ) { g = g.replace('@', chr); i = i ? i.replace('@', chr) : null }
        return tuple(mode, g, i, params.sample, chr)
      } else if( mode=='bed' ) {
        def p = params.bed
        if( bed_has_at ) p = p.replace('@', chr)
        return tuple(mode, p, null, null, chr)
      } else {
        def p = params.pgen
        if( pgen_has_at ) p = p.replace('@', chr)
        return tuple(mode, p, null, null, chr)
      }
    }
    shard_geno.set { GENO_STEP2 }

    def KEEP_CH    = params.keep    ? Channel.value(file(params.keep))       : Channel.value(null)
    def REMOVE_CH  = params.remove  ? Channel.value(file(params.remove))     : Channel.value(null)
    def EXTRACT_CH = params.extract ? Channel.value(file(params.extract))    : Channel.value(null)
    def EXCLUDE_CH = params.exclude ? Channel.value(file(params.exclude))    : Channel.value(null)

    // Step 1
    REGENIE_STEP1(
      GENO_STEP1,
      file(params.phenoFile),
      file(params.covarFile),
      KEEP_CH, REMOVE_CH, EXTRACT_CH, EXCLUDE_CH,
      params.pheno_col_list ?: '', params.event_col_list ?: ''
    )

    // Fix pred list
    FIX_PRED(REGENIE_STEP1.out.pred_list)

    // Step 2 per chr
    REGENIE_STEP2_SNV(
      GENO_STEP2,
      file(params.phenoFile),
      file(params.covarFile),
      FIX_PRED.out.pred_list_fixed,
      KEEP_CH, REMOVE_CH, EXTRACT_CH, EXCLUDE_CH
    )

    emit:
    versions = ch_collated_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
