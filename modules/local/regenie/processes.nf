process REGENIE_STEP1 {
  tag "step1"
  cpus params.step1_cpus ?: 8
  memory params.step1_mem ?: '32 GB'
  errorStrategy 'retry'
  maxRetries 2
  conda "${projectDir}/envs/regenie.yml"

  input:
  tuple val(geno_mode), val(geno_main), val(bgi), val(sample)
  path pheno
  path covar
  path keepfile, optional: true
  path removefile, optional: true
  path extractfile, optional: true
  path excludefile, optional: true
  val  pheno_cols
  val  event_cols

  output:
  path "*_pred.list", emit: pred_list
  path "*.loco*", optional: true, emit: loco
  path "*.prs*",  optional: true, emit: prs
  path "*.firth*", optional: true, emit: firth
  path "*.log", emit: logs

  script:
  def traitFlag = params.bt ? '--bt' : (params.t2e ? '--t2e' : '--qt')
  def refFirst  = params.ref_first ? "--ref-first" : ""
  def genoFlag  = (geno_mode=='bgen') ? "--bgen ${geno_main} ${bgi ? "--bgi ${bgi}" : ''} ${sample ? "--sample ${sample}" : ''} ${refFirst}"
                 : (geno_mode=='bed') ? "--bed ${geno_main} ${refFirst}"
                 : "--pgen ${geno_main}"
  def phenoCols = pheno_cols ? "--phenoColList ${pheno_cols}" : ""
  def phenoExcl = params.pheno_exclude_list ? "--phenoExcludeList ${params.pheno_exclude_list}" : ""
  def eventCols = (params.t2e && event_cols) ? "--eventColList ${event_cols}" : ""
  def covCols   = params.covar_col_list ? "--covarColList ${params.covar_col_list}" : ""
  def catCovars = params.cat_covar_list ? "--catCovarList ${params.cat_covar_list}" : ""
  def covExcl   = params.covar_exclude_list ? "--covarExcludeList ${params.covar_exclude_list}" : ""
  def strict    = params.strict ? "--strict" : ""
  def rint      = params.apply_rint ? "--apply-rint" : ""
  def cvopt     = params.cv ? "--cv ${params.cv}" : ""
  def loocv     = params.loocv ? "--loocv" : ""
  def nbopt     = params.nb ? "--nb ${params.nb}" : ""
  def minCase   = params.minCaseCount ? "--minCaseCount ${params.minCaseCount}" : ""
  def lowmem    = params.lowmem ? "--lowmem ${params.lowmem_prefix ? "--lowmem-prefix ${params.lowmem_prefix}" : "--lowmem-prefix step1_tmp"}" : ""
  def keep      = keepfile ? "--keep ${keepfile}" : ""
  def remove    = removefile ? "--remove ${removefile}" : ""
  def extract   = extractfile ? "--extract ${extractfile}" : ""
  def exclude   = excludefile ? "--exclude ${excludefile}" : ""
  def threads   = "--threads ${task.cpus}"
  def bsize     = "--bsize ${params.bsize1}"
  def force1    = params.force_step1 ? "--force-step1" : ""
  def printprs  = params.print_prs ? "--print-prs" : ""
  def wnull     = params.write_null_firth ? "--write-null-firth" : ""
  def compall   = params.compute_all ? "--compute-all" : ""
  def useRel    = params.use_relative_pred ? "--use-relative-path" : ""
  def tpheno    = params.tpheno_file ? "--tpheno-file ${params.tpheno_file}" : ""
  def tphIdx    = params.tpheno_indexCol ? "--tpheno-indexCol ${params.tpheno_indexCol}" : ""
  def tphIgn    = params.tpheno_ignoreCols ? "--tpheno-ignoreCols ${params.tpheno_ignoreCols}" : ""
  def iidOnly   = params.iid_only ? "--iid-only" : ""

  """
  set -euo pipefail
  ./regenie \
    --step 1 \
    ${genoFlag} \
    --covarFile ${covar} \
    --phenoFile ${pheno} \
    ${traitFlag} ${phenoCols} ${phenoExcl} ${eventCols} ${covCols} ${catCovars} ${covExcl} ${tpheno} ${tphIdx} ${tphIgn} ${iidOnly} \
    ${keep} ${remove} ${extract} ${exclude} \
    ${bsize} ${threads} ${strict} ${rint} ${cvopt} ${loocv} ${nbopt} ${minCase} ${lowmem} ${force1} ${printprs} ${wnull} ${compall} ${useRel} \
    --out fit_step1
  """
}

process FIX_PRED {
  tag "fix-pred-list"
  conda "${projectDir}/envs/regenie.yml"

  input:
  path pred_list

  output:
  path "fit_step1_pred.fixed.list", emit: pred_list_fixed
  path "*.loco*", optional: true
  path "*.prs*",  optional: true

  script:
  """
  set -euo pipefail
  awk '{print}' ${pred_list} | while read -r pheno fname; do
    base=\$(basename "\${fname}")
    cp -f "\${fname}" "\${base}"
    echo "\${pheno} \${base}"
  done > fit_step1_pred.fixed.list
  """
}

process REGENIE_STEP2_SNV {
  tag { "chr:${chr}" }
  cpus params.step2_cpus ?: 8
  memory params.step2_mem ?: '24 GB'
  errorStrategy 'retry'
  maxRetries 2
  conda "${projectDir}/envs/regenie.yml"

  input:
  tuple val(geno_mode), val(geno_shard), val(bgi), val(sample), val(chr)
  path pheno
  path covar
  path pred_list
  path keepfile, optional: true
  path removefile, optional: true
  path extractfile, optional: true
  path excludefile, optional: true

  output:
  path "*.regenie*", emit: sumstats
  path "*.ids", optional: true, emit: ids
  path "*.log", emit: logs

  script:
  def traitFlag  = params.bt ? '--bt' : (params.t2e ? '--t2e' : '--qt')
  def genoFlag   = (geno_mode=='bgen') ? "--bgen ${geno_shard} ${bgi ? "--bgi ${bgi}" : ''} ${sample ? "--sample ${sample}" : ''}"
                  : (geno_mode=='bed') ? "--bed ${geno_shard}" : "--pgen ${geno_shard}"
  def keep       = keepfile ? "--keep ${keepfile}" : ""
  def remove     = removefile ? "--remove ${removefile}" : ""
  def extract    = extractfile ? "--extract ${extractfile}" : ""
  def exclude    = excludefile ? "--exclude ${excludefile}" : ""
  def threads    = "--threads ${task.cpus}"
  def bsize      = "--bsize ${params.bsize2}"
  def firth      = params.firth ? "--firth ${params.approx ? '--approx' : ''}" : ""
  def spa        = params.spa ? "--spa" : ""
  def pthr       = "--pThresh ${params.pThresh}"
  def miss       = params.force_impute ? "--force-impute" : ""
  def gz         = params.gz ? "--gz" : ""
  def nosplit    = params.no_split ? "--no-split" : ""
  def minmac     = params.minMAC ? "--minMAC ${params.minMAC}" : ""
  def mininfo    = params.minINFO ? "--minINFO ${params.minINFO}" : ""
  def testmode   = params.test_mode ? "--test ${params.test_mode}" : ""
  def writesamp  = params.write_samples ? "--write-samples ${params.print_pheno ? '--print-pheno' : ''}" : ""
  def chrArg     = "--chr ${chr}"
  def sexspec    = params.sex_specific ? "--sex-specific ${params.sex_specific}" : ""
  def afcc       = params.af_cc ? "--af-cc" : ""
  def rangeOpt   = params.range ? "--range ${params.range}" : ""
  def startBlk   = params.starting_block ? "--starting-block ${params.starting_block}" : ""
  def htpOpt     = params.htp ? "--htp ${params.htp}" : ""
  def exactp     = params.exact_p ? "--exact-p" : ""
  def useNull    = params.use_null_firth ? "--use-null-firth ${params.use_null_firth}" : ""
  def firthse    = params.firth_se ? "--firth-se" : ""
  def forceqt    = params.force_qt ? "--force-qt" : ""
  def ignorepred = params.ignore_pred ? "--ignore-pred" : ""
  def niterOpt   = params.niter ? "--niter ${params.niter}" : ""
  def maxitNull  = params.maxiter_null ? "--maxiter-null ${params.maxiter_null}" : ""
  def maxstepN   = params.maxstep_null ? "--maxstep-null ${params.maxstep_null}" : ""
  def nautoOpt   = params.nauto ? "--nauto ${params.nauto}" : ""
  def maxCat     = params.maxCatLevels ? "--maxCatLevels ${params.maxCatLevels}" : ""
  def parReg     = params.par_region ? "--par-region ${params.par_region}" : ""
  def condList   = params.condition_list ? "--condition-list ${params.condition_list}" : ""
  def condFile   = params.condition_file ? "--condition-file ${params.condition_file}" : ""
  def interVar   = params.interaction ? "--interaction ${params.interaction}" : ""
  def interSnp   = params.interaction_snp ? "--interaction-snp ${params.interaction_snp}" : ""
  def interFile  = params.interaction_file ? "--interaction-file ${params.interaction_file}" : ""
  def interFileRF= params.interaction_file_reffirst ? "--interaction-file-reffirst" : ""
  def usePrs     = params.use_prs ? "--use-prs" : ""

  """
  set -euo pipefail
  ./regenie \
    --step 2 \
    ${genoFlag} \
    --covarFile ${covar} \
    --phenoFile ${pheno} \
    ${traitFlag} \
    --pred ${pred_list} \
    ${keep} ${remove} ${extract} ${exclude} \
    ${chrArg} ${bsize} ${threads} ${firth} ${spa} ${pthr} ${miss} \
    ${minmac} ${mininfo} ${testmode} ${writesamp} ${gz} ${nosplit} \
    ${sexspec} ${afcc} ${rangeOpt} ${startBlk} ${htpOpt} ${exactp} ${useNull} ${firthse} ${forceqt} ${ignorepred} \
    ${niterOpt} ${maxitNull} ${maxstepN} ${nautoOpt} ${maxCat} ${parReg} \
    ${condList} ${condFile} ${interVar} ${interSnp} ${interFile} ${interFileRF} ${usePrs} \
    --out test_step2_${chr}
  """
}

process REGENIE_STEP1_SPLIT_INIT {
  tag "step1-split-init"
  cpus params.step1_cpus ?: 8
  memory params.step1_mem ?: '32 GB'
  errorStrategy 'retry'
  maxRetries 2
  conda "${projectDir}/envs/regenie.yml"

  input:
  tuple val(geno_mode), val(geno_main), val(bgi), val(sample)
  path pheno
  path covar
  path keepfile, optional: true
  path removefile, optional: true
  path extractfile, optional: true
  path excludefile, optional: true
  val  pheno_cols
  val  event_cols

  output:
  path "*.log", emit: logs
  val { params.split_l0_prefix + '.master' }, emit: master

  script:
  def traitFlag = params.bt ? '--bt' : (params.t2e ? '--t2e' : '--qt')
  def refFirst  = params.ref_first ? "--ref-first" : ""
  def genoFlag  = (geno_mode=='bgen') ? "--bgen ${geno_main} ${bgi ? "--bgi ${bgi}" : ''} ${sample ? "--sample ${sample}" : ''} ${refFirst}"
                 : (geno_mode=='bed') ? "--bed ${geno_main} ${refFirst}"
                 : "--pgen ${geno_main}"
  def phenoCols = pheno_cols ? "--phenoColList ${pheno_cols}" : ""
  def phenoExcl = params.pheno_exclude_list ? "--phenoExcludeList ${params.pheno_exclude_list}" : ""
  def eventCols = (params.t2e && event_cols) ? "--eventColList ${event_cols}" : ""
  def covCols   = params.covar_col_list ? "--covarColList ${params.covar_col_list}" : ""
  def catCovars = params.cat_covar_list ? "--catCovarList ${params.cat_covar_list}" : ""
  def covExcl   = params.covar_exclude_list ? "--covarExcludeList ${params.covar_exclude_list}" : ""
  def strict    = params.strict ? "--strict" : ""
  def rint      = params.apply_rint ? "--apply-rint" : ""
  def cvopt     = params.cv ? "--cv ${params.cv}" : ""
  def loocv     = params.loocv ? "--loocv" : ""
  def nbopt     = params.nb ? "--nb ${params.nb}" : ""
  def minCase   = params.minCaseCount ? "--minCaseCount ${params.minCaseCount}" : ""
  def lowmem    = params.lowmem ? "--lowmem ${params.lowmem_prefix ? "--lowmem-prefix ${params.lowmem_prefix}" : "--lowmem-prefix step1_tmp"}" : ""
  def keep      = keepfile ? "--keep ${keepfile}" : ""
  def remove    = removefile ? "--remove ${removefile}" : ""
  def extract   = extractfile ? "--extract ${extractfile}" : ""
  def exclude   = excludefile ? "--exclude ${excludefile}" : ""
  def threads   = "--threads ${task.cpus}"
  def bsize     = "--bsize ${params.bsize1}"
  def force1    = params.force_step1 ? "--force-step1" : ""
  def printprs  = params.print_prs ? "--print-prs" : ""
  def wnull     = params.write_null_firth ? "--write-null-firth" : ""
  def compall   = params.compute_all ? "--compute-all" : ""
  def useRel    = params.use_relative_pred ? "--use-relative-path" : ""
  def tpheno    = params.tpheno_file ? "--tpheno-file ${params.tpheno_file}" : ""
  def tphIdx    = params.tpheno_indexCol ? "--tpheno-indexCol ${params.tpheno_indexCol}" : ""
  def tphIgn    = params.tpheno_ignoreCols ? "--tpheno-ignoreCols ${params.tpheno_ignoreCols}" : ""
  def iidOnly   = params.iid_only ? "--iid-only" : ""
  def split     = (params.split_l0_prefix && params.split_l0_n) ? "--split-l0 ${params.split_l0_prefix},${params.split_l0_n}" : ""

  """
  set -euo pipefail
  ./regenie \
    --step 1 \
    ${genoFlag} \
    --covarFile ${covar} \
    --phenoFile ${pheno} \
    ${traitFlag} ${phenoCols} ${phenoExcl} ${eventCols} ${covCols} ${catCovars} ${covExcl} ${tpheno} ${tphIdx} ${tphIgn} ${iidOnly} \
    ${keep} ${remove} ${extract} ${exclude} \
    ${bsize} ${threads} ${strict} ${rint} ${cvopt} ${loocv} ${nbopt} ${minCase} ${lowmem} ${force1} ${printprs} ${wnull} ${compall} ${useRel} \
    ${split} \
    --out fit_step1_split
  """
}

process REGENIE_STEP1_RUN_L0 {
  tag { "run-l0:${k}" }
  cpus params.step1_cpus ?: 4
  memory params.step1_mem ?: '16 GB'
  errorStrategy 'retry'
  maxRetries 2
  conda "${projectDir}/envs/regenie.yml"

  input:
  path master
  val  k

  output:
  path "*.log", emit: logs

  script:
  """
  set -euo pipefail
  ./regenie --run-l0 ${master} ${k} --out runl0_${k}
  """
}

process REGENIE_STEP1_RUN_L1 {
  tag "run-l1"
  cpus params.step1_cpus ?: 8
  memory params.step1_mem ?: '32 GB'
  errorStrategy 'retry'
  maxRetries 2
  conda "${projectDir}/envs/regenie.yml"

  input:
  path master

  output:
  path "*_pred.list", emit: pred_list
  path "*.loco*", optional: true, emit: loco
  path "*.prs*",  optional: true, emit: prs
  path "*.firth*", optional: true, emit: firth
  path "*.log", emit: logs

  script:
  def keepL0   = params.keep_l0 ? "--keep-l0" : ""
  def l1ph     = params.l1_pheno_list ? "--l1-phenoList ${params.l1_pheno_list}" : ""
  """
  set -euo pipefail
  ./regenie --run-l1 ${master} ${l1ph} ${keepL0} --out fit_step1
  """
}
