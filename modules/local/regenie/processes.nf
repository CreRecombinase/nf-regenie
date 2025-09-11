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
  path keepfile
  path removefile
  path extractfile
  path excludefile
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
  def genoFlag  = (geno_mode=='bgen') ? "--bgen ${geno_main} ${bgi ? "--bgi ${bgi}" : ''} ${sample ? "--sample ${sample}" : ''}"
                 : (geno_mode=='bed') ? "--bed ${geno_main}"
                 : "--pgen ${geno_main}"
  def phenoCols = pheno_cols ? "--phenoColList ${pheno_cols}" : ""
  def eventCols = (params.t2e && event_cols) ? "--eventColList ${event_cols}" : ""
  def strict    = params.strict ? "--strict" : ""
  def rint      = params.apply_rint ? "--apply-rint" : ""
  def lowmem    = params.lowmem ? "--lowmem --lowmem-prefix step1_tmp" : ""
  def keep      = keepfile ? "--keep ${keepfile}" : ""
  def remove    = removefile ? "--remove ${removefile}" : ""
  def extract   = extractfile ? "--extract ${extractfile}" : ""
  def exclude   = excludefile ? "--exclude ${excludefile}" : ""
  def threads   = "--threads ${task.cpus}"
  def bsize     = "--bsize ${params.bsize1}"

  """
  set -euo pipefail
  ./regenie \
    --step 1 \
    ${genoFlag} \
    --covarFile ${covar} \
    --phenoFile ${pheno} \
    ${traitFlag} ${phenoCols} ${eventCols} \
    ${keep} ${remove} ${extract} ${exclude} \
    ${bsize} ${threads} ${strict} ${rint} ${lowmem} \
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
  awk '{print}' ${pred_list} | while read -r pheno fname; do \
    base=\$(basename \"\${fname}\"); \
    cp -f \"\${fname}\" \"\${base}\"; \
    echo \"\${pheno} \${base}\"; \
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
  path keepfile
  path removefile
  path extractfile
  path excludefile

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
    --out test_step2_${chr}
  """
}

