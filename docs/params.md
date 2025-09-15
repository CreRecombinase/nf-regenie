# nf-regenie parameters

Parameters for the nf-core REGENIE pipeline. Defines grouped inputs and options for Step 1/2, phenotypes/covariates, genotype sources, and I/O.

## Input/output

Define where the pipeline should save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `outdir` | Output directory for pipeline results. Pipeline writes logs, reports and REGENIE outputs under this path. | `string` | results |  |  |

## Genotype input

Specify genotype inputs for Step 1/2. Use '@' in paths to expand per-chromosome shards.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `bgen` | Input BGEN genotype file path or pattern. Use '@' placeholder to expand per-chromosome sharding (e.g. s3://.../chr@.bgen). | `['string', 'null']` |  |  |  |
| `bgi` | Index .bgi file corresponding to BGEN input. If 'bgen' uses '@', this may also use '@' to expand per-chromosome. | `['string', 'null']` |  |  |  |
| `sample` | Sample file corresponding to BGEN input ('.sample'). | `['string', 'null']` |  |  |  |
| `bed` | PLINK bed/bim/fam prefix. Provide path prefix without extension; pipeline expects .bed/.bim/.fam. | `['string', 'null']` |  |  |  |
| `pgen` | PLINK2 pgen/pvar/psam prefix. Provide path prefix without extension; pipeline expects .pgen/.pvar/.psam. | `['string', 'null']` |  |  |  |
| `input_mode` | Force genotype source for Step 1/2 ('bgen'|'bed'|'pgen'). If null, Step 1 prefers BED/PGEN then BGEN; Step 2 prefers BGEN then BED/PGEN. | `['string', 'null']` | None |  |  |
| `ref_first` | Use first allele as reference for BGEN or PLINK BED input (--ref-first). | `boolean` |  |  |  |
| `keep` | Inclusion file listing FID IID pairs to retain. | `['string', 'null']` | None |  |  |
| `remove` | Exclusion file listing FID IID pairs to remove. | `['string', 'null']` | None |  |  |
| `extract` | Variant inclusion list (IDs) to keep. | `['string', 'null']` | None |  |  |
| `exclude` | Variant exclusion list (IDs) to remove. | `['string', 'null']` | None |  |  |

## Phenotypes/covariates

Phenotype and covariate inputs for Step 1 and 2.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `phenoFile` | Phenotype file with header 'FID IID Y1 ...'. Missing values must be 'NA'. Used in Step 1 and Step 2. | `string` |  |  |  |
| `covarFile` | Covariate file with header 'FID IID C1 ...'. Samples with missing covariates are removed. Use same file in Step 1 and 2. | `string` |  |  |  |
| `pheno_col_list` | Comma-separated list of phenotype names to include (equivalent to --phenoColList). | `string` |  |  |  |
| `event_col_list` | Comma-separated list of event time columns for time-to-event analyses (equivalent to --eventColList). | `string` |  |  |  |
| `pheno_exclude_list` | Comma-separated list of phenotypes to exclude (--phenoExcludeList). | `['string', 'null']` | None |  |  |
| `covar_col_list` | Comma-separated list of covariates to include (--covarColList). | `['string', 'null']` | None |  |  |
| `cat_covar_list` | Comma-separated list of categorical covariates (--catCovarList). | `['string', 'null']` | None |  |  |
| `covar_exclude_list` | Comma-separated list of covariates to exclude (--covarExcludeList). | `['string', 'null']` | None |  |  |
| `tpheno_file` | Transposed phenotype file to use instead of standard format (--tpheno-file). | `['string', 'null']` | None |  |  |
| `tpheno_indexCol` | Index of phenotype name column in transposed phenotype file (--tpheno-indexCol). | `['integer', 'null']` | None |  |  |
| `tpheno_ignoreCols` | Comma-separated indices of columns to ignore in transposed phenotype file (--tpheno-ignoreCols). | `['string', 'null']` | None |  |  |
| `iid_only` | Transposed phenotype header only contains IID (assumes FID=IID) (--iid-only). | `boolean` |  |  |  |

## Step 1 options

Whole-genome model fitting and prediction settings.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `bt` | Analyze binary traits in Step 1/2 (--bt). | `boolean` |  |  |  |
| `t2e` | Analyze time-to-event traits (--t2e). | `boolean` |  |  |  |
| `bsize1` | Block size for Step 1 (--bsize). | `integer` | 100 |  |  |
| `lowmem` | Enable low-memory mode in Step 1 to use disk-backed temporary files (--lowmem). | `boolean` | True |  |  |
| `apply_rint` | Apply rank-based inverse normal transform to QTs in Step 1 (--apply-rint). | `boolean` |  |  |  |
| `strict` | Remove samples with missing values at any phenotype in Step 1 and Step 2 (--strict). | `boolean` |  |  |  |
| `cv` | K-fold cross validation for Step 1 when applicable (--cv). | `['integer', 'null']` | None |  |  |
| `loocv` | Use leave-one-out cross validation in Step 1 (--loocv). | `boolean` |  |  |  |
| `nb` | Number of blocks for cross validation when applicable (--nb). | `['integer', 'null']` | None |  |  |
| `minCaseCount` | Minimum case count required for binary traits in Step 1 (--minCaseCount). | `['integer', 'null']` | None |  |  |
| `lowmem_prefix` | Prefix for low-memory temporary files in Step 1 (--lowmem-prefix). | `['string', 'null']` | None |  |  |
| `force_step1` | Force Step 1 to run even if predictions already exist (--force-step1). | `boolean` |  |  |  |
| `print_prs` | Write whole-genome PRS (non-LOCO) predictions in Step 1 (--print-prs). | `boolean` |  |  |  |
| `write_null_firth` | Write null estimates for approximate Firth in Step 1/2 (--write-null-firth). | `boolean` |  |  |  |
| `compute_all` | Write null Firth estimates for all chromosomes regardless of genotype file (--compute-all). | `boolean` |  |  |  |
| `use_relative_pred` | Use relative paths in Step 1 prediction list to ease staging in Nextflow. | `boolean` | True |  |  |
| `split_l0_prefix` | Prefix for split-L0 mode. Generates '<prefix>.master' and L0/L1 outputs. | `['string', 'null']` | None |  |  |
| `split_l0_n` | Number of split-L0 jobs (K) to run in parallel. | `['integer', 'null']` | None |  |  |
| `l1_pheno_list` | Comma-separated subset of phenotypes to analyze at L1 when using split-L0. | `['string', 'null']` | None |  |  |
| `keep_l0` | Keep L0 predictions on disk after L1 completes (--keep-l0). | `boolean` |  |  |  |

## Step 2 options

Association testing and summary statistics settings.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `chr_list` | Comma-separated list of chromosomes to test in Step 2 (passed to --chr per value). | `string` | 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X |  |  |
| `bsize2` | Block size for Step 2 (--bsize). | `integer` | 200 |  |  |
| `firth` | Use Firth LRT fallback for small p-values in Step 2 (--firth). | `boolean` | True |  |  |
| `approx` | Use approximate Firth LRT when --firth is enabled (--approx). | `boolean` | True |  |  |
| `spa` | Use Saddlepoint approximation fallback for small p-values (--spa). | `boolean` |  |  |  |
| `pThresh` | P-value threshold below which to apply Firth/SPA corrections (--pThresh). | `number` | 0.05 |  |  |
| `minMAC` | Minimum minor allele count threshold when testing variants (--minMAC). | `integer` | 5 |  |  |
| `minINFO` | Minimum imputation info score for variant inclusion in Step 2 (--minINFO). | `['number', 'null']` | None |  |  |
| `test_mode` | Genetic model to test in Step 2 (--test). If null, additive is used. | `['string', 'null']` | None |  |  |
| `gz` | Gzip-compress REGENIE outputs (--gz). | `boolean` | True |  |  |
| `no_split` | Write summary statistics for all traits to a single file (--no-split). | `boolean` |  |  |  |
| `force_impute` | Keep and impute missing QT observations in Step 2 (--force-impute). | `boolean` |  |  |  |
| `sex_specific` | Perform sex-specific analyses for 'male' or 'female' (--sex-specific). | `['string', 'null']` | None |  |  |
| `af_cc` | Output A1FREQ in cases and controls separately in Step 2 results (--af-cc). | `boolean` |  |  |  |
| `range` | Chromosome region to test in Step 2 (CHR:MINPOS-MAXPOS) (--range). | `['string', 'null']` | None |  |  |
| `starting_block` | Start Step 2 at a specific block/set number (--starting-block). | `['integer', 'null']` | None |  |  |
| `htp` | Write Step 2 summary stats in HTP format (--htp). | `['string', 'null']` | None |  |  |
| `exact_p` | Write exact p-values for logistic models where applicable (--exact-p). | `boolean` |  |  |  |
| `use_null_firth` | List file of null estimates for approximate Firth to reuse in Step 2 (--use-null-firth). | `['string', 'null']` | None |  |  |
| `firth_se` | Compute SE from effect size and LRT p-value under Firth correction (--firth-se). | `boolean` |  |  |  |
| `force_qt` | Force QT model for binary traits (--force-qt). | `boolean` |  |  |  |
| `ignore_pred` | Ignore Step 1 predictions in Step 2 (--ignore-pred). | `boolean` |  |  |  |
| `niter` | Maximum number of iterations for logistic regression (--niter). | `['integer', 'null']` | None |  |  |
| `maxiter_null` | Max iterations for logistic model with Firth penalty under the null (--maxiter-null). | `['integer', 'null']` | None |  |  |
| `maxstep_null` | Max step size for logistic model with Firth penalty under the null (--maxstep-null). | `['integer', 'null']` | None |  |  |
| `nauto` | Number of autosomal chromosomes for non-human studies (--nauto). | `['integer', 'null']` | None |  |  |
| `maxCatLevels` | Maximum number of levels allowed for categorical covariates (--maxCatLevels). | `['integer', 'null']` | None |  |  |
| `par_region` | Build code or bounds for chrX PAR regions (e.g., b36/b37/b38 or 'start,end' of non-PAR) (--par-region). | `['string', 'null']` | None |  |  |
| `condition_list` | Comma-separated list of variant IDs to condition on (--condition-list). | `['string', 'null']` | None |  |  |
| `condition_file` | File listing variant IDs to condition on (--condition-file). | `['string', 'null']` | None |  |  |
| `interaction` | Covariate name for genotype-covariate interaction testing (--interaction). | `['string', 'null']` | None |  |  |
| `interaction_snp` | Variant ID to test genotype-by-variant interaction (--interaction-snp). | `['string', 'null']` | None |  |  |
| `interaction_file` | File with variant-covariate interaction pairs (--interaction-file). | `['string', 'null']` | None |  |  |
| `interaction_file_reffirst` | Interpret first allele as reference in interaction file (--interaction-file-reffirst). | `boolean` |  |  |  |
| `use_prs` | Use PRS predictions from Step 1 instead of LOCO (--use-prs). | `boolean` |  |  |  |

## Resources

Optional resource overrides for Nextflow processes.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `step1_cpus` | Override CPU threads for Step 1 processes in Nextflow (maps to --threads). | `['integer', 'null']` | None |  |  |
| `step1_mem` | Override memory for Step 1 processes in Nextflow (e.g. '32 GB'). | `['string', 'null']` | None |  |  |
| `step2_cpus` | Override CPU threads for Step 2 processes in Nextflow (maps to --threads). | `['integer', 'null']` | None |  |  |
| `step2_mem` | Override memory for Step 2 processes in Nextflow (e.g. '24 GB'). | `['string', 'null']` | None |  |  |
