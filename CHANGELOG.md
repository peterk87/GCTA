# Changelog

All notable changes to [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
(historical GCTA tags often used `beta` suffixes).

Release notes through **1.95.3** are adapted from the official
[Download / Update log](https://yanglab.westlake.edu.cn/software/gcta/#Download).
Where possible, entries link to matching commits in
[JianYang-Lab/GCTA](https://github.com/JianYang-Lab/GCTA) (via `GCTA_VERSION` bumps in `src/config.h` or commit messages).
Git tags were not present upstream at the time this file was created; compare URLs use commit SHAs until tags are published.

---

## [Unreleased]

Work on branch `dev` (fork [`peterk87/GCTA`](https://github.com/peterk87/GCTA)), based on upstream `v1.95.3`.

### Added

- COJO characterization / golden-test harness (`scripts/ci/golden_test.sh`, `compare_cojo.py`, `numdiff.py`) with **exact** and **`--tol`** modes.
- Deterministic COJO fixtures under `tests/fixtures/` and blessed goldens under `tests/golden/` (`scripts/ci/make_fixtures.py`).
- GoogleTest unit tests for `CommFunc` / `StrFunc` / `StatFunc` and mirrored COJO inverse-update oracles (`tests/unit/`).
- GitHub Actions CI (lint, unit tests, MKL build, golden job) under `.github/workflows/`.
- Mostly-static `gcta64` via `-DGCTA_STATIC_EXE=ON` (static MKL/GSL/libgcc/libstdc++); CI uploads `gcta64-linux-x86_64-static` and comments the download URL on PRs.
- This `CHANGELOG.md`.

### Changed

- **COJO / PLINK I/O:** `--extract-region-bp` is applied during BIM read (`set_bim_region_filter`) so only in-region SNPs are indexed; BED rows are sought via `_snp_bed_row` (avoids loading a full-chromosome BIM per LD-block process).
- **COJO deferred genotype load:** for single-bfile `--cojo-slct` / `--cojo-joint` / `--cojo-cond` / `--cojo-sblup` (without `--update-freq` / dosage Rsq filters), BIM records BED row indices and genotypes are decoded only for `.ma`∩panel SNPs inside `init_massoc`, then `--maf` / `--max-maf` are applied — avoids materializing the full-chromosome BED for full-chr COJO.
- **COJO MA match:** with the BIM region filter active, `.ma` lines outside the region index are skipped early (phenotypic-variance median still uses all lines for parity).
- **COJO stepwise:** `insert_B_and_Z` / `erase_B_and_Z` maintain dense `_B_i` / `_B_N_i` with rank-1 Schur updates (append-then-permute for sorted SNP order) instead of rebuilding via `SimplicialLDLT` each step; incremental insert keeps stock’s three guards (Schur/`denom>0`, LDLT `cond(D)>30`, per-SNP collinearity).
- **COJO LD fill:** `_bp_order` + window bounds restrict `init_Z` / insert-Z genotype dots to SNPs within `--cojo-wind`.
- **COJO LD cache:** dense `_Z_cache` / `_Z_N_cache` with OpenMP window row fills; insert/erase append or drop one row instead of rebuilding sparse `_Z`.
- **COJO genotype dots:** optional centered genotype cache (`_cojo_X`, capped ~320MB per process) plus bit-packed keep-subset planes and popcnt kernel in `cojo_snp_dot`; cache-off makex path still hoists the outer vector.
- **PLINK BED I/O:** `read_bedfile` uses 64 MiB block-buffered reads and OpenMP SNP decode, including when seeking via `_snp_bed_row`.
- **OpenMP:** `mainV1` / `gcta64` / GCTA libs link `OpenMP::OpenMP_CXX` on all platforms (previously Apple-only); `option()` always sets `OMP_NUM_THREADS` and `omp_set_num_threads` from `--thread-num` / `--threads`.

### Fixed

- Golden harness no longer aborts the whole run on the first `DIFFERENT` comparator under `set -euo pipefail` (exit codes captured via `if txt=$(…)`).
- `erase_B_and_Z` always refreshes inverses (even when `_Z` was never built) — intentional semantics fix vs stock’s early return when `_Z_N` was empty (e.g. backward-only `slct_stay`).
- Incremental COJO insert restores non-PD (`denom≤0`) and `cond>30` guards that were dropped when switching off per-step LDLT solve.
- `cojo_snp_dot` cache-disabled fallback no longer rematerializes both genotype vectors every pair.
- Legacy `option()` path never called `omp_set_num_threads`, so `#pragma omp` in COJO Z fills was a no-op even when `--thread-num` was set.
- CI / CTest primary golden gate uses `--tol 1e-6`; exact byte mode kept as `cojo_golden_exact`.
- CI lint no longer py_compiles gitignored `sandbox/` scripts; golden tests run against the mostly-static `gcta64` artifact (avoids MKL/GSL shared-library version skew between jobs).

---

## [1.95.3] - 2026-07-10 — [`0dc78f0`](https://github.com/JianYang-Lab/GCTA/commit/0dc78f0)

### Added

- Integrate MCP (Model Context Protocol) service to enable AI tools such as Claude Codex to invoke and interact with the program.

### Commits

- [`0dc78f0`](https://github.com/JianYang-Lab/GCTA/commit/0dc78f0) — bump GCTA_VERSION to v1.95.3
- [`8798554`](https://github.com/JianYang-Lab/GCTA/commit/8798554) — add AI-MCP


## [1.95.2] - 2026-06-09 — [`4ee76be`](https://github.com/JianYang-Lab/GCTA/commit/4ee76be)

### Changed

- --heidi-thresh 0 now disables HEIDI (equivalent to 0 0) with a notice, and the default value 0.01 is applied and reported when --heidi-thresh is omitted.
- HEIDI is also automatically disabled with a warning when a trait has 5 or fewer clumped IVs.

### Commits

- [`4ee76be`](https://github.com/JianYang-Lab/GCTA/commit/4ee76be) — optimize --heidi-thresh


## [1.95.1] - 2026-02-09 — [`69d42e2`](https://github.com/JianYang-Lab/GCTA/commit/69d42e2)

### Fixed

- Fixed a bug causing fastGWA to falsely report a missing --pheno argument when a phenotype file is correctly specified.

### Commits

- [`69d42e2`](https://github.com/JianYang-Lab/GCTA/commit/69d42e2) — fix: keep --pheno when used with fastGWA
- [`d8e91d8`](https://github.com/JianYang-Lab/GCTA/commit/d8e91d8) — chore: update README
- [`9095ce1`](https://github.com/JianYang-Lab/GCTA/commit/9095ce1) — chore: update macos package script


## [1.95.0] - 2025-07-22 — [`3dbff08`](https://github.com/JianYang-Lab/GCTA/commit/3dbff08)

### Changed

- The structure of the executable files has been updated.

### Fixed

- Fixed a bug related to the --bgen error message.

### Commits

- [`af09f7c`](https://github.com/JianYang-Lab/GCTA/commit/af09f7c) — fix bgen error messages
- [`3dbff08`](https://github.com/JianYang-Lab/GCTA/commit/3dbff08) — updata version → v1.95.0
- [`299cca7`](https://github.com/JianYang-Lab/GCTA/commit/299cca7) — Fix missing so
- [`a66e042`](https://github.com/JianYang-Lab/GCTA/commit/a66e042) — Fix missing so
- [`7738ed2`](https://github.com/JianYang-Lab/GCTA/commit/7738ed2) — Add AppImage packaging
- [`4708b1e`](https://github.com/JianYang-Lab/GCTA/commit/4708b1e) — Support Mac
- [`405cf94`](https://github.com/JianYang-Lab/GCTA/commit/405cf94) — Add windows support


## [1.94.4] - 2025-03-13 — [`3ecd226`](https://github.com/JianYang-Lab/GCTA/commit/3ecd226)

*Website release 13 Mar 2025; `GCTA_VERSION` set to v1.94.4 in fork commit `3ecd226` (2025-05-30). GRM-parts fix is `9ab283a`.*

### Fixed

- Fixed a bug involving large numbers of GRM parts.

### Commits

- [`9ab283a`](https://github.com/JianYang-Lab/GCTA/commit/9ab283a) — fixed GRM part number issue
- [`3ecd226`](https://github.com/JianYang-Lab/GCTA/commit/3ecd226) — Fork from original GCTA (config → v1.94.4)


## [1.94.3] - 2024-12-17

### Added

- Added the --reml-est-fix-varcov option. This option displays the variance-covariance matrix of the estimated fixed effects, which is used to calculate standard errors in mixed-model analyses.

### Fixed

- Fixed several minor bugs.


## [1.94.2] - 2023-10-13

### Changed

- Before using the mac version of the software, you need to configure environment variables. The configuration manual is located in the downloaded software directory.

### Notes

- The Mac version of the software is updated to be compatible with the ARM architeture.


## [1.94.1] - 2022-08-01 — [`0bd4c0f`](https://github.com/JianYang-Lab/GCTA/commit/0bd4c0f)

### Added

- Added mBAT module.
- Added a QC step in fastBAT.

### Fixed

- Fixed bugs.

### Notes

- Download: Linux gcta-1.94.1-linux-kernel-3-x86_64.zip, MacOS gcta-1.94.1-macOS-x86_64.zip, Windows gcta-1.94.1-Win-x86_64.zip.

### Commits

- [`a92785e`](https://github.com/JianYang-Lab/GCTA/commit/a92785e) — add mBAT module
- [`5536e12`](https://github.com/JianYang-Lab/GCTA/commit/5536e12) — add fastBAT with QC
- [`0bd4c0f`](https://github.com/JianYang-Lab/GCTA/commit/0bd4c0f) — introduce src/config.h as v1.94.1


## [1.94.0beta] - 2022-01-04

### Changed

- Changed the buffer size of line field from 50 to 512 bytes in ACAT.

### Notes

- A major update of the software to be compatible with ARM architecture (credits to the openEuler BIO-SIG).
- Proofreading of all the error and warning messages.
- Download: Linux gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip, MacOS gcta_v1.94.0Beta_macOS.zip, Windows gcta_v1.94.0Beta_windows_x86_64.zip.


## [1.93.3beta2] - 2021-08-17

### Fixed

- Fixed a bug when running fastGWA with genotype files in BGEN format.
- Fixed an issue when running fastGWA with a stringent filtering which leads to not sufficient null SNPs (i.e., < 100) to estimate the parameter gamma.

### Notes

- Download: Linux gcta_1.93.3beta2.zip.


## [1.93.3beta] - 2021-06-01 — [`0fbe2d3`](https://github.com/JianYang-Lab/GCTA/commit/0fbe2d3)

### Added

- Added a new module fastGWA-GLMM (a resource-efficient generalized linear mixed model association tool for biobank-scale data).
- Added a new module fastGWA-BB (a set-based burden test for binary traits based on the framework of fastGWA-GLMM).
- Added a new module ACAT-V (a very efficient summary-level set-based test that only requires GWAS summary statistics, originally proposed by Liu et al, 2019).

### Commits

- [`0fbe2d3`](https://github.com/JianYang-Lab/GCTA/commit/0fbe2d3) — v1.93.3


## [1.93.2beta] - 2020-05-08 — [`9571d1a`](https://github.com/JianYang-Lab/GCTA/commit/9571d1a)

### Fixed

- Fixed a bug in --bgen when there are missing missing genotypes.
- Fixed a bug in --update-sex.

### Notes

- Download: Linux gcta_1.93.2beta.zip, macOS gcta_1.93.2beta_mac.zip, Windows gcta_1.93.2beta_win.zip.

### Commits

- [`9571d1a`](https://github.com/JianYang-Lab/GCTA/commit/9571d1a) — v1.93.2


## [1.93.1beta] - 2020-04-01

### Added

- Added two new flags in fastGWA: --model-only to save the estimated fastGWA model parameters and --load-model to load the saved estimates for association tests.

### Changed

- Changed the allele frequency calculation method for ChrX to coordinate with the corresponding changes in PLINK2.
- Updated fastGWA to use "-9" as a missing value symbol in phenotype or covariate data.

### Fixed

- Fixed a bug in fastGWA when dealing with ChrX, and a bug when the number of SNPs to calibrate the gamma parameter is too small.
- Fixed a bug in --update-freq.
- Fixed a bug in COJO when performed in conjunction with the --diff-freq flag.
- Fixed a bug in a few analyses when no gender information is present in the .fam file.
- Fixed a bug when loading sample information using the --sample flag.
- Fixed a bug in reporting the "Illegal instruction" error for old versions of CPU.
- Fixed a bug in computing the likelihood value in the within-family REML analysis.
- Fixed a bug in checking case/control data in REML.


## [1.93.0beta] - 2019-12-09

### Added

- Added a flag --bgen to input genotype data (including imputed dosage data) in bgen format (>=v1.2).
- Added a flag --mbgen to input genotype data in multiple bgen files.
- Added a flag --pfile to input genotype data in pgen format.
- Added a flag --mpfile to input genotype data in multiple pgen files.
- Added a flag --bpfile to input genotype data in hybrid pgen format (i.e., *.pgen, *.bim and *.fam).
- Added a flag --mbpfile to input genotype data in multiple sets of hybrid pgen files.
- Added a flag --geno to filter out individuals based on genotype missingness rate.
- Added a flag --info to filter out SNPs based on imputation INFO score.
- Added a flag --recodet to output a transposed matrix of the genotypes.
- Allows duplicated SNP IDs in genotype data.
- Added a flag --save-fastGWA-mlm-residual to output fastGWA residuals.

### Fixed

- Fixed a bug in fastGWA for chromosome X.
- Fixed a bug in --make-grm-xchr.

### Notes

- Amended the flag --dc to support fastGWA analysis of SNPs on chromosome X based on different dosage compensation models.


## [1.92.4beta2] - 2019-11-25

### Changed

- Changed the flag --fastGWA-lmm to --fastGWA-mlm.

### Fixed

- Fixed a bug in GSMR (HEIDI outlier test).

### Notes

- Download: Linux gcta_1.92.4beta2.zip, Windows gcta_1.92.4beta2_win.zip, Mac gcta_1.92.4beta2_mac.zip


## [1.92.4beta] - 2019-09-23

### Fixed

- Fixed a bug in --make-grm when the MAFs of some SNPs are 0.
- Fixed a bug in --make-grm-d and --make-grm-alg.


## [1.92.3beta3] - 2019-08-20

### Fixed

- Fixed a bug related to the convergence of fastGWA-REML in some rare scenarios.


## [1.92.3beta2] - 2019-08-12

### Fixed

- Fixed a bug in --maf. The bug only occurred in rare scenarios with specific numbers of variants.


## [1.92.3beta] - 2019-08-09

### Changed

- Updated the fastGWA module with fastGWA-REML as the default method for variance component estimation and the GRAMMAR-GAMMA approximation as the default method to compute test-statistics. The updated version is ~10 times faster than the previous version.
- Improved the speed of --make-grm and --make-grm-part by ~3-fold.


## [1.92.2beta] - 2019-06-18

### Changed

- Modified mtCOJO to accept LD score files with 4 columns.
- Changed flag --gsmr-beta to --gsmr2-beta.
- Improved the performance of fastGWA.

### Fixed

- Fixed a bug in GWAS simulation.
- Fixed a bug in fastGWA if a p-value is extremely small.


## [1.92.1beta6] - 2019-04-13

### Fixed

- Fixed a bug in --reml-bivar due to the update of Linux compiler.


## [1.92.1beta5] - 2019-04-01

### Added

- Added a flag '--reml-res-diag' to specify the diagonal elements of the residual correlation matrix in REML.
- Added a new module fastGWA (an extremely resource-efficient tool for mixed linear model association analysis of biobank-scale GWAS data).


## [1.92.0beta3] - 2019-02-01

### Fixed

- Fixed a bug in COJO for some circumstances where the standard errors of SNP effects are extremely small.


## [1.92.0beta] - 2019-01-22

### Added

- Added a flag '--gsmr-beta' to use a testing version of the HEIDI-outlier method.

### Changed

- Modified '--make-grm-alg' so that it can be used in combination with –make-grm-part.


## [1.91.7beta] - 2018-10-08

### Added

- Added a flag (--mtcojo-bxy) in the mtCOJO analysis to read the effects of covariates on trait from a user-specified file.
- Added a multi-SNP-based HEIDI-outlier test in the HEIDI-outlier analysis.
- Added a function in mtCOJO to compute the effects of covariates on trait from a genetic correlation analysis if there are not enough SNPs to perform the GSMR analysis.

### Changed

- Modified mtCOJO and GSMR to read summary data from compressed text files.
- Modified HEIDI-outlier to save the removed pleiotropic SNPs in text file.
- Modified GSMR to sort SNPs by chi-squared values in the clumping analysis.
- Changed the flag --mlma-no-adj-covar to --mlma-no-preadj-covar.

### Fixed

- Fixed a bug in the GSMR effect plot.
- Fixed a bug in --reml-bivar.


## [1.91.6beta] - 2018-08-17

### Added

- Added a function to check the consistency of allele frequency between the GWAS summary data and the reference sample in COJO (the --diff-freq flag).
- Added a new flag --unify-grm to unify the order of the IDs in multiple GRM files.

### Changed

- Changed the criterion of selecting the top associated SNP by p-value in COJO to that by chi-squared value to avoid the problem of having extremely small p-values (e.g. those = 1e-300).

### Fixed

- Fixed a bug in COJO in the Windows version.
- Fixed an issue related to allele frequency in COJO when the first allele differs between the GWAS summary data and the LD reference sample.
- Fixed a bug when manipulating the GRM in PCA.


## [1.91.5beta] - 2018-07-07

### Added

- Added a flag --diff-freq to check difference in allele frequency between data sets in the GSMR and mtCOJO analyses.
- Added a flag --mbfile to merge multiple BED files (e.g. genotype data of each chromsome saved in a separate BED file) into a single BED file.

### Changed

- Removed flags --clump-p1 and --heidi-snp from the GSMR and mtCOJO analyses
- Improved compatibility with old Linux version.

### Fixed

- Fixed a bug in GSMR when there is a very small number of SNPs used to run an HEIDI-outlier analysis.
- Fixed a memory issue with the flag --make-grm-x.
- Fixed a build stack issue in the Windows version.
- Fixed a rare thread freezing with --make-grm.

### Notes

- The flag --gsmr-snp has been superseded by --gsmr-snp-min.


## [1.91.4beta] - 2018-04-17

### Added

- Added --mbfile in GRM functions to proceed genotypes stored in multiple PLINK files.
- Added an additional option --threads to specify the number of threads (the same as --thread-num). The number of threads will be obtained from standard OpenMP environment variable OMP_NUM_THREADS if --thread-num or --threads is not specified.

### Changed

- Improved the speed and memory usage of --make-grm-xchr, and added an option --make-grm-xchr-part to reduce the memory usage further.
- Updated the options --update-sex, --update-ref-allele and --update-freq to be compatible with the new GRM functions.

### Fixed

- Fixed a bug in GSMR when there are multiple outcome variables.
- Fixed a bug in COJO when the standard error is extremely small.
- Fixed a bug of reporting "Illegal instruction" error for old CPUs (earlier than 2009).


## [1.91.3beta] - 2018-03-14 — [`51b86e8`](https://github.com/JianYang-Lab/GCTA/commit/51b86e8)

### Added

- Added a flag (--effect-plot) in GSMR for visualization.

### Changed

- Speeded up dominance GRM and added a flag --make-grm-d-part to partition the computation.

### Fixed

- Fixed a bug in REML, REML bivar, MLMA and LD when the number of threads (specified by --thread-num) is larger than 1.
- Fixed a bug in COJO for the X chromosome when there is no gender information in the .fam file.
- Fixed a bug in mtCOJO.

### Notes

- Redirected the log output to both screen and .log file.

### Commits

- [`51b86e8`](https://github.com/JianYang-Lab/GCTA/commit/51b86e8) — Merge mtcojo plot v1.91.3


## [1.91.2beta] - 2018-02-02 — [`c50935e`](https://github.com/JianYang-Lab/GCTA/commit/c50935e)

### Added

- Added a new module GSMR.
- Added a flag (--mbfile) to read multiple PLINK binary files for GSMR and mtCOJO.

### Fixed

- Fixed a bug in SBLUP, and improved the speed by 40%.
- Fixed a bug in MLMA when dealing with the individuals' ID.
- Fixed unreadable characters in the output of some computer clusters.

### Commits

- [`c50935e`](https://github.com/JianYang-Lab/GCTA/commit/c50935e) — Fix SBLUP NA problem, v1.91.2


## [1.91.1beta] - 2017-11-25

### Changed

- Changed to use the shared library glibc avoid segmentation fault in higher versions of Linux kernel.

### Fixed

- Fixed a bug in --mtcojo.
- Fixed a memory issue in REML analysis and improved the speed by 3 times in the Linux version.


## [1.91.0beta] - 2017-10-21

### Added

- Added a new module mtCOJO

### Fixed

- Fixed an issue of file path in the Windows version


## [1.90.2beta] - 2017-09-24

### Changed

- Removed the VC++ runtime dependency in the Windows version.

### Fixed

- Fixed a bug in --mlma-loco with the --mlma-no-adj-covar option.
- Fixed a bug in --make-grm-part when the sample size of one partition is larger than 69K.
- Fixed the performance issue in reading the PLINK .fam file.
- Fixed an issue with --autosome-num.


## [1.90.1beta] - 2017-09-13

### Changed

- Removed --grm-no-relative and added --grm-singleton to get singleton subjects from a sample.

### Fixed

- Fixed a bug in estimating allele frequency in some occasions.
- Fixed a bug in computing a GRM occasionally in small sample.
- Fixed an issue in computing a GRM including rare variants.
- Fixed an issue to run Linux binary in the Linux subsystem on Windows 10.
- Fixed a memory issue in the Windows version.
- Fixed a memory issue in --make-bK.


## [1.90.0beta] - 2017-08-08

### Added

- Added a new option --make-grm-part to partition the GRM computation into a large of parts to facilitate the analysis in large data set such as the UK Biobank.
- Added the --grm-no-relative option to extract the GRM of a subset of individuals who do not have any close relative in the sample.
- Added an option --cojo-sblup to perform a summary-data-based BLUP prediction analysis.
- Added the Haseman-Elston regression analysis to estimate the SNP-based heritability for a trait and genetic correlation between traits.
- Added the Mac and Windows versions.

### Changed

- Improved the speed and memory usage of GRM computation by orders of magnitude.
- Improved the memory usage of the --grm-cutoff option.
- Improved the speed and memory usage of --freq by orders of magnitude.
- Improved the approximation accuracy of the COJO analysis.
- Improved the speed of the bivariate GREML analysis (5X faster than original version).

### Fixed

- Fixed the memory issue when the sample size exceeds 500K in some functions (e.g. bivariate GREML and reading the GRM in gz format).

### Notes

- Update the package dependencies to the latest, such as Intel MKL and Eigen. This improved the performance by ~40%.


## [1.26.0] - 2016-06-22

### Added

- Added a new module (GCTA-fastBAT) for a set- or gene-based association analysis using GWAS summary data.

### Fixed

- Fixed a bug in MLMA.


## [1.25.3] - 2016-04-27

### Fixed

- Fixed a memory leaking issue in --mlma


## [1.25.2] - 2015-12-22

### Added

- Added a new option (--mlma-subtract-grm) for MLMA-LOCO with large data sets.
- Added a new option (--make-grm-inbred) to compute GRM for an inbred population (e.g. inbred mice or crops).
- Added a new option (--recode-std) to output standardised SNP genotypes.

### Fixed

- Fst calculation has been changed to that based on a random model. The previous version was based on a fixed model. The difference is trivial for small Fst values but the random model has a good property that Fst is bounded at 1 for the most extreme allele frequency difference.

### Notes

- A much more memory-efficient version of MLMA.


## [1.25.1] - 2015-12-08

### Added

- Added an option --reml-bendV


## [1.25.0] - 2015-10-30

### Added

- Added an option to calculate an unbiased estimate of LD score for LDSC regression analysis (see gcta.freeforums.net/thread/177/gcta-lds-calculating-score-snp); Added an option to calculate multi-component LD score following Finucane et al. (2015 Nat Genet).
- Added options to extract or exclude a region.
- Add the --reml-bivar-no-constrain option to the bivariate GREML analysis.
- Add an option to select a fixed number of top associated SNPs (taking LD into account) from GWAS.

### Fixed

- Fixed a bug in --imp-rsq

### Notes

- We have implemented the Zaitlen et al. method in GCTA which allows to estimate SNP-based h2 in family data without having to remove related individuals.


## [1.24.7] - 2015-06-11

### Notes

- Mixed linear model association (MLMA) analysis with multiple GRMs
- Fst calculation
- Haseman-Elston regression
- LD score calculation


## [1.24.4] - 2014-07-29

### Fixed

- changed the syntax for the conditional and joint analysis; fixed memory leak issues in mixed linear model based association analysis and bivariate GREML analysis with multiple GRMs; enabled the function converting dosage data to PLINK best guess.


## [1.24.3] - 2014-06-05

### Added

- allows you to transform variance explained by all SNPs on the observed scale to that on the underlying scale in a bivariate analysis of a case-control study and a quantitative trait; pca

### Notes

- only the top eigenvalues will be printed out.


## [1.24.2] - 2014-03-12

### Fixed

- fixed a bug in the conditonal and joint analysis (GCTA-COJO) when doing a backward model selection.


## [1.24.1] - 2014-03-06

### Notes

- a small change that allows you to use "Rsq" or "Rsq_hat" as the header for the last column of the *.mlinfo file from MACH imputation.


## [1.24] - 2014-01-08

### Fixed

- fixed a bug in REML analysis as a result of a change made in v1.23 in transforming the estimate of genetic variance on the observed scale to that on the underlying scale; fixed a bug in GWAS simulation where the reported variance explained by a causal variant in the *.par file was incorrect.


## [1.23] - 2013-12-18

### Fixed

- changed --dosage-mach option and added a new option --dosage-mach-gz; fixed a bug in the --cojo-cond option when two SNPs are in very high LD and their allele frequencies are consistently higher in the reference sample than those in the discovery sample.


## [1.22] - 2013-10-31

### Fixed

- fixed a bug in the --dosage-mach option when used in combined with the --imput-rsq option.


## [1.21] - 2013-10-16

### Fixed

- fixed a bug in bivariate analysis including covariates; re-wrote the code for the option --dosage-mach; added a new option and changed syntax for the mixed linear model association analysis.


## [1.20] - 2013-08-23

### Added

- added a new module mixed linear model association analysis; fixed a few bugs; made a few improvements.


## [1.13] - 2013-03-19

### Fixed

- fixed a bug for the --make-grm-bin option.


## [1.11] - 2013-02-14

### Fixed

- fixed a bug for the --mgrm-bin option and added the option to test for genetic correlation = 0 or 1 in a bivariate analysis.


## [1.1] - 2013-02-10

### Notes

- a much faster version which allows multi-thread computing (new option --thread-num); added new options --make-grm-bin and --grm-bin to more efficiently read and write the GRM files.


## [1.04] - 2012-09-13

### Added

- added a new option to convert Minimac dosage data to PLINK binary PED format.


## [1.03] - 2012-08-30

### Fixed

- fixed a few bugs and added a new option to convert MACH dosage data to PLINK binary PED format.
- fixed 2 bugs.
- fixed a few bugs.

### Notes

- version 1.0 released!
- latest version (version 0.93.9) of source codes released.


## [0.93.9] - 2011-11-18

### Fixed

- modified the --dosage-mach option to be compatiable with the latest MACH version; fixed a bug with the option --ld.


## [0.93.8] - 2011-09-30

### Fixed

- fixed a bug for the option --grm-adj when the genotype data of some individuals are completely missing.


## [0.93.7] - 2011-09-10

### Fixed

- fixed a bug when the option --ibc is used in combined with the option --keep or --remove, which causes wrong IDs in the *.ibc fie; fixed a bug in --gxe option when there are missing values for the environmental factor; and modified the function for converting Illumina raw genotype data to that in PLINK format.


## [0.93.6] - 2011-08-28

### Fixed

- fixed a bug in the new option --reml-lrt which caused memory leak.


## [0.93.5] - 2011-08-26

### Added

- added an option to turn off the LRT and fixed a bug in the case that the IDs of multiple GRM files are not in the same order.


## [0.93.4] - 2011-08-15

### Added

- added a function to calculate the LRT for the REML analysis.


## [0.93.2] - 2011-07-18

### Fixed

- fixed a bug in the matrix bending subroutine.


## [0.93.1] - 2011-07-12

### Changed

- improved the efficiency of reading PLINK binary data.


## [0.93.0] - 2011-07-08

### Added

- added a subroutine to deal with the issue when the variance-covariance matrix V is negative-definite; changed the default number of maximum REML iterations from 30 to 100; changed the method of calculating the diagonal elements of GRM to be the same as that for the off-diagonal elements; modified REML procedure to allow some elements of the GRM to be missing (printing a warning on the screen in stead of an error message).
- added a few new functions, e.g. convert the raw genotype data into PLINK binary format.

### Changed

- modified the output of LD estimation and the input format of GWAS simulation

### Fixed

- fixed a bug in GWAS simulation.
- fixed a bug in a REML analysis, i.e. the estimate may be stuck at zero if the true parameter is very small.
- fixed a few bugs.
- fixed a bug in reading the PLINK FAM file.
- fixed a bug in transforming the estimate of variance explained by the SNPs on the observed scale to that on the underlying scale for a case-control study.
- fixed a bug in the estimation of LD and compiled the program statically (more compatible

### Notes

- source codes released.
- MacOS version released.
- first release.


---

[Unreleased]: https://github.com/peterk87/GCTA/compare/0dc78f0...HEAD
[1.95.3]: https://github.com/JianYang-Lab/GCTA/compare/4ee76be...0dc78f0
[1.95.2]: https://github.com/JianYang-Lab/GCTA/compare/69d42e2...4ee76be
[1.95.1]: https://github.com/JianYang-Lab/GCTA/compare/3dbff08...69d42e2
[1.95.0]: https://github.com/JianYang-Lab/GCTA/compare/3ecd226...3dbff08
[1.94.4]: https://github.com/JianYang-Lab/GCTA/compare/0bd4c0f...3ecd226
[1.94.1]: https://github.com/JianYang-Lab/GCTA/compare/0fbe2d3...0bd4c0f
[1.93.3beta]: https://github.com/JianYang-Lab/GCTA/compare/9571d1a...0fbe2d3
[1.93.2beta]: https://github.com/JianYang-Lab/GCTA/compare/51b86e8...9571d1a
[1.91.3beta]: https://github.com/JianYang-Lab/GCTA/compare/c50935e...51b86e8
[1.91.2beta]: https://github.com/JianYang-Lab/GCTA/commit/c50935e

Compare links for older releases are omitted until matching tags or version-bump commits are identified.

Upstream remote: `https://github.com/JianYang-Lab/GCTA.git`. Fork origin: `git@github.com:peterk87/GCTA.git`.
