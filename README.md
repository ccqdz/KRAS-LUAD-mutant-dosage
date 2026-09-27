# Orthogonal genomic mechanisms of KRAS mutant dosage in lung adenocarcinoma

Reproducibility repository for the manuscript:

**Orthogonal genomic mechanisms of KRAS mutant dosage define distinct biological states in lung adenocarcinoma**

Authors: Yuanyuan Zhou; Jingzhang Li (corresponding author)
Affiliation: Liuzhou People's Hospital
Correspondence: lzryljz@163.com

## Scope
This repository contains the consolidated R workflow and locked audit files used to reproduce the manuscript analyses and figures. The workflow separates KRAS mutant dosage into two locus-normalised axes: copy-number expansion (KRAS total copy number/ploidy) and allelic dominance (fraction of alleles mutated, FAM).

## Main analysis script
`R/00_KRAS_LUAD_MANUSCRIPT_FINAL_MASTER_V3_CENTERED.R`

Run from a project directory containing the required `input/` files:

```r
setwd("/path/to/KRAS_LUAD_FINAL_CLEAN")
source("R/00_KRAS_LUAD_MANUSCRIPT_FINAL_MASTER_V3_CENTERED.R")
```

The original local workflow used:

```r
setwd("/Users/googleinsulin/Documents/KRAS_LUAD_FINAL_CLEAN")
source("00_KRAS_LUAD_MANUSCRIPT_FINAL_MASTER_V3_CENTERED.R")
```

For a public repository, place the script under `R/` and adapt the working directory to your own environment.

## Required input files
See `INPUT_MANIFEST.tsv`. The minimum reproducible input set comprises 12 cohort/support files. Public or controlled-access source data are not redistributed here unless their licenses and access conditions permit it.

## Outputs
The workflow creates:
- `results/` analysis tables
- `figures/` manuscript figures
- `audit/FINAL_LOCKED_AUDIT.tsv`
- `audit/FINAL_MANUSCRIPT_PREFLIGHT.tsv`
- `audit/06_TCGA_RNA_PROVENANCE_AUDIT.tsv`
- `audit/SESSION_INFO.txt`

## Locked analysis conventions
- Primary TCGA mechanistic cohort: n=138.
- TCGA clinical endpoint cohort: n=135; adjusted complete-case cohort: n=132.
- TCGA progression endpoint is PFI, not PFS.
- TCGA RNA-evaluable cohort: n=63.
- MSK external allelic metric is a FACETS/ASCN-based proxy and is not mathematically identical to TCGA FAM.
- FAM >=0.66 was prespecified; categorical Kaplan-Meier results are descriptive, while adjusted Cox models are primary.
- CPTAC results are interpreted as proteogenomic concordance rather than causal mediation.

## Audit status
The final locked audit passed 32/32 checks. The exact audit artifacts are included under `audit/`.

## Software environment
See `audit/SESSION_INFO.txt` for the exact R and package versions used in the locked analysis.

## Data availability
TCGA-LUAD data are available through the NCI Genomic Data Commons and associated PanCancer Atlas resources. TCGA clinical endpoints derive from the TCGA Clinical Data Resource. OncoSG and CPTAC data are available through their respective study resources. MSK-IMPACT data should be obtained from the exact public or institutional source used in the manuscript and subject to its access terms. This repository does not redistribute controlled-access patient-level data.

## Citation
If you use this repository, cite the associated manuscript and the archived release DOI once assigned. Citation metadata are provided in `CITATION.cff`.

## License
No software license has been asserted automatically. Add the license approved by the authors/institution before public release.
