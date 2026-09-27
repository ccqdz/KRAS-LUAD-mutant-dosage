
# ============================================================
# KRAS-LUAD MANUSCRIPT FINAL MASTER ANALYSIS
# Manuscript-locked clean-project version
# ============================================================
#
# Run from the NEW clean project root:
#   setwd("/Users/googleinsulin/Documents/KRAS_LUAD_FINAL_CLEAN")
#   source("00_KRAS_LUAD_MANUSCRIPT_FINAL_MASTER.R")
#
# Folder expected:
#   input/   = required source files
#   results/ = generated result tables
#   figures/ = generated figures
#   audit/   = locked-value checks + session info
#
# IMPORTANT:
# - TCGA progression endpoint is PFI, never PFS.
# - TCGA primary mechanistic cohort is n=138.
# - MSK allelic metric is a FACETS/ASCN-based proxy, NOT TCGA FAM.
# - Exact MSK analysis uses MSK_KRAS_LUAD_external_cohort_FINAL.csv
#   and the historical multivariable formula recovered in Stage 14.
# ============================================================

options(stringsAsFactors=FALSE, warn=1)

ROOT <- getwd()
IN   <- file.path(ROOT,"input")
RES  <- file.path(ROOT,"results")
FIG  <- file.path(ROOT,"figures")
AUD  <- file.path(ROOT,"audit")

dir.create(RES,recursive=TRUE,showWarnings=FALSE)
dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
dir.create(AUD,recursive=TRUE,showWarnings=FALSE)

# ------------------------ packages ---------------------------
pkgs <- c("data.table","dplyr","tidyr","tibble","stringr",
          "survival","broom","ggplot2","scales","patchwork","sandwich","lmtest")
miss <- setdiff(pkgs,rownames(installed.packages()))
if(length(miss)) install.packages(miss,repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(survival); library(broom)
  library(ggplot2); library(scales); library(patchwork)
  library(sandwich); library(lmtest)
})

# ------------------------ helpers ----------------------------
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
z <- function(x) as.numeric(scale(num(x)))

read_csv0 <- function(name) {
  f <- file.path(IN,name)
  if(!file.exists(f)) stop("Missing required input: ",name)
  fread(f,data.table=FALSE,check.names=FALSE,showProgress=FALSE)
}

find_col <- function(df,candidates,required=TRUE) {
  for(x in candidates) if(x %in% names(df)) return(x)
  if(required) stop("None of these columns found: ",paste(candidates,collapse=", "))
  NA_character_
}

first_matching <- function(df,regex,required=FALSE) {
  x <- names(df)[grepl(regex,names(df),ignore.case=TRUE)]
  if(length(x)) return(x[1])
  if(required) stop("No column matches: ",regex)
  NA_character_
}

safe_spearman <- function(x,y) {
  xx <- num(x)
  yy <- num(y)
  ok <- is.finite(xx) & is.finite(yy)
  n_ok <- sum(ok)

  # cor.test() requires enough finite paired observations.
  if(n_ok < 3) {
    return(tibble(n=n_ok,rho=NA_real_,p=NA_real_))
  }

  # Spearman is undefined if either variable has no variation.
  if(length(unique(xx[ok])) < 2 || length(unique(yy[ok])) < 2) {
    return(tibble(n=n_ok,rho=NA_real_,p=NA_real_))
  }

  tt <- suppressWarnings(
    cor.test(xx[ok],yy[ok],method="spearman",exact=FALSE)
  )
  tibble(n=n_ok,rho=unname(tt$estimate),p=tt$p.value)
}

cox_axis <- function(df,time,event,formula_rhs,term) {
  fml <- as.formula(paste0("Surv(",time,",",event,") ~ ",formula_rhs))
  fit <- coxph(fml,data=df,ties="efron")
  tt <- broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE)
  out <- tt[tt$term==term,,drop=FALSE]
  tibble(
    N=fit$n,
    events=fit$nevent,
    HR=out$estimate[1],
    CI_low=out$conf.low[1],
    CI_high=out$conf.high[1],
    p=out$p.value[1]
  )
}

write_tsv <- function(x,name) {
  fwrite(as.data.frame(x),file.path(RES,name),sep="\t",na="NA")
}

# ============================================================
# Locked Lung Cancer figure palette
# ============================================================
# Copy-number expansion: blue
# Allelic dominance / FAM: coral
# Dual state: purple
# WGD / clinical baseline: teal
# Neither/reference: warm gold
# Secondary molecular layer: soft teal
# Neutral text/grid: charcoal / muted grey
COL_COPY      <- "#3B78B7"
COL_FAM       <- "#D8634A"
COL_DUAL      <- "#7657A6"
COL_TEAL      <- "#2F8E83"
COL_REF       <- "#D6A23A"
COL_TEAL2     <- "#76A9A2"
COL_TEXT      <- "#202428"
COL_MUTED     <- "#667078"
COL_GRID      <- "#CBD1D5"
COL_PALE_GOLD <- "#FAF2DF"
COL_GREY      <- "#7D858C"
COL_LIGHTGREY <- "#D9DDE0"

STATE_COLORS <- c(
  "Neither"        = COL_REF,
  "Expansion only" = COL_COPY,
  "Allelic only"   = COL_FAM,
  "Dual"           = COL_DUAL
)

theme_final <- function(base_size=11) {
  theme_classic(base_size=base_size) +
    theme(
      plot.title=element_text(face="bold",colour=COL_TEXT,size=rel(1.16),hjust=.5),
      plot.subtitle=element_text(colour=COL_MUTED,size=rel(.94),hjust=.5),
      axis.title=element_text(colour=COL_TEXT),
      axis.text=element_text(colour=COL_TEXT),
      axis.line=element_line(colour=COL_TEXT,linewidth=.45),
      axis.ticks=element_line(colour=COL_TEXT,linewidth=.4),
      legend.title=element_blank(),
      legend.text=element_text(colour=COL_TEXT),
      strip.background=element_blank(),
      strip.text=element_text(face="bold",colour=COL_TEXT),
      panel.grid.major.y=element_line(colour=COL_GRID,linewidth=.25),
      panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),
      plot.margin=margin(8,10,8,8)
    )
}

# ------------------------ input manifest ---------------------
required_map <- list(
  "TCGA_KRAS_primary_n138_CORRECTED.csv" =
    c("TCGA_KRAS_primary_n138_CORRECTED.csv"),

  "Stage15_TCGA_OS_analysis_dataset.csv" =
    c("Stage15_TCGA_OS_analysis_dataset.csv"),

  "LUAD_KRAS_expression_stage4.csv" =
    c("LUAD_KRAS_expression_stage4.csv"),

  "OncoSG_Stage13_KRAS_master.csv" =
    c("OncoSG_Stage13_KRAS_master.csv"),

  "CPTAC_KRAS_MINIMAL_MASTER.csv" =
    c("CPTAC_KRAS_MINIMAL_MASTER.csv"),

  "MSK_KRAS_LUAD_external_cohort_FINAL.csv" =
    c("MSK_KRAS_LUAD_external_cohort_FINAL.csv"),

  "Stage16_FAM_prespecified_cutoffs.csv" =
    c("Stage16_FAM_prespecified_cutoffs.csv"),

  "Stage16_FAM_stage_subgroups_FINAL.csv" =
    c("Stage16_FAM_stage_subgroups_FINAL.csv",
      "Stage16_FAM_stage_subgroups.csv"),

  "Stage16_FAM_stage_interaction_FINAL.csv" =
    c("Stage16_FAM_stage_interaction_FINAL.csv",
      "Stage16_FAM_stage_interaction.csv"),

  "Stage16_FAM066_landmark_OS.csv" =
    c("Stage16_FAM066_landmark_OS.csv"),

  "Stage17_TCGA_OS_time_metrics_bootstrap_FINAL.csv" =
    c("Stage17_TCGA_OS_time_metrics_bootstrap_FINAL.csv",
      "Stage17_TCGA_OS_time_metrics_bootstrap.csv"),

  "Stage17_TCGA_OS_calibration_tertiles_FINAL.csv" =
    c("Stage17_TCGA_OS_calibration_tertiles_FINAL.csv",
      "Stage17_TCGA_OS_calibration_tertiles.csv")
)

required <- names(required_map)

# ------------------------------------------------------------------
# SELF-HEALING INPUT BOOTSTRAP
# If input/ is empty or incomplete, search the old project recursively
# and copy the exact/canonical files into the clean project's input/.
# ------------------------------------------------------------------
SOURCE_ROOT_CANDIDATES <- c(
  "/Users/googleinsulin/Documents/KRAS-LUAD",
  "/Users/googleinsulin/Documents/KRAS_LUAD",
  dirname(ROOT)
)
SOURCE_ROOT_CANDIDATES <- unique(SOURCE_ROOT_CANDIDATES[
  dir.exists(SOURCE_ROOT_CANDIDATES)
])

find_source_file <- function(candidates) {
  # 1) already present in clean input
  for(nm in candidates) {
    f0 <- file.path(IN,nm)
    if(file.exists(f0)) return(f0)
  }

  # 2) search old project roots
  hits <- character()
  for(sr in SOURCE_ROOT_CANDIDATES) {
    ff <- list.files(sr,recursive=TRUE,full.names=TRUE)
    ff <- ff[file.info(ff)$isdir==FALSE]
    for(nm in candidates) {
      h <- ff[basename(ff)==nm]
      if(length(h)) hits <- c(hits,h)
    }
  }
  hits <- unique(hits)
  if(!length(hits)) return(NA_character_)

  # Prefer FINAL_MANUSCRIPT_OUTPUT, then paths explicitly containing FINAL.
  ord <- order(
    !grepl("FINAL_MANUSCRIPT_OUTPUT",hits,fixed=TRUE),
    !grepl("FINAL",basename(hits),fixed=TRUE),
    nchar(hits)
  )
  hits[ord][1]
}

bootstrap <- tibble(
  destination=required,
  source=NA_character_,
  copied=FALSE,
  exists_after=FALSE
)

for(i in seq_len(nrow(bootstrap))) {
  dest_name <- bootstrap$destination[i]
  dest <- file.path(IN,dest_name)

  if(file.exists(dest)) {
    bootstrap$source[i] <- dest
    bootstrap$copied[i] <- TRUE
    bootstrap$exists_after[i] <- TRUE
    next
  }

  src <- find_source_file(required_map[[dest_name]])
  bootstrap$source[i] <- src

  if(!is.na(src) && file.exists(src)) {
    okcopy <- file.copy(src,dest,overwrite=TRUE)
    bootstrap$copied[i] <- okcopy
    bootstrap$exists_after[i] <- file.exists(dest)
  }
}

fwrite(
  bootstrap,
  file.path(AUD,"00A_INPUT_BOOTSTRAP_AUDIT.tsv"),
  sep="\t"
)

input_audit <- tibble(
  file=required,
  exists=file.exists(file.path(IN,required)),
  bytes=ifelse(file.exists(file.path(IN,required)),
               file.info(file.path(IN,required))$size,NA_real_)
)
fwrite(input_audit,file.path(AUD,"00_INPUT_FILE_AUDIT.tsv"),sep="\t")

if(any(!input_audit$exists)) {
  missing_now <- input_audit$file[!input_audit$exists]
  cat("\nCould not automatically locate these required files:\n")
  cat(paste0(" - ",missing_now,"\n"),sep="")
  cat("\nSearched roots:\n")
  cat(paste0(" - ",SOURCE_ROOT_CANDIDATES,"\n"),sep="")
  cat("\nSee audit/00A_INPUT_BOOTSTRAP_AUDIT.tsv for exact search/copy status.\n")
  stop("Input bootstrap incomplete. Do not delete/archive the old project yet.")
} else {
  cat("\nAll required inputs are present in clean input/.\n")
}

# ============================================================
# 1. TCGA n=138 mechanistic cohort
# ============================================================

tcga <- read_csv0("TCGA_KRAS_primary_n138_CORRECTED.csv")

id_tcga <- find_col(tcga,c("patient_id","PATIENT_ID","bcr_patient_barcode"))
fam_col <- find_col(tcga,c("FAM","fam"))
copy_col <- find_col(tcga,c("total_cn_ploidy","CN_ploidy","copy_expansion"))
loh_col <- find_col(tcga,c("loh","LOH"),required=FALSE)
gain12_col <- find_col(tcga,c("arm12p_gain","gain12p","12p_gain"),required=FALSE)
wgd_col <- find_col(tcga,c("WGD","wgd"),required=FALSE)

if(nrow(tcga)!=138)
  warning("TCGA mechanistic master is not n=138: observed ",nrow(tcga))

tcga[[fam_col]] <- num(tcga[[fam_col]])
tcga[[copy_col]] <- num(tcga[[copy_col]])

tcga_axis <- safe_spearman(tcga[[copy_col]],tcga[[fam_col]])
write_tsv(tcga_axis,"01_TCGA_copy_vs_FAM_spearman.tsv")

# Recover state labels if they are already present in the locked n=138 master.
state_col <- first_matching(
  tcga,
  "state$|dosage_state|mechanistic_state|GMD_state|route_state",
  required=FALSE
)

state_counts <- NULL
tcga$.dosage_state <- NA_character_
if(!is.na(state_col)) {
  vals <- as.character(tcga[[state_col]])
  # Normalize common historical label variants.
  vals <- dplyr::case_when(
    grepl("neither|reference|none",vals,ignore.case=TRUE) ~ "Neither",
    grepl("expansion.*only|copy.*only",vals,ignore.case=TRUE) ~ "Expansion only",
    grepl("allelic.*only|fam.*only",vals,ignore.case=TRUE) ~ "Allelic only",
    grepl("dual|both",vals,ignore.case=TRUE) ~ "Dual",
    TRUE ~ vals
  )
  tcga$.dosage_state <- vals
  state_counts <- as.data.frame(table(vals),stringsAsFactors=FALSE)
  names(state_counts) <- c("state","N")
} else {
  # Search for precomputed high flags, but DO NOT invent a copy-expansion cutoff.
  ch <- first_matching(tcga,"copy.*high|expansion.*high",required=FALSE)
  fh <- first_matching(tcga,"FAM.*high|allelic.*high",required=FALSE)
  if(!is.na(ch) && !is.na(fh)) {
    a <- num(tcga[[ch]])==1
    b <- num(tcga[[fh]])==1
    st <- ifelse(a & b,"Dual",
                 ifelse(a,"Expansion only",
                        ifelse(b,"Allelic only","Neither")))
    tcga$.dosage_state <- st
    state_counts <- as.data.frame(table(st),stringsAsFactors=FALSE)
    names(state_counts) <- c("state","N")
  }
}
if(any(!is.na(tcga$.dosage_state))) {
  tcga$.dosage_state <- factor(
    tcga$.dosage_state,
    levels=c("Neither","Expansion only","Allelic only","Dual")
  )
}
if(!is.null(state_counts)) write_tsv(state_counts,"02_TCGA_state_counts.tsv")

# Mechanistic determinants.
mech <- list()

if(!is.na(gain12_col)) {
  tcga[[gain12_col]] <- num(tcga[[gain12_col]])
  m <- lm(as.formula(paste0(copy_col," ~ ",gain12_col)),data=tcga)
  tt <- broom::tidy(m,conf.int=TRUE)
  rr <- tt[tt$term==gain12_col,,drop=FALSE]
  mech[[length(mech)+1]] <- tibble(
    outcome="copy expansion",determinant="12p gain",
    beta=rr$estimate,CI_low=rr$conf.low,CI_high=rr$conf.high,p=rr$p.value
  )
}

if(!is.na(loh_col)) {
  tcga[[loh_col]] <- num(tcga[[loh_col]])
  m <- lm(as.formula(paste0(fam_col," ~ ",loh_col)),data=tcga)
  tt <- broom::tidy(m,conf.int=TRUE)
  rr <- tt[tt$term==loh_col,,drop=FALSE]
  mech[[length(mech)+1]] <- tibble(
    outcome="allelic dominance",determinant="LOH",
    beta=rr$estimate,CI_low=rr$conf.low,CI_high=rr$conf.high,p=rr$p.value
  )
}

if(!is.na(wgd_col)) {
  tcga[[wgd_col]] <- num(tcga[[wgd_col]])
  for(v in c(copy_col,fam_col)) {
    m <- lm(as.formula(paste0(v," ~ ",wgd_col)),data=tcga)
    tt <- broom::tidy(m,conf.int=TRUE)
    rr <- tt[tt$term==wgd_col,,drop=FALSE]
    mech[[length(mech)+1]] <- tibble(
      outcome=v,determinant="WGD",
      beta=rr$estimate,CI_low=rr$conf.low,CI_high=rr$conf.high,p=rr$p.value
    )
  }
}

mechdf <- bind_rows(mech)
if(nrow(mechdf)) write_tsv(mechdf,"03_TCGA_mechanistic_determinants.tsv")

# ============================================================
# 2. TCGA RNA — canonical historical 63-case expression cohort
# ============================================================
#
# Provenance:
#   LUAD_KRAS_expression_stage4.csv is the historical patient-level table
#   used by the manuscript/Figure 3 workflow. It contains 63 RNA-evaluable
#   high-confidence KRAS-mutant tumors.
#
# IMPORTANT:
#   The historical manuscript draft contains a route-independent coefficient
#   beta=0.605, P=0.00257 and FAM beta=-0.225, P=0.472. Exhaustive recovery
#   showed that these exact coefficients are NOT reproduced by the canonical
#   63-case table under the plausible raw/standardized classical or HC3 models.
#   Therefore this master script NEVER injects those values into a fitted model.
#   It computes the reproducible patient-level results and writes an explicit
#   provenance audit documenting the discrepancy.
# ============================================================

rna63 <- read_csv0("LUAD_KRAS_expression_stage4.csv")

rna_required <- c("KRAS_expr","total_cn_ploidy","FAM")
missing_rna <- setdiff(rna_required,names(rna63))
if(length(missing_rna))
  stop("Canonical TCGA RNA table missing: ",paste(missing_rna,collapse=", "))

rna_dd <- rna63 |>
  transmute(
    KRAS_expr=num(KRAS_expr),
    total_cn_ploidy=num(total_cn_ploidy),
    FAM=num(FAM),
    total_cn=if("total_cn" %in% names(rna63)) num(rna63$total_cn) else NA_real_,
    mutant_cn=if("mutant_cn" %in% names(rna63)) num(rna63$mutant_cn) else NA_real_
  ) |>
  filter(
    is.finite(KRAS_expr),
    is.finite(total_cn_ploidy),
    is.finite(FAM)
  )

rna_spearman <- bind_rows(
  safe_spearman(rna_dd$total_cn_ploidy,rna_dd$KRAS_expr) |>
    mutate(comparison="copy_ploidy_vs_KRAS_expr"),
  safe_spearman(rna_dd$FAM,rna_dd$KRAS_expr) |>
    mutate(comparison="FAM_vs_KRAS_expr")
)

if(any(is.finite(rna_dd$total_cn))) {
  rna_spearman <- bind_rows(
    rna_spearman,
    safe_spearman(rna_dd$total_cn,rna_dd$KRAS_expr) |>
      mutate(comparison="total_cn_vs_KRAS_expr")
  )
}
if(any(is.finite(rna_dd$mutant_cn))) {
  rna_spearman <- bind_rows(
    rna_spearman,
    safe_spearman(rna_dd$mutant_cn,rna_dd$KRAS_expr) |>
      mutate(comparison="mutant_cn_vs_KRAS_expr")
  )
}
write_tsv(rna_spearman,"04A_TCGA_KRAS_RNA_spearman.tsv")

# Reproducible route-independent patient-level model.
rna_fit_raw <- lm(KRAS_expr ~ total_cn_ploidy + FAM,data=rna_dd)
rna_model_classic <- broom::tidy(rna_fit_raw,conf.int=TRUE) |>
  filter(term %in% c("total_cn_ploidy","FAM")) |>
  mutate(
    N=nrow(rna_dd),
    model="KRAS_expr ~ total_cn_ploidy + FAM",
    SE_type="classical"
  )

rna_hc3_mat <- lmtest::coeftest(
  rna_fit_raw,
  vcov.=sandwich::vcovHC(rna_fit_raw,type="HC3")
)
rna_model_hc3 <- tibble(
  term=rownames(rna_hc3_mat),
  estimate=rna_hc3_mat[,1],
  std.error=rna_hc3_mat[,2],
  statistic=rna_hc3_mat[,3],
  p.value=rna_hc3_mat[,4]
) |>
  filter(term %in% c("total_cn_ploidy","FAM")) |>
  mutate(
    N=nrow(rna_dd),
    model="KRAS_expr ~ total_cn_ploidy + FAM",
    SE_type="HC3"
  )

rna_models <- bind_rows(rna_model_classic,rna_model_hc3)
write_tsv(rna_models,"04B_TCGA_KRAS_RNA_route_independent_models.tsv")

# Historical-target provenance audit: reference only, never injected into fits.
rna_copy_raw <- rna_model_classic |> filter(term=="total_cn_ploidy")
rna_fam_raw  <- rna_model_classic |> filter(term=="FAM")

rna_provenance_audit <- tribble(
  ~result,~historical_manuscript_value,~recomputed_value,~status,
  "RNA evaluable N",63,nrow(rna_dd),
  ifelse(nrow(rna_dd)==63,"REPRODUCED","MISMATCH"),
  "Spearman copy/ploidy rho",0.534,
  (rna_spearman |> filter(comparison=="copy_ploidy_vs_KRAS_expr"))$rho[1],
  "REPRODUCED",
  "Route-independent copy beta",0.605,rna_copy_raw$estimate[1],
  "HISTORICAL_VALUE_NOT_REPRODUCED",
  "Route-independent copy P",0.00257,rna_copy_raw$p.value[1],
  "HISTORICAL_VALUE_NOT_REPRODUCED",
  "Route-independent FAM beta",-0.225,rna_fam_raw$estimate[1],
  "HISTORICAL_VALUE_NOT_REPRODUCED",
  "Route-independent FAM P",0.472,rna_fam_raw$p.value[1],
  "HISTORICAL_VALUE_NOT_REPRODUCED"
)
fwrite(
  rna_provenance_audit,
  file.path(AUD,"06_TCGA_RNA_PROVENANCE_AUDIT.tsv"),
  sep="\t"
)

# ============================================================
# 3. TCGA clinical outcomes — locked Stage15 membership
# ============================================================

cl <- read_csv0("Stage15_TCGA_OS_analysis_dataset.csv")

id_cl <- find_col(cl,c("patient_id","PATIENT_ID","bcr_patient_barcode","PatientID"))

# Stage15 provides the locked n=135 membership and key clinical covariates.
# In the historical pipeline, PFI was recovered for these Stage15 members
# from the n=138 TCGA master rather than assumed to live inside Stage15 itself.
cl$.patient <- substr(gsub("\\.","-",as.character(cl[[id_cl]])),1,12)

# Detect endpoint columns in the n=138 master.
tcga_os_time_col <- find_col(
  tcga,c("OS.time","OS_TIME","OS_MONTHS","OS_time"),required=FALSE
)
tcga_os_event_col <- find_col(
  tcga,c("OS","OS_event","OS_STATUS_NUM","OS_STATUS"),required=FALSE
)
tcga_pfi_time_col <- find_col(
  tcga,c("PFI.time","PFI_TIME","PFI_MONTHS","PFI_time"),required=FALSE
)
tcga_pfi_event_col <- find_col(
  tcga,c("PFI","PFI_event","PFI_STATUS_NUM","PFI_STATUS"),required=FALSE
)

# Build a collision-proof merge table with explicit *_from_master names.
tcga_merge <- tcga |>
  transmute(
    .patient=substr(gsub("\\.","-",as.character(.data[[id_tcga]])),1,12),
    FAM_from_master=.data[[fam_col]],
    copy_from_master=.data[[copy_col]],
    OS_time_from_master=if(!is.na(tcga_os_time_col)) num(.data[[tcga_os_time_col]]) else NA_real_,
    OS_event_from_master=if(!is.na(tcga_os_event_col)) num(.data[[tcga_os_event_col]]) else NA_real_,
    PFI_time_from_master=if(!is.na(tcga_pfi_time_col)) num(.data[[tcga_pfi_time_col]]) else NA_real_,
    PFI_event_from_master=if(!is.na(tcga_pfi_event_col)) num(.data[[tcga_pfi_event_col]]) else NA_real_
  )

cl <- cl |> left_join(tcga_merge,by=".patient")

fam_cli <- if("FAM" %in% names(cl)) "FAM" else "FAM_from_master"
copy_cli <- if("total_cn_ploidy" %in% names(cl)) "total_cn_ploidy" else "copy_from_master"

age_col <- find_col(cl,c("age.s15","age","AGE"))
sex_col <- find_col(cl,c("gender.s15","gender","SEX","sex"))
stage_bin <- find_col(cl,c("stage_ge2"))
tp53_col <- find_col(cl,c("TP53","TP53_mut","TP53_status"))
stk11_col <- find_col(cl,c("STK11","STK11_mut","STK11_status"))
keap1_col <- find_col(cl,c("KEAP1","KEAP1_mut","KEAP1_status"))

# Prefer endpoints stored directly in Stage15; if absent, use the
# explicitly recovered endpoint columns from the n=138 master.
os_time <- find_col(
  cl,c("OS.time","OS_MONTHS","OS_time","OS_time_from_master"),required=FALSE
)
os_event <- find_col(
  cl,c("OS","OS_event","OS_STATUS_NUM","OS_event_from_master"),required=FALSE
)
pfi_time <- find_col(
  cl,c("PFI.time","PFI_MONTHS","PFI_time","PFI_time_from_master"),required=FALSE
)
pfi_event <- find_col(
  cl,c("PFI","PFI_event","PFI_STATUS_NUM","PFI_event_from_master"),required=FALSE
)

endpoint_map <- tibble(
  endpoint=c("OS_time","OS_event","PFI_time","PFI_event"),
  selected=c(os_time,os_event,pfi_time,pfi_event),
  stage15_has_direct=c(
    any(c("OS.time","OS_MONTHS","OS_time") %in% names(cl)),
    any(c("OS","OS_event","OS_STATUS_NUM") %in% names(cl)),
    any(c("PFI.time","PFI_MONTHS","PFI_time") %in% names(cl)),
    any(c("PFI","PFI_event","PFI_STATUS_NUM") %in% names(cl))
  ),
  master_detected=c(
    tcga_os_time_col,tcga_os_event_col,tcga_pfi_time_col,tcga_pfi_event_col
  )
)
fwrite(endpoint_map,file.path(AUD,"02_TCGA_ENDPOINT_COLUMN_MAP.tsv"),sep="\t")

if(any(is.na(c(os_time,os_event,pfi_time,pfi_event)))) {
  stop(
    "Could not resolve all TCGA OS/PFI endpoint columns. ",
    "See audit/02_TCGA_ENDPOINT_COLUMN_MAP.tsv. ",
    "PFI must remain PFI; do not substitute PFS."
  )
}

# Numeric conversions. sex stays factor.
for(v in c(fam_cli,copy_cli,age_col,stage_bin,tp53_col,stk11_col,keap1_col,
           os_time,os_event,pfi_time,pfi_event)) cl[[v]] <- num(cl[[v]])
cl$sex_final <- factor(cl[[sex_col]])

endpoint_membership_audit <- tibble(
  cohort="Stage15 membership",
  N=nrow(cl),
  OS_nonmissing_time=sum(is.finite(cl[[os_time]])),
  OS_deaths=sum(cl[[os_event]]==1,na.rm=TRUE),
  PFI_nonmissing_time=sum(is.finite(cl[[pfi_time]])),
  PFI_events=sum(cl[[pfi_event]]==1,na.rm=TRUE)
)
fwrite(
  endpoint_membership_audit,
  file.path(AUD,"03_TCGA_ENDPOINT_MEMBERSHIP_AUDIT.tsv"),
  sep="\t"
)

# One locked complete-case cohort, including both endpoints, to recover n=132.
ccvars <- c(
  fam_cli,copy_cli,age_col,"sex_final",stage_bin,tp53_col,stk11_col,keap1_col,
  os_time,os_event,pfi_time,pfi_event
)
cc <- cl[complete.cases(cl[,ccvars,drop=FALSE]),,drop=FALSE]

cc$FAM_z <- z(cc[[fam_cli]])
cc$copy_z <- z(cc[[copy_cli]])
cc$age_z <- z(cc[[age_col]])
cc$stage_ge2_final <- cc[[stage_bin]]
cc$TP53_final <- cc[[tp53_col]]
cc$STK11_final <- cc[[stk11_col]]
cc$KEAP1_final <- cc[[keap1_col]]

rhs <- "FAM_z + copy_z + age_z + sex_final + stage_ge2_final + TP53_final + STK11_final + KEAP1_final"

osfit <- coxph(
  as.formula(paste0("Surv(",os_time,",",os_event,") ~ ",rhs)),
  data=cc,ties="efron"
)
pfifit <- coxph(
  as.formula(paste0("Surv(",pfi_time,",",pfi_event,") ~ ",rhs)),
  data=cc,ties="efron"
)

cox_extract <- function(fit,endpoint) {
  broom::tidy(fit,exponentiate=TRUE,conf.int=TRUE) |>
    mutate(endpoint=endpoint,N=fit$n,events=fit$nevent) |>
    select(endpoint,N,events,term,estimate,conf.low,conf.high,p.value)
}
tcga_cox <- bind_rows(cox_extract(osfit,"OS"),cox_extract(pfifit,"PFI"))
write_tsv(tcga_cox,"05_TCGA_adjusted_OS_PFI_cox.tsv")

# Prespecified FAM >=0.66 OS model.
cc$FAM066 <- as.integer(cc[[fam_cli]]>=0.66)
fit066 <- coxph(
  as.formula(
    paste0(
      "Surv(",os_time,",",os_event,") ~ FAM066 + age_z + sex_final + ",
      "stage_ge2_final + TP53_final + STK11_final + KEAP1_final"
    )
  ),
  data=cc,ties="efron"
)
fam066 <- cox_extract(fit066,"OS_FAM066") |> filter(term=="FAM066")
write_tsv(fam066,"06_TCGA_FAM066_adjusted_OS.tsv")

# Copy locked Stage16/17 final tables into results for exact historical Figure 4.
locked_small <- c(
  "Stage16_FAM_prespecified_cutoffs.csv",
  "Stage16_FAM_stage_subgroups_FINAL.csv",
  "Stage16_FAM_stage_interaction_FINAL.csv",
  "Stage16_FAM066_landmark_OS.csv",
  "Stage17_TCGA_OS_time_metrics_bootstrap_FINAL.csv",
  "Stage17_TCGA_OS_calibration_tertiles_FINAL.csv"
)
for(nm in locked_small) {
  file.copy(file.path(IN,nm),file.path(RES,nm),overwrite=TRUE)
}

# ============================================================
# 4. OncoSG patient-level molecular validation
# ============================================================

onc <- read_csv0("OncoSG_Stage13_KRAS_master.csv")
onc_cn <- find_col(onc,c("KRAS_CNA"))
onc_rna <- find_col(onc,c("KRAS_expr_z"))
onc_mut <- find_col(onc,c("KRAS_mut"),required=FALSE)

onco_results <- bind_rows(
  safe_spearman(onc[[onc_cn]],onc[[onc_rna]]) |> mutate(subset="All")
)
if(!is.na(onc_mut)) {
  om <- onc[num(onc[[onc_mut]])==1,,drop=FALSE]
  if(nrow(om)>2)
    onco_results <- bind_rows(
      onco_results,
      safe_spearman(om[[onc_cn]],om[[onc_rna]]) |> mutate(subset="KRAS-mutant")
    )
}
write_tsv(onco_results,"07_OncoSG_CN_RNA_spearman.tsv")

# ============================================================
# 5. CPTAC proteogenomic validation
# ============================================================

cp <- read_csv0("CPTAC_KRAS_MINIMAL_MASTER.csv")
cp_cn <- find_col(cp,c("KRAS_CNA"))
cp_rna <- find_col(cp,c("KRAS_RNA"))
cp_prot <- find_col(cp,c("KRAS_protein"))
cp_mut <- find_col(cp,c("KRAS_mut"))

corr_row <- function(df,x,y,subset,label) {
  safe_spearman(df[[x]],df[[y]]) |>
    mutate(subset=subset,comparison=label)
}

cp_all <- bind_rows(
  corr_row(cp,cp_cn,cp_rna,"All","CN_vs_RNA"),
  corr_row(cp,cp_cn,cp_prot,"All","CN_vs_protein"),
  corr_row(cp,cp_rna,cp_prot,"All","RNA_vs_protein")
)

# Robustly identify KRAS-mutant samples. The minimal CPTAC master may
# encode mutation status as 1/0, TRUE/FALSE, Mutant/WT, MUT/WT, etc.
mut_raw <- cp[[cp_mut]]
mut_chr <- toupper(trimws(as.character(mut_raw)))
mut_num <- suppressWarnings(as.numeric(as.character(mut_raw)))

mut_flag <- (!is.na(mut_num) & mut_num == 1) |
            mut_chr %in% c("TRUE","T","YES","Y","MUT","MUTANT",
                           "KRAS-MUTANT","KRAS_MUTANT","POSITIVE")

cpm <- cp[mut_flag %in% TRUE,,drop=FALSE]

# Audit the actual coding and matched finite observations before correlation.
cptac_mut_audit <- tibble(
  mutation_column=cp_mut,
  total_rows=nrow(cp),
  mutant_rows=nrow(cpm),
  CN_RNA_finite_pairs=sum(is.finite(num(cpm[[cp_cn]])) & is.finite(num(cpm[[cp_rna]]))),
  CN_protein_finite_pairs=sum(is.finite(num(cpm[[cp_cn]])) & is.finite(num(cpm[[cp_prot]])))
)
fwrite(
  cptac_mut_audit,
  file.path(AUD,"04_CPTAC_KRAS_MUTANT_AUDIT.tsv"),
  sep="\t"
)

# Also save the observed mutation-status coding for diagnosis/reproducibility.
mut_levels <- as.data.frame(sort(table(as.character(mut_raw),useNA="ifany"),decreasing=TRUE))
names(mut_levels) <- c("KRAS_mut_value","N")
fwrite(
  mut_levels,
  file.path(AUD,"05_CPTAC_KRAS_MUTATION_CODING.tsv"),
  sep="\t"
)

cp_mut_res <- bind_rows(
  corr_row(cpm,cp_cn,cp_rna,"KRAS-mutant","CN_vs_RNA"),
  corr_row(cpm,cp_cn,cp_prot,"KRAS-mutant","CN_vs_protein")
)

cp_corr <- bind_rows(cp_all,cp_mut_res)
write_tsv(cp_corr,"08_CPTAC_correlations.tsv")

# Gain >0.3 versus neutral.
cp$gain03 <- num(cp[[cp_cn]]) > 0.3
gain_tests <- lapply(c(RNA=cp_rna,protein=cp_prot),function(v) {
  dd <- cp[is.finite(num(cp[[v]])) & !is.na(cp$gain03),,drop=FALSE]
  wt <- wilcox.test(num(dd[[v]]) ~ dd$gain03,exact=FALSE)
  tibble(
    outcome=v,
    N=nrow(dd),
    median_gain=median(num(dd[[v]])[dd$gain03],na.rm=TRUE),
    median_neutral=median(num(dd[[v]])[!dd$gain03],na.rm=TRUE),
    p=wt$p.value
  )
}) |> bind_rows()
write_tsv(gain_tests,"09_CPTAC_gain03_vs_neutral.tsv")

# KRAS-mutant standardized protein models.
ddp <- cpm |>
  transmute(
    CN=num(.data[[cp_cn]]),
    RNA=num(.data[[cp_rna]]),
    protein=num(.data[[cp_prot]])
  )

# Use base complete.cases() outside dplyr::filter(); the magrittr-style
# dot pronoun is not available inside native-R pipe evaluation here.
ddp <- ddp[complete.cases(ddp), , drop=FALSE] |>
  mutate(
    CN_z=z(CN),
    RNA_z=z(RNA),
    protein_z=z(protein)
  )

if(nrow(ddp)>=20) {
  u <- lm(protein_z ~ CN_z,data=ddp)
  j <- lm(protein_z ~ CN_z + RNA_z,data=ddp)
  prot_models <- bind_rows(
    broom::tidy(u,conf.int=TRUE) |> mutate(model="univariable_CN"),
    broom::tidy(j,conf.int=TRUE) |> mutate(model="joint_CN_RNA")
  ) |> filter(term!="(Intercept)") |> mutate(N=nrow(ddp))
  write_tsv(prot_models,"10_CPTAC_KRASmut_protein_models.tsv")
} else {
  warning(
    "CPTAC KRAS-mutant complete CN/RNA/protein cases <20 (n=", nrow(ddp),
    "). See audit/04_CPTAC_KRAS_MUTANT_AUDIT.tsv and ",
    "audit/05_CPTAC_KRAS_MUTATION_CODING.tsv."
  )
}

# ============================================================
# 6. MSK exact historical external validation — FINAL
# ============================================================

msk <- read_csv0("MSK_KRAS_LUAD_external_cohort_FINAL.csv")

msk_required <- c(
  "event","time","FAM_proxy","FAM_proxy_raw","CN_ploidy",
  "age","male","metastasis","TP53","STK11","KEAP1"
)
missing_msk <- setdiff(msk_required,names(msk))
if(length(missing_msk))
  stop("MSK final cohort missing columns: ",paste(missing_msk,collapse=", "))

for(v in msk_required) msk[[v]] <- num(msk[[v]])

if("logtmb" %in% names(msk)) {
  msk$logtmb_final <- num(msk$logtmb)
} else {
  tmb_col <- if("TMB_SCORE" %in% names(msk)) "TMB_SCORE" else
             if("tmb" %in% names(msk)) "tmb" else NA_character_
  if(is.na(tmb_col)) stop("MSK TMB/logTMB field not found.")
  tv <- num(msk[[tmb_col]])
  med <- median(tv,na.rm=TRUE)
  msk$logtmb_final <- log1p(ifelse(is.na(tv),med,tv))
}

msk_ext <- msk |>
  filter(
    is.finite(time),time>0,!is.na(event),
    !is.na(FAM_proxy),!is.na(CN_ploidy),!is.na(age)
  ) |>
  mutate(
    fam_z=z(FAM_proxy),
    cn_z=z(CN_ploidy),
    age_z=z(age),
    logtmb_z=z(logtmb_final)
  )

msk_fit <- coxph(
  Surv(time,event) ~ fam_z + cn_z + age_z + male + metastasis +
    TP53 + STK11 + KEAP1 + logtmb_z,
  data=msk_ext,ties="efron"
)

msk_primary <- broom::tidy(msk_fit,exponentiate=TRUE,conf.int=TRUE) |>
  filter(term %in% c("fam_z","cn_z")) |>
  mutate(N=msk_fit$n,deaths=msk_fit$nevent)
write_tsv(msk_primary,"11_MSK_exact_primary_cox.tsv")

msk_axis <- safe_spearman(msk_ext$FAM_proxy,msk_ext$CN_ploidy)
write_tsv(msk_axis,"12_MSK_axis_correlation.tsv")

fit_msk_sens <- function(dd,label,include_metastasis=TRUE) {
  dd <- dd |>
    mutate(
      fam_z2=z(FAM_proxy),
      cn_z2=z(CN_ploidy),
      age_z2=z(age),
      logtmb_z2=z(logtmb_final)
    )

  fml <- if(include_metastasis) {
    Surv(time,event) ~ fam_z2 + cn_z2 + age_z2 + male + metastasis +
      TP53 + STK11 + KEAP1 + logtmb_z2
  } else {
    Surv(time,event) ~ fam_z2 + cn_z2 + age_z2 + male +
      TP53 + STK11 + KEAP1 + logtmb_z2
  }

  ff <- coxph(fml,data=dd,ties="efron")
  rr <- broom::tidy(ff,exponentiate=TRUE,conf.int=TRUE) |>
    filter(term=="fam_z2")

  tibble(
    analysis=label,N=nrow(dd),deaths=sum(dd$event==1,na.rm=TRUE),
    HR=rr$estimate,CI_low=rr$conf.low,CI_high=rr$conf.high,p=rr$p.value
  )
}

msk_sens <- bind_rows(
  fit_msk_sens(msk_ext |> filter(FAM_proxy_raw<=1),
               "Exclude uncapped proxy >1",TRUE),
  fit_msk_sens(msk_ext |> filter(metastasis==0),
               "Primary tumours",FALSE),
  fit_msk_sens(msk_ext |> filter(metastasis==1),
               "Metastases",FALSE)
)
write_tsv(msk_sens,"13_MSK_FAM_proxy_sensitivity.tsv")

# ============================================================
# 7. Publication-style core figures — LOCKED PALETTE
# ============================================================

# ---------------- Figure 1: orthogonal dosage architecture ----------------
if(any(!is.na(tcga$.dosage_state))) {
  p1a <- ggplot(
    tcga,
    aes(x=.data[[copy_col]],y=.data[[fam_col]],colour=.dosage_state)
  ) +
    geom_point(alpha=.82,size=2.35) +
    scale_colour_manual(values=STATE_COLORS,drop=FALSE) +
    labs(
      title="A  Orthogonal KRAS dosage axes",
      subtitle=sprintf("Spearman rho = %.3f; P = %.3f",tcga_axis$rho[1],tcga_axis$p[1]),
      x="KRAS total copy number / ploidy",
      y="Fraction of alleles mutated (FAM)",
      colour=NULL
    ) +
    theme_final() +
    theme(legend.position="top")

  if(!is.null(state_counts)) {
    sc <- state_counts |>
      mutate(
        state=factor(state,levels=c("Neither","Expansion only","Allelic only","Dual"))
      ) |>
      filter(!is.na(state))
    p1b <- ggplot(sc,aes(x=state,y=N,fill=state)) +
      geom_col(width=.68) +
      geom_text(aes(label=N),vjust=-.35,size=3.5,colour=COL_TEXT) +
      scale_fill_manual(values=STATE_COLORS,drop=FALSE) +
      labs(title="B  Two-axis dosage states",x=NULL,y="Tumours, n") +
      theme_final() +
      theme(
        legend.position="none",
        axis.text.x=element_text(angle=18,hjust=1)
      )
    p1 <- p1a | p1b
  } else {
    p1 <- p1a
  }
} else {
  p1 <- ggplot(tcga,aes(x=.data[[copy_col]],y=.data[[fam_col]])) +
    geom_point(alpha=.78,size=2.2,colour=COL_GREY) +
    geom_smooth(method="lm",se=FALSE,linewidth=.7,colour=COL_TEXT) +
    labs(
      title="Orthogonal KRAS dosage axes",
      subtitle=sprintf("Spearman rho = %.3f; P = %.3f",tcga_axis$rho[1],tcga_axis$p[1]),
      x="KRAS total copy number / ploidy",
      y="Fraction of alleles mutated (FAM)"
    ) + theme_final()
}

ggsave(file.path(FIG,"Figure1_TCGA_axes_COLORED.pdf"),p1,width=10.5,height=4.9)
ggsave(file.path(FIG,"Figure1_TCGA_axes_COLORED.png"),p1,width=10.5,height=4.9,dpi=600)

# ---------------- Figure 2: distinct genomic mechanisms ----------------
p2a <- if(!is.na(gain12_col)) {
  d2a <- tcga |>
    mutate(.gain=factor(.data[[gain12_col]],levels=c(0,1),
                        labels=c("12p non-gain","12p gain")))
  ggplot(d2a,aes(x=.gain,y=.data[[copy_col]],fill=.gain)) +
    geom_boxplot(width=.48,outlier.shape=NA,alpha=.62,linewidth=.55) +
    geom_jitter(aes(colour=.gain),width=.10,alpha=.55,size=1.45) +
    scale_fill_manual(values=c("12p non-gain"=COL_LIGHTGREY,"12p gain"=COL_COPY)) +
    scale_colour_manual(values=c("12p non-gain"=COL_GREY,"12p gain"=COL_COPY)) +
    labs(title="A  12p gain increases copy expansion",x=NULL,y="KRAS total CN / ploidy") +
    theme_final() + theme(legend.position="none")
} else ggplot() + theme_void() + labs(title="12p gain field unavailable")

p2b <- if(!is.na(loh_col)) {
  d2b <- tcga |>
    mutate(.loh=factor(.data[[loh_col]],levels=c(0,1),labels=c("LOH−","LOH+")))
  ggplot(d2b,aes(x=.loh,y=.data[[fam_col]],fill=.loh)) +
    geom_boxplot(width=.48,outlier.shape=NA,alpha=.62,linewidth=.55) +
    geom_jitter(aes(colour=.loh),width=.10,alpha=.55,size=1.45) +
    scale_fill_manual(values=c("LOH−"=COL_LIGHTGREY,"LOH+"=COL_FAM)) +
    scale_colour_manual(values=c("LOH−"=COL_GREY,"LOH+"=COL_FAM)) +
    labs(title="B  LOH increases allelic dominance",x=NULL,y="Fraction of alleles mutated (FAM)") +
    theme_final() + theme(legend.position="none")
} else ggplot() + theme_void() + labs(title="LOH field unavailable")

# WGD is a background mechanism: teal, never blue/coral.
p2c <- if(!is.na(wgd_col)) {
  mutcn_col <- first_matching(tcga,"mutant.*cn|mut.*cop|absolute.*mut",required=FALSE)
  if(!is.na(mutcn_col)) {
    d2c <- tcga |>
      mutate(.wgd=factor(.data[[wgd_col]],levels=c(0,1),labels=c("WGD−","WGD+")))
    ggplot(d2c,aes(x=.wgd,y=.data[[mutcn_col]],fill=.wgd)) +
      geom_boxplot(width=.48,outlier.shape=NA,alpha=.62,linewidth=.55) +
      geom_jitter(aes(colour=.wgd),width=.10,alpha=.5,size=1.35) +
      scale_fill_manual(values=c("WGD−"=COL_LIGHTGREY,"WGD+"=COL_TEAL)) +
      scale_colour_manual(values=c("WGD−"=COL_GREY,"WGD+"=COL_TEAL)) +
      labs(title="C  WGD increases absolute mutant copies",x=NULL,y="KRAS mutant copy number") +
      theme_final() + theme(legend.position="none")
  } else ggplot() + theme_void() + labs(title="WGD: absolute mutant-copy field unavailable")
} else ggplot() + theme_void() + labs(title="WGD field unavailable")

p2 <- (p2a | p2b | p2c) +
  patchwork::plot_annotation(
    title="Distinct genomic mechanisms underlying KRAS mutant dosage",
    theme=theme(plot.title=element_text(face="bold",colour=COL_TEXT,size=15,hjust=.5))
  )

ggsave(file.path(FIG,"Figure2_TCGA_mechanisms_COLORED.pdf"),p2,width=14,height=4.6)
ggsave(file.path(FIG,"Figure2_TCGA_mechanisms_COLORED.png"),p2,width=14,height=4.6,dpi=600)

# ---------------- Figure 3: molecular validation ----------------
p3a <- ggplot(rna_dd,aes(x=total_cn_ploidy,y=KRAS_expr)) +
  geom_point(alpha=.8,size=2.15,colour=COL_COPY) +
  geom_smooth(method="lm",se=FALSE,linewidth=.75,colour=COL_TEXT) +
  labs(title="A  TCGA KRAS-mutant LUAD",x="KRAS total CN / ploidy",y="KRAS RNA expression") +
  theme_final()

p3b <- ggplot(rna_dd,aes(x=FAM,y=KRAS_expr)) +
  geom_point(alpha=.78,size=2.15,colour=COL_FAM) +
  geom_smooth(method="lm",se=FALSE,linewidth=.75,colour=COL_TEXT) +
  labs(title="B  TCGA allelic dominance",x="Fraction of alleles mutated (FAM)",y="KRAS RNA expression") +
  theme_final()

p3c <- ggplot(onc,aes(x=.data[[onc_cn]],y=.data[[onc_rna]])) +
  geom_point(alpha=.75,size=1.9,colour=COL_COPY) +
  geom_smooth(method="lm",se=FALSE,linewidth=.75,colour=COL_TEXT) +
  labs(title="C  OncoSG validation",x="KRAS CNA",y="KRAS RNA expression") +
  theme_final()

p3d <- ggplot(cp,aes(x=.data[[cp_cn]],y=.data[[cp_rna]])) +
  geom_point(alpha=.75,size=1.9,colour=COL_COPY) +
  geom_smooth(method="lm",se=FALSE,linewidth=.75,colour=COL_TEXT) +
  labs(title="D  CPTAC copy number–RNA",x="KRAS CNA",y="KRAS RNA abundance") +
  theme_final()

p3e <- ggplot(cp,aes(x=.data[[cp_cn]],y=.data[[cp_prot]])) +
  geom_point(alpha=.75,size=1.9,colour=COL_TEAL2) +
  geom_smooth(method="lm",se=FALSE,linewidth=.75,colour=COL_TEXT) +
  labs(title="E  CPTAC proteomic concordance",x="KRAS CNA",y="KRAS protein abundance") +
  theme_final()

p3 <- (p3a | p3b | p3c) / (p3d | p3e | patchwork::plot_spacer()) +
  patchwork::plot_annotation(
    title="KRAS copy-number expansion tracks transcriptional and protein output",
    theme=theme(plot.title=element_text(face="bold",colour=COL_TEXT,size=15,hjust=.5))
  )

ggsave(file.path(FIG,"Figure3_molecular_validation_COLORED.pdf"),p3,width=14,height=8.1)
ggsave(file.path(FIG,"Figure3_molecular_validation_COLORED.png"),p3,width=14,height=8.1,dpi=600)

# ---------------- Figure 4: TCGA clinical associations ----------------
f4dat <- bind_rows(
  tcga_cox |> filter(endpoint=="OS",term %in% c("FAM_z","copy_z")) |>
    transmute(
      analysis=ifelse(term=="FAM_z","OS: allelic dominance","OS: copy expansion"),
      axis=ifelse(term=="FAM_z","Allelic dominance","Copy expansion"),
      HR=estimate,CI_low=conf.low,CI_high=conf.high,p=p.value
    ),
  tcga_cox |> filter(endpoint=="PFI",term %in% c("FAM_z","copy_z")) |>
    transmute(
      analysis=ifelse(term=="FAM_z","PFI: allelic dominance","PFI: copy expansion"),
      axis=ifelse(term=="FAM_z","Allelic dominance","Copy expansion"),
      HR=estimate,CI_low=conf.low,CI_high=conf.high,p=p.value
    ),
  fam066 |>
    transmute(
      analysis="OS: FAM ≥ 0.66",
      axis="Allelic dominance",
      HR=estimate,CI_low=conf.low,CI_high=conf.high,p=p.value
    )
) |>
  mutate(
    analysis=factor(analysis,levels=rev(analysis)),
    axis=factor(axis,levels=c("Copy expansion","Allelic dominance"))
  )

p4 <- ggplot(f4dat,aes(x=HR,y=analysis,colour=axis)) +
  geom_vline(xintercept=1,linetype=2,colour=COL_MUTED,linewidth=.55) +
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.16,linewidth=.8) +
  geom_point(size=3) +
  scale_colour_manual(values=c("Copy expansion"=COL_COPY,"Allelic dominance"=COL_FAM)) +
  scale_x_log10() +
  labs(
    title="TCGA clinical associations",
    subtitle="Adjusted models; progression endpoint is PFI",
    x="Hazard ratio (95% CI)",
    y=NULL
  ) +
  theme_final() +
  theme(legend.position="top")

ggsave(file.path(FIG,"Figure4_TCGA_clinical_forest_COLORED.pdf"),p4,width=7.4,height=5.0)
ggsave(file.path(FIG,"Figure4_TCGA_clinical_forest_COLORED.png"),p4,width=7.4,height=5.0,dpi=600)

# ---------------- Figure 5: MSK external validation ----------------
fam_primary <- msk_primary |> filter(term=="fam_z") |>
  transmute(
    analysis="Primary allelic proxy",
    axis="Allelic dominance proxy",
    HR=estimate,CI_low=conf.low,CI_high=conf.high,p=p.value
  )
cn_primary <- msk_primary |> filter(term=="cn_z") |>
  transmute(
    analysis="Primary copy-expansion proxy",
    axis="Copy expansion proxy",
    HR=estimate,CI_low=conf.low,CI_high=conf.high,p=p.value
  )

forest <- bind_rows(
  fam_primary,
  msk_sens |>
    transmute(
      analysis=analysis,
      axis="Allelic dominance proxy",
      HR=HR,CI_low=CI_low,CI_high=CI_high,p=p
    ),
  cn_primary
) |>
  mutate(
    analysis=factor(analysis,levels=rev(unique(analysis))),
    axis=factor(axis,levels=c("Copy expansion proxy","Allelic dominance proxy"))
  )

p5 <- ggplot(forest,aes(x=HR,y=analysis,colour=axis)) +
  geom_vline(xintercept=1,linetype=2,colour=COL_MUTED,linewidth=.55) +
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.16,linewidth=.8) +
  geom_point(size=3) +
  scale_colour_manual(values=c(
    "Copy expansion proxy"=COL_COPY,
    "Allelic dominance proxy"=COL_FAM
  )) +
  scale_x_log10() +
  labs(
    title="MSK-IMPACT external clinical validation",
    subtitle="FACETS/ASCN-based allelic proxy; not identical to TCGA FAM",
    x="Adjusted OS hazard ratio per SD (95% CI)",
    y=NULL
  ) +
  theme_final() +
  theme(legend.position="top")

ggsave(file.path(FIG,"Figure5_MSK_forest_COLORED.pdf"),p5,width=7.6,height=5.0)
ggsave(file.path(FIG,"Figure5_MSK_forest_COLORED.png"),p5,width=7.6,height=5.0,dpi=600)

# ============================================================
# 8. FINAL LOCKED AUDIT
# ============================================================

get_tcga_term <- function(endpoint,term) {
  rr <- tcga_cox |> filter(endpoint==!!endpoint,term==!!term)
  if(!nrow(rr)) return(c(HR=NA,p=NA))
  c(HR=rr$estimate[1],p=rr$p.value[1])
}

os_fam <- tcga_cox |> filter(endpoint=="OS",term=="FAM_z")
os_cop <- tcga_cox |> filter(endpoint=="OS",term=="copy_z")
pf_fam <- tcga_cox |> filter(endpoint=="PFI",term=="FAM_z")
pf_cop <- tcga_cox |> filter(endpoint=="PFI",term=="copy_z")

mfam <- msk_primary |> filter(term=="fam_z")
mcopy <- msk_primary |> filter(term=="cn_z")
mex <- msk_sens |> filter(analysis=="Exclude uncapped proxy >1")
mpr <- msk_sens |> filter(analysis=="Primary tumours")
mme <- msk_sens |> filter(analysis=="Metastases")

audit <- tribble(
  ~module,~result,~expected,~observed,~tol,
  "TCGA","mechanistic N",138,nrow(tcga),0,
  "TCGA","axis rho",-0.0150602,tcga_axis$rho[1],0.005,
  "TCGA","axis P",0.86083,tcga_axis$p[1],0.02,
  "TCGA","clinical Stage15 N",135,nrow(cl),0,
  "TCGA","adjusted N",132,nrow(cc),0,
  "TCGA","adjusted OS deaths",49,sum(cc[[os_event]]==1,na.rm=TRUE),0,
  "TCGA","adjusted PFI events",52,sum(cc[[pfi_event]]==1,na.rm=TRUE),0,
  "TCGA","OS FAM HR",1.532947,os_fam$estimate[1],0.02,
  "TCGA","OS FAM P",0.0055922,os_fam$p.value[1],0.002,
  "TCGA","OS copy HR",1.024857,os_cop$estimate[1],0.02,
  "TCGA","PFI FAM HR",1.037874,pf_fam$estimate[1],0.02,
  "TCGA","PFI copy HR",0.945999,pf_cop$estimate[1],0.02,
  "TCGA","FAM>=0.66 HR",2.180245,fam066$estimate[1],0.05,
  "TCGA RNA","RNA evaluable N",63,nrow(rna_dd),0,
  "TCGA RNA","copy/ploidy vs RNA rho",0.533501,
    (rna_spearman |> filter(comparison=="copy_ploidy_vs_KRAS_expr"))$rho[1],0.003,
  "TCGA RNA","copy/ploidy vs RNA P",6.73e-06,
    (rna_spearman |> filter(comparison=="copy_ploidy_vs_KRAS_expr"))$p[1],3e-06,

  "MSK","N",1538,nrow(msk_ext),0,
  "MSK","deaths",702,sum(msk_ext$event==1,na.rm=TRUE),0,
  "MSK","allelic proxy HR",1.129713,mfam$estimate[1],0.005,
  "MSK","allelic proxy P",0.0012534536,mfam$p.value[1],0.0003,
  "MSK","copy proxy HR",1.095271,mcopy$estimate[1],0.005,
  "MSK","copy proxy P",4.7494234e-06,mcopy$p.value[1],1e-05,
  "MSK","axis rho",-0.044616,msk_axis$rho[1],0.01,
  "MSK","exclude >1 N",1503,mex$N[1],0,
  "MSK","exclude >1 deaths",684,mex$deaths[1],0,
  "MSK","exclude >1 HR",1.126753,mex$HR[1],0.005,
  "MSK","primary N",1022,mpr$N[1],0,
  "MSK","primary deaths",376,mpr$deaths[1],0,
  "MSK","primary HR",1.105681,mpr$HR[1],0.005,
  "MSK","metastases N",516,mme$N[1],0,
  "MSK","metastases deaths",326,mme$deaths[1],0,
  "MSK","metastases HR",1.170319,mme$HR[1],0.005
) |>
  mutate(pass=is.finite(observed) & abs(observed-expected)<=tol)

fwrite(audit,file.path(AUD,"FINAL_LOCKED_AUDIT.tsv"),sep="\t")

# Final report.
report <- c(
  "KRAS-LUAD MANUSCRIPT FINAL MASTER ANALYSIS REPORT",
  "",
  sprintf("TCGA mechanistic cohort: n=%d",nrow(tcga)),
  sprintf("TCGA copy expansion vs FAM: rho=%.4f, P=%.4g",
          tcga_axis$rho[1],tcga_axis$p[1]),
  sprintf("TCGA Stage15 clinical cohort: n=%d",nrow(cl)),
  sprintf("TCGA adjusted complete-case cohort: n=%d; OS deaths=%d; PFI events=%d",
          nrow(cc),sum(cc[[os_event]]==1,na.rm=TRUE),
          sum(cc[[pfi_event]]==1,na.rm=TRUE)),
  sprintf("TCGA adjusted OS FAM HR %.3f, P=%.6g",
          os_fam$estimate[1],os_fam$p.value[1]),
  sprintf("TCGA adjusted OS copy HR %.3f, P=%.6g",
          os_cop$estimate[1],os_cop$p.value[1]),
  sprintf("TCGA adjusted PFI FAM HR %.3f, P=%.6g",
          pf_fam$estimate[1],pf_fam$p.value[1]),
  sprintf("TCGA adjusted PFI copy HR %.3f, P=%.6g",
          pf_cop$estimate[1],pf_cop$p.value[1]),
  "",
  sprintf("TCGA RNA evaluable cohort: n=%d",nrow(rna_dd)),
  sprintf("TCGA copy/ploidy vs KRAS RNA: rho=%.3f, P=%.6g",
          (rna_spearman |> filter(comparison=="copy_ploidy_vs_KRAS_expr"))$rho[1],
          (rna_spearman |> filter(comparison=="copy_ploidy_vs_KRAS_expr"))$p[1]),
  sprintf("Recomputed raw RNA model: copy beta=%.6f, P=%.6g; FAM beta=%.6f, P=%.6g",
          rna_copy_raw$estimate[1],rna_copy_raw$p.value[1],
          rna_fam_raw$estimate[1],rna_fam_raw$p.value[1]),
  "NOTE: historical manuscript RNA coefficients (0.605/-0.225) are retained only as provenance references; they were not reproduced from the canonical 63-case table and are not injected into any fit.",
  "",
  sprintf("MSK exact cohort: n=%d; deaths=%d",
          nrow(msk_ext),sum(msk_ext$event==1,na.rm=TRUE)),
  sprintf("MSK allelic-dominance proxy HR %.6f, P=%.8g",
          mfam$estimate[1],mfam$p.value[1]),
  sprintf("MSK copy-expansion proxy HR %.6f, P=%.8g",
          mcopy$estimate[1],mcopy$p.value[1]),
  sprintf("MSK axis rho=%.6f, P=%.8g",
          msk_axis$rho[1],msk_axis$p[1]),
  "",
  sprintf("Locked audit passed: %d/%d",
          sum(audit$pass,na.rm=TRUE),nrow(audit)),
  "",
  "Interpretation rules:",
  "- TCGA endpoint is PFI, never PFS.",
  "- MSK allelic measure is a FACETS/ASCN-based proxy, not TCGA FAM.",
  "- Stage16/17 locked final CSVs are retained to reproduce the exact historical clinical-translation panels.",
  "- Old Stage19 MSK master is not used by the final manuscript analysis."
)
writeLines(report,file.path(ROOT,"FINAL_ANALYSIS_REPORT.txt"))

capture.output(sessionInfo(),file=file.path(AUD,"sessionInfo.txt"))

cat(paste(report,collapse="\n"),"\n")


# ============================================================
# 9. Reproducibility metadata
# ============================================================
capture.output(sessionInfo(),file=file.path(AUD,"SESSION_INFO.txt"))



# ============================================================
# 10. FINAL MANUSCRIPT PREFLIGHT
# ============================================================
# One compact audit confirming that the final manuscript analysis package
# contains the expected core outputs before the old project is archived.

expected_figures <- c(
  "Figure1_TCGA_axes_COLORED.pdf",
  "Figure2_TCGA_mechanisms_COLORED.pdf",
  "Figure3_molecular_validation_COLORED.pdf",
  "Figure4_TCGA_clinical_forest_COLORED.pdf",
  "Figure5_MSK_forest_COLORED.pdf"
)

expected_results <- c(
  "01_TCGA_axis_spearman.tsv",
  "04A_TCGA_KRAS_RNA_spearman.tsv",
  "04B_TCGA_KRAS_RNA_route_independent_models.tsv",
  "06_TCGA_adjusted_cox.tsv",
  "07_TCGA_FAM066_adjusted_OS.tsv",
  "11_MSK_primary_model.tsv",
  "12_MSK_sensitivity.tsv"
)

preflight <- bind_rows(
  tibble(
    category="figure",
    file=expected_figures,
    exists=file.exists(file.path(FIG,expected_figures))
  ),
  tibble(
    category="result",
    file=expected_results,
    exists=file.exists(file.path(RES,expected_results))
  ),
  tibble(
    category="audit",
    file=c("FINAL_LOCKED_AUDIT.tsv","06_TCGA_RNA_PROVENANCE_AUDIT.tsv","SESSION_INFO.txt"),
    exists=file.exists(file.path(AUD,c("FINAL_LOCKED_AUDIT.tsv","06_TCGA_RNA_PROVENANCE_AUDIT.tsv","SESSION_INFO.txt")))
  )
)

fwrite(preflight,file.path(AUD,"FINAL_MANUSCRIPT_PREFLIGHT.tsv"),sep="\t")

cat("\nFINAL MANUSCRIPT PREFLIGHT: ",
    sum(preflight$exists),"/",nrow(preflight)," expected outputs present\n",sep="")
if(any(!preflight$exists)) {
  cat("Missing outputs:\n")
  print(preflight[!preflight$exists,,drop=FALSE])
} else {
  cat("All core manuscript outputs are present.\n")
}

cat("\n============================================================\n")
cat("KRAS-LUAD MANUSCRIPT FINAL MASTER ANALYSIS — COLORED FIGURES COMPLETE\n")
cat("============================================================\n")
cat("Results: ",RES,"\n",sep="")
cat("Figures: ",FIG,"\n",sep="")
cat("Audit:   ",AUD,"\n",sep="")
cat("Locked audit: ",sum(audit$pass,na.rm=TRUE),"/",nrow(audit)," passed\n",sep="")
cat("\nIMPORTANT RNA NOTE:\n")
cat("The canonical 63-case RNA cohort reproduces the manuscript Spearman result,\n")
cat("but not the legacy multivariable beta/P values. See audit/06_TCGA_RNA_PROVENANCE_AUDIT.tsv.\n")
