suppressMessages({library(readxl); library(dplyr); library(tidyr); library(writexl)})
setwd("/Users/to909/Desktop/Meg/Dose Leyre/Body comp data")
d  <- suppressWarnings(read_excel("Data/OncoGFR1200_MGH_chemo_CJASN.xlsx"))
fr <- suppressWarnings(read_excel("Data/Frailty data to MGH.xlsx"))   # cols: protocal, fried_classification

reg <- data.frame(
  equation = c("CKDEPIcr21","EKFCcr","CKDEPIcys","EKFCcys","CKDEPIcrcys21","EKFCcrcys","CG","CKDEPI4mark"),
  type     = c("eGFRcr","eGFRcr","eGFRcys","eGFRcys","eGFRcr-cys","eGFRcr-cys","ClCr","Panel"),
  label    = c("CKD-EPI 2021","EKFC 2021","CKD-EPI 2012","EKFC 2023","CKD-EPI 2021","EKFC 2023","Cockcroft-Gault","CKD-EPI 4-marker"),
  idxcol   = c("CKDEPIcr21index","EKFCcrindex","CKDEPIcysindex","EKFCcysindex","CKDEPIcrcys21index","EKFCcrcysindex","CGindex","CKDEPI4mark"),
  stringsAsFactors = FALSE)

# --- flags per patient x equation (indexed GFR, vs measured GFR) ---
long <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
  g <- d[[reg$idxcol[i]]]; m <- d$mGFRindex; pred <- 5 * (g + 25)
  data.frame(
    protocol = d$protocol,
    eGFR_type = reg$type[i], eGFR_equation = reg$label[i], equation_id = reg$equation[i],
    egfr_index = round(g, 1), mgfr_index = round(m, 1),
    cisplatin_eligible = as.integer(g >= 40),
    carbo15_over  = as.integer(pred > 5.75 * (m + 25)),
    carbo15_under = as.integer(pred < 4.25 * (m + 25)),
    carbo20_over  = as.integer(pred > 6.00 * (m + 25)),
    carbo20_under = as.integer(pred < 4.00 * (m + 25)),
    stringsAsFactors = FALSE)
}))
# merge Fried by protocol (fried file key is misspelled 'protocal')
long$fried_classification <- fr$fried_classification[match(long$protocol, fr$protocal)]
long$fried_label <- c("Robust","Pre-frail","Frail")[long$fried_classification + 1]

# --- wide: one row per patient, flags for each equation as columns ---
wide <- long %>%
  select(protocol, equation_id, cisplatin_eligible, carbo15_over, carbo15_under, carbo20_over, carbo20_under) %>%
  pivot_wider(names_from = equation_id,
              values_from = c(cisplatin_eligible, carbo15_over, carbo15_under, carbo20_over, carbo20_under),
              names_glue = "{equation_id}_{.value}") %>%
  arrange(match(protocol, d$protocol))
wide$fried_classification <- fr$fried_classification[match(wide$protocol, fr$protocal)]
wide$fried_label <- c("Robust","Pre-frail","Frail")[wide$fried_classification + 1]

codebook <- data.frame(Variable = c(
  "protocol","eGFR_type","eGFR_equation","equation_id","egfr_index","mgfr_index",
  "cisplatin_eligible","carbo15_over","carbo15_under","carbo20_over","carbo20_under",
  "fried_classification","fried_label"),
  Definition = c(
  "Patient study ID; matches OncoGFR1200 and the Frailty file",
  "eGFR equation family (eGFRcr / eGFRcys / eGFRcr-cys / ClCr / Panel)",
  "eGFR equation (readable label)","eGFR equation (short id)",
  "This equation's INDEXED eGFR (mL/min/1.73m2)","Measured GFR, indexed (mL/min/1.73m2)",
  "1 = equation GFR >= 40 -> cisplatin-eligible (not 'avoid'); 0 = avoid",
  "1 = Calvert dose >15% ABOVE the mGFR AUC-5 dose (overdose), else 0",
  "1 = Calvert dose >15% BELOW the mGFR AUC-5 dose (underdose), else 0",
  "1 = Calvert dose >20% ABOVE the mGFR AUC-5 dose (overdose), else 0",
  "1 = Calvert dose >20% BELOW the mGFR AUC-5 dose (underdose), else 0",
  "Fried frailty class from 'Frailty data to MGH.xlsx' (0/1/2); NA if not in that file (592 of 1200 have it)",
  "Assumed label for fried_classification: 0=Robust, 1=Pre-frail, 2=Frail (CONFIRM mapping)"),
  stringsAsFactors = FALSE)

out <- "Data/OncoGFR1200_dose_flags_fried.xlsx"
write_xlsx(list(codebook = codebook, flags_by_equation = long, wide_by_patient = wide), out)
cat("WROTE:", out, "\n")
cat("sheets: codebook, flags_by_equation (", nrow(long), "rows =", nrow(d), "pts x", nrow(reg),
    "eqs), wide_by_patient (", nrow(wide), "rows )\n")
cat("fried merged (non-NA rows, long):", sum(!is.na(long$fried_classification)),
    " ; patients with fried (wide):", sum(!is.na(wide$fried_classification)), "\n")
cat("wide columns:", ncol(wide), "\n")
