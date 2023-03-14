library(dplyr); library(tableone)
library(ggplot2); library(htmlTable)
library(ipw); library(stargazer)
library(survival); library(survminer)
library(finalfit); library(table1)
library(geepack); library(forcats)
library(forestplot); library(DiagrammeR)

library(lme4)
library(gtsummary)
library(patchwork)

source('./R/poisson_ipwtm.R')

dt_forest_plot <- function(OM_fit, BCM_fit, label) {
  datd.uniquerx <- data.frame(cbind(
    HR.A=round(exp(OM_fit$coefficients),2), 
    round(exp(confint(OM_fit)),2),
    pval.A=ifelse(round(summary(OM_fit)$coefficients[,6],3)<0.001,
                  summary(OM_fit)$coefficients[,6],
                  round(summary(OM_fit)$coefficients[,6],2)),
    HR.B=round(exp(BCM_fit$coefficients),2), 
    round(exp(confint(BCM_fit)),2),
    pval.B=ifelse(round(summary(BCM_fit)$coefficients[,6],3)<0.001,
                  summary(BCM_fit)$coefficients[,6],
                  round(summary(BCM_fit)$coefficients[,6],2)),
    fullHR.A = paste(round(exp(OM_fit$coefficients),2)," (",
                     round(exp(confint(OM_fit))[,1],2),", ",
                     round(exp(confint(OM_fit))[,2],2),")",sep=''),
    fullHR.B = paste(round(exp(BCM_fit$coefficients),2)," (",
                     round(exp(confint(BCM_fit))[,1],2),", ",
                     round(exp(confint(BCM_fit))[,2],2),")",sep=''))) %>%
    mutate(vars=c(
      label,
      'Age at diagnosis (per 1 year)',
      'Hispanic (vs. White)',
      'Asian/Pacific Islander (vs. White)',
      'Black (vs. White)',
      'Socioeconomic status quintile (vs. lowest)',
      'Stage II (vs. I)',
      'Stage III (vs. I)',
      'Tumor grade 2 (vs. 1)',
      'Tumor grade 3 (vs. 1)',
      'Tumor grade unknown',
      'Received chemotherapy (vs. no chemotherapy)',
      'Received radiotherapy (vs. no radiotherapy)',
      'Ever used growth factor support (vs. never)',
      'Bilateral mastectomy (vs. lumpectomy)',
      'Unilateral mastectomy (vs. lumpectomy)',
      'No surgery'),
      HR.A=as.numeric(HR.A),
      HR.B=as.numeric(HR.B),
      X2.5..=as.numeric(X2.5..),
      X97.5..=as.numeric(X97.5..),
      X2.5...1=as.numeric(X2.5...1),
      X97.5...1=as.numeric(X97.5...1)) %>%
    rename(LCL.A=X2.5.., UCL.A=X97.5.., LCL.B=X2.5...1, UCL.B=X97.5...1)
  
  # # Removing variables suggested by George
  # datd.uniquerx <- datd.uniquerx %>% 
  #   filter(stringr::str_detect(vars, 'antibiotic|stage|grade|G-CSF|chemo|radio'))
  
  # # Omit unknown tumor grade and no surgery
  # datd.uniquerx <- datd.uniquerx %>%
  #   filter(!(vars %in% c('Tumor grade unknown', 'No surgery')))
  
  # Omit all covariates
  datd.uniquerx <- datd.uniquerx %>%
    filter(vars == label)
  
  return(datd.uniquerx)
}

d <- readRDS('./data/TNBC_data_baseline.rds')
dlong <- readRDS('./data/TNBC_data_longitudinal.rds')


# TABLE 1 -----------------------------------------------------------------

labels <- list(
  variables=list(
    dx_age_cat='Age at diagnosis',
    ses_cat='Socioeconomic status quintile',
    race_eth='Race and ethnicity',
    bmi_cat='Body mass index',
    PUBLIC_UNINSURED_EVER='Ever publicly insured',
    INSTITUTION='Place of care',
    stage='Stage',
    grade='Tumor grade',
    surg_type='Surgery type',
    CT_1yr='Received chemotherapy',
    RT='Received radiotherapy',
    brca3='Germline BRCA1/2 pathogenic variant status',
    neutropenic_post='Ever low absolute neutrophil count (<1K/uL) post-diagnosis',
    min.ANC='Minimum absolute neutrophil count (K/uL) post-diagnosis',
    lymphopenic_post='Ever low absolute lymphocyte count (<1K/uL) post diagnosis',
    min.ALC='Minimum absolute lymphocyte count (K/uL) post-diagnosis',
    GF_ever='Ever used growth factor support',
    n_GFuse='Number of growth factor uses',
    time2death='Follow-up time (months)'),
  groups=list("Antibiotic usage", ""))

tab1 <- d %>%
  select(ever_anti, labels$variables %>% names()) %>%
  labelled::set_variable_labels(.labels = labels$variables) %>%
  tbl_summary(by = ever_anti,
              type = list(n_GFuse ~ "continuous")) %>% 
  add_overall(last = T) %>% 
  add_p() %>% 
  modify_spanning_header(c('stat_1', 'stat_2') ~ '**Antibiotic usage**')

saveRDS(tab1, './output/tab1.rds')



# ANY ANTIBIOTICS ---------------------------------------------------------

# Marginal Structural Model -----------------------------------------------
dlong1 <- dlong %>%
  filter(!is.na(LYMAB))

temp1 <- ipwtm(exposure = anti_start, family = "survival",
               numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
               denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
               id = ANON_ID, tstart = tstart, timevar = tstop, type = "first", data = dlong1)

ipwtm_alc_ever <- tbl_regression(temp1$den.mod, exponentiate = T)
saveRDS(ipwtm_alc_ever, './output/ipwtm_alc_ever.rds')

# Adjusted Associations with Survival for MSM
OMa_ever <- coxph(Surv(tstart, tstop, death) ~ anti_start + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp1$ipw.weights, robust=T)
BCMa_ever <- coxph(Surv(tstart, tstop, bc_death) ~ anti_start + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp1$ipw.weights, robust=T)

# Adjusting for disease severity
temp1b <- ipwtm(exposure = anti_start, family = "survival",
                numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness,
                id = ANON_ID, tstart = tstart, timevar = tstop, type = "first", data = dlong1)

OMa_ever2 <- coxph(Surv(tstart, tstop, death) ~ anti_start + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness +
                     cluster(ANON_ID), data = dlong1, weights = temp1b$ipw.weights, robust=T)

BCMa_ever2 <- coxph(Surv(tstart, tstop, bc_death) ~ anti_start + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness +
                      cluster(ANON_ID), data = dlong1, weights = temp1b$ipw.weights, robust=T)


# Weighted KMs ------------------------------------------------------------
sfit <- survfit(Surv(tstart, tstop, death) ~ anti_start, data=dlong1, weights = temp1$ipw.weights)
sfit2 <- survfit(Surv(tstart, tstop, bc_death) ~ anti_start, data=dlong1, weights = temp1$ipw.weights)

splots <- list()
splots[[1]] <- ggsurvplot(
  sfit,  legend.title = "Any antimicrobial exposure", legend.labs = c("No", "Yes"), risk.table.fontsize=4,
  censor=F, ggtheme = theme_classic(), break.time.by = 24, # break time axis by 24
  risk.table = TRUE, risk.table.y.text.col = TRUE, ylab='Overall survival', risk.table.pos = "in",
  palette = "Dark2",  conf.int = F, xlab='Months from breast cancer diagnosis', legend = c(0.5, 0.25)) +
  guides(colour = guide_legend(nrow = 1))

splots[[2]] <- ggsurvplot(
  sfit2,  legend.title = "Any antimicrobial exposure", legend.labs = c("No", "Yes"), risk.table.fontsize=4,
  censor=F, ggtheme = theme_classic(), break.time.by = 24, # break time axis by 24
  risk.table = TRUE, risk.table.y.text.col = TRUE, ylab='Breast cancer-specific survival', risk.table.pos = "in",
  palette = "Dark2",  conf.int = F, xlab='Months from breast cancer diagnosis', legend = c(0.5, 0.25)) +
  guides(colour = guide_legend(nrow = 1))

ggkm <- splots


# Save regression tables --------------------------------------------------
tbE_om1 <- tbl_regression(
  OMa_ever, exponentiate = T,
  label = list(anti_start ~ 'Any antimicrobial exposure')) %>% 
  modify_header(n_event = '**N Events**')

saveRDS(tbE_om1, './output/aux_cox_reg.rds')

tbE_bcm1 <- tbl_regression(
  BCMa_ever, exponentiate = T,
  label = list(anti_start ~ 'Any antimicrobial exposure')) %>%
  modify_header(n_event = '**N Events**')

adj_msm_ever <- tbl_merge(
  list(tbE_om1, tbE_bcm1),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label)))

saveRDS(adj_msm_ever, './output/adj_msm_ever.rds')

dt_main_fig_1 <- dt_forest_plot(OMa_ever, BCMa_ever, label = 'Any antimicrobial exposure')
saveRDS(dt_main_fig_1, './output/dt_main_fig_1.rds')


# Unweighted Cox ----------------------------------------------------------
OMu_ever <- update(OMa_ever, weights = rep(1, length(OMa_ever$weights)))
BCMu_ever <- update(BCMa_ever, weights = rep(1, length(BCMa_ever$weights)))

tb_OM_ever_uwt <- tbl_merge(list(
  tbl_regression(OMu_ever, exponentiate = T) %>% 
    modify_header(estimate = '**Unweighted**'),
  tbl_regression(OMa_ever, exponentiate = T) %>% 
    modify_header(estimate = '**MSM**')))

tb_BCM_ever_uwt <- tbl_merge(list(
  tbl_regression(BCMu_ever, exponentiate = T) %>% 
    modify_header(estimate = '**Unweighted**'),
  tbl_regression(BCMa_ever, exponentiate = T) %>% 
    modify_header(estimate = '**MSM**')))

adj_unwt_ever_alc <- tbl_merge(
  list(tb_OM_ever_uwt, tb_BCM_ever_uwt),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>% 
  modify_column_hide(starts_with('p.value')) %>% 
  modify_footnote(everything() ~ NA,abbreviation = TRUE) %>% 
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No', ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label),
    label = ifelse(label == 'anti_start', 'Any antimicrobial exposure', label),
    label = ifelse(label == 'stage', 'Cancer stage', label),
    label = ifelse(label == 'ses', 'Socioeconomic status quintile', label)))

saveRDS(adj_unwt_ever_alc, './output/adj_unwt_ever_alc.rds')


# TOTAL ANTIBIOTICS -------------------------------------------------------
dlong1$tot_rx_grp = factor(case_when(
  dlong1$n_tot_anti==0~'None',
  dlong1$n_tot_anti>0 & dlong1$n_tot_anti<=3 ~'1-3',
  dlong1$n_tot_anti>3 & dlong1$n_tot_anti<=7 ~'4-7',
  dlong1$n_tot_anti>7~'8+'), levels=c('None','1-3','4-7','8+'))

### with fewer groups
dlong1$tot_rx_grp1 = factor(case_when(
  dlong1$n_tot_anti==0~'None',
  dlong1$n_tot_anti>0 & dlong1$n_tot_anti<=5 ~'1-5',
  dlong1$n_tot_anti>5 ~ '6+'), levels=c('None','1-5','6+'))

dlong1$lymph = factor(ifelse(dlong1$LYMAB<1,'Low ALC','Normal ALC'))

temp.ntot <- poisson_ipwtm(exposure = n_tot_anti, family = "poisson", 
                   numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type, 
                   denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                   id = ANON_ID, timevar = tstop, type = "all", data = dlong1, corstr="exchangeable", trunc=0.05)

ipwtm_alc_tot <- tbl_regression(temp.ntot$den.mod)
saveRDS(ipwtm_alc_tot, './output/ipwtm_alc_tot.rds')

#### Adjusted Associations with Survival for MSM
OMa_tot <- coxph(Surv(tstart, tstop, death) ~ n_tot_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp.ntot$weights.trunc, robust=T)

#### Adjusted Associations with Survival for MSM
BCMa_tot <- coxph(Surv(tstart, tstop, bc_death) ~ n_tot_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp.ntot$weights.trunc, robust=T)

#### Adjusting for disease severity ####
temp.ntotb <- poisson_ipwtm(exposure = n_tot_anti, family = "poisson", 
                    numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness, 
                    denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness,
                    id = ANON_ID, timevar = tstop, type = "all", data = dlong1, corstr="exchangeable", trunc=0.01)

OMa_tot2 <- coxph(Surv(tstart, tstop, death) ~ n_tot_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness +
                    cluster(ANON_ID), data = dlong1, weights = temp.ntotb$weights.trunc, robust=T)

#### Adjusted Associations with Survival for MSM
BCMa_tot2 <- coxph(Surv(tstart, tstop, bc_death) ~ n_tot_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness +
                     cluster(ANON_ID), data = dlong1, weights = temp.ntotb$weights.trunc, robust=T)

# Weighted KMs ------------------------------------------------------------
sfit <- survfit(Surv(tstart, tstop, death) ~ tot_rx_grp, data=dlong1, weights = temp.ntot$weights.trunc)
sfit2 <- survfit(Surv(tstart, tstop, bc_death) ~ tot_rx_grp, data=dlong1, weights = temp.ntot$weights.trunc)

splots <- list()
splots[[1]] <- ggsurvplot(
  sfit,  legend.title = "Total antimicrobial exposures",
  legend.labs = c('None','1-3','4-7','8+'), risk.table.fontsize=4,
  censor=F, ggtheme = theme_classic(), break.time.by = 24, limits=c(0,216), # break time axis by 24
  risk.table = TRUE, risk.table.y.text.col = TRUE, ylab='Overall survival', risk.table.pos = "in",
  palette = "Dark2",  conf.int = F, xlab='Months from breast cancer diagnosis', legend = c(0.5, 0.35)) +
  guides(colour = guide_legend(nrow = 1))

splots[[2]] <- ggsurvplot(
  sfit2,  legend.title = "Total antimicrobial exposures",
  legend.labs = c('None','1-3','4-7','8+'), risk.table.fontsize=4,
  censor=F, ggtheme = theme_classic(), break.time.by = 24, limits=c(0,216), # break time axis by 24
  risk.table = TRUE, risk.table.y.text.col = TRUE, ylab='Breast cancer-specific survival', risk.table.pos = "in",
  palette = "Dark2",  conf.int = F, xlab='Months from breast cancer diagnosis', legend = c(0.5, 0.35)) +
  guides(colour = guide_legend(nrow = 1))

ggkm <- c(ggkm, splots)

hist1 <- ggplot(dlong1) +
  geom_histogram(aes(x = n_tot_anti, y = ..density.., weight = temp.ntot$weights.trunc),
                 color = 'black', fill = 'white') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Number of total antimicrobial prescriptions', y = 'Density') +
  theme_classic()

tbE_om1 <- tbl_regression(
  OMa_tot, exponentiate = T,
  label = list(n_tot_anti ~ 'Total antimicrobial exposures')) %>% 
  modify_header(n_event = '**N Events**')

tbE_bcm1 <- tbl_regression(
  BCMa_tot, exponentiate = T,
  label = list(n_tot_anti ~ 'Total antimicrobial exposures')) %>% 
  modify_header(n_event = '**N Events**')

adj_msm_tot <- tbl_merge(
  list(tbE_om1, tbE_bcm1),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label)))

saveRDS(adj_msm_tot, './output/adj_msm_tot.rds')

dt_main_fig_2 <- dt_forest_plot(OMa_tot, BCMa_tot, label = 'Total antimicrobial exposures')
saveRDS(dt_main_fig_2, './output/dt_main_fig_2.rds')


OMu_tot <- update(OMa_tot, weights = rep(1, length(OMa_tot$weights)))
BCMu_tot <- update(BCMa_tot, weights = rep(1, length(BCMa_tot$weights)))

tbE_om1u <- tbl_regression(
  OMu_tot, exponentiate = T,
  label = list(n_tot_anti ~ 'Total antimicrobial exposures')) %>% 
  modify_header(n_event = '**N Events**')

tbE_bcm1u <- tbl_regression(
  BCMu_tot, exponentiate = T,
  label = list(n_tot_anti ~ 'Total antimicrobial exposures')) %>% 
  modify_header(n_event = '**N Events**')

adj_unwt_tot <- tbl_merge(
  list(tbE_om1u, tbE_bcm1u),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label)))

saveRDS(adj_unwt_tot, './output/adj_unwt_tot.rds')


# UNIQUE ANTIBIOTICS -----------------------------------------------------

dlong1$uniq_rx_grp = factor(case_when(
  dlong1$n_uniq_anti==0~'None',
  dlong1$n_uniq_anti>0 & dlong1$n_uniq_anti<=2 ~'1-2',
  dlong1$n_uniq_anti>2 & dlong1$n_uniq_anti<=4 ~'3-4',
  dlong1$n_uniq_anti>4~'5+'), levels=c('None','1-2','3-4','5+'))

temp.nuni <- poisson_ipwtm(exposure = n_uniq_anti, family = "poisson", 
                   numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type, 
                   denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                   id = ANON_ID, timevar = tstop, type = "all", data = dlong1, corstr="exchangeable", trunc=0.01)

ipwtm_alc_uni <- tbl_regression(temp.nuni$den.mod)
saveRDS(ipwtm_alc_uni, './output/ipwtm_alc_uni.rds')

#### Adjusted Associations with Survival for MSM
OMa_uni <- coxph(Surv(tstart, tstop, death) ~ n_uniq_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp.nuni$weights.trunc, robust=T)

BCMa_uni <- coxph(Surv(tstart, tstop, bc_death) ~ n_uniq_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp.nuni$weights.trunc, robust=T)

#### Adjusting for disease severity ####
temp.nuni2 <- poisson_ipwtm(exposure = n_uniq_anti, family = "poisson",
                   numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness, 
                   denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness,
                   id = ANON_ID, timevar = tstop, type = "all", data = dlong1, corstr="exchangeable", trunc=0.01)

OMa_uni2 <- coxph(Surv(tstart, tstop, death) ~ n_uniq_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness +
                   cluster(ANON_ID), data = dlong1, weights = temp.nuni2$weights.trunc, robust=T)

BCMa_uni2 <- coxph(Surv(tstart, tstop, bc_death) ~ n_uniq_anti + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + acute_illness +
                    cluster(ANON_ID), data = dlong1, weights = temp.nuni2$weights.trunc, robust=T)



tbE_om2 <- tbl_regression(
  OMa_ever2, exponentiate = T,
  label = list(anti_start ~ 'Any antimicrobial exposure')) %>% 
  modify_header(n_event = '**N Events**')

tbE_bcm2 <- tbl_regression(
  BCMa_ever2, exponentiate = T,
  label = list(anti_start ~ 'Any antimicrobial exposure')) %>% 
  modify_header(n_event = '**N Events**')

adj_msm_ever_severe <- tbl_merge(
  list(tbE_om2, tbE_bcm2),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label))))

saveRDS(adj_msm_ever_severe, './output/adj_msm_ever_severe.rds')

tbE_om2 <- tbl_regression(
  OMa_tot2, exponentiate = T,
  label = list(n_tot_anti ~ 'Total antimicrobial exposures')) %>% 
  # add_n(location = "level") %>% 
  # add_nevent(location = "level")
  modify_header(n_event = '**N Events**')

tbE_bcm2 <- tbl_regression(
  BCMa_tot2, exponentiate = T,
  label = list(n_tot_anti ~ 'Total antimicrobial exposures')) %>% 
  # add_n(location = "level") %>% 
  # add_nevent(location = "level")
  modify_header(n_event = '**N Events**')

adj_msm_tot_severe <- tbl_merge(
  list(tbE_om2, tbE_bcm2),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label))))

saveRDS(adj_msm_tot_severe, './output/adj_msm_tot_severe.rds')

tbE_om2 <- tbl_regression(
  OMa_uni2, exponentiate = T,
  label = list(n_uniq_anti ~ 'Unique antimicrobial exposures')) %>% 
  # add_n(location = "level") %>% 
  # add_nevent(location = "level")
  modify_header(n_event = '**N Events**')

tbE_bcm2 <- tbl_regression(
  BCMa_uni2, exponentiate = T,
  label = list(n_uniq_anti ~ 'Unique antimicrobial exposures')) %>% 
  # add_n(location = "level") %>% 
  # add_nevent(location = "level")
  modify_header(n_event = '**N Events**')

adj_msm_uni_severe <- tbl_merge(
  list(tbE_om2, tbE_bcm2),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label))))

saveRDS(adj_msm_uni_severe, './output/adj_msm_uni_severe.rds')

sfit <- survfit(Surv(tstart, tstop, death) ~ uniq_rx_grp, data=dlong1, weights = temp.nuni$weights.trunc)
sfit2 <- survfit(Surv(tstart, tstop, bc_death) ~ uniq_rx_grp, data=dlong1, weights = temp.nuni$weights.trunc)

splots <- list()
splots[[1]] <- ggsurvplot(
  sfit, legend.title = "Unique antimicrobial exposures",
  legend.labs = c('None','1-2','3-4','5+'), risk.table.fontsize=4,
  censor=F, ggtheme = theme_classic(), break.time.by = 24, # break time axis by 24
  risk.table = TRUE, risk.table.y.text.col = TRUE, ylab='Overall survival', risk.table.pos = "in",
  palette = "Dark2",  conf.int = F, xlab='Months from breast cancer diagnosis', legend = c(0.5, 0.35)) +
  guides(colour = guide_legend(nrow = 1))
splots[[2]] <- ggsurvplot(
  sfit2,  legend.title = "Unique antimicrobial exposures",
  legend.labs = c('None','1-2','3-4','5+'), risk.table.fontsize=4,
  censor=F, ggtheme = theme_classic(), break.time.by = 24, # break time axis by 24
  risk.table = TRUE, risk.table.y.text.col = TRUE, ylab='Breast cancer-specific survival', risk.table.pos = "in",
  palette = "Dark2",  conf.int = F, xlab='Months from breast cancer diagnosis', legend = c(0.5, 0.35)) +
  guides(colour = guide_legend(nrow = 1))

ggkm <- c(ggkm, splots)
saveRDS(ggkm, './output/ggsurvplots.rds')

hist2 <- ggplot(dlong1) +
  geom_histogram(aes(x = n_uniq_anti, y = ..density.., weight = temp.nuni$ipw.weights),
                 color = 'black', fill = 'white') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Number of unique antimicrobial prescriptions', y = 'Density') +
  theme_classic()

hist <- hist1/hist2 # + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')', tag_sep = '')
ggsave('./output/gghist.png', hist)

tbE_om1 <- tbl_regression(
  OMa_uni, exponentiate = T,
  label = list(n_uniq_anti ~ 'Unique antimicrobial exposures')) %>% 
  modify_header(n_event = '**N Events**')

tbE_bcm1 <- tbl_regression(
  BCMa_uni, exponentiate = T,
  label = list(n_uniq_anti ~ 'Unique antimicrobial exposures')) %>% 
  modify_header(n_event = '**N Events**')

adj_msm_uni <- tbl_merge(
  list(tbE_om1, tbE_bcm1),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label)))

saveRDS(adj_msm_uni, './output/adj_msm_uni.rds')

dt_main_fig_3 <- dt_forest_plot(OMa_uni, BCMa_uni, label = 'Unique antimicrobial exposures')
saveRDS(dt_main_fig_3, './output/dt_main_fig_3.rds')

OMu_uni <- update(OMa_uni, weights = rep(1, length(OMa_uni$weights)))
BCMu_uni <- update(BCMa_uni, weights = rep(1, length(BCMa_uni$weights)))

# tables suggested by Esther
tbE_om1u <- tbl_regression(
  OMu_uni, exponentiate = T,
  label = list(n_uniq_anti ~ 'Unique antimicrobial exposures')) %>% 
  # add_n(location = "level") %>% 
  # add_nevent(location = "level")
  modify_header(n_event = '**N Events**')

tbE_bcm1u <- tbl_regression(
  BCMu_uni, exponentiate = T,
  label = list(n_uniq_anti ~ 'Unique antimicrobial exposures')) %>% 
  # add_n(location = "level") %>% 
  # add_nevent(location = "level")
  modify_header(n_event = '**N Events**')

adj_unwt_uni <- tbl_merge(
  list(tbE_om1u, tbE_bcm1u),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>%
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No',
                   ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label)))

saveRDS(adj_unwt_uni, './output/adj_unwt_uni.rds')

# ANC exposure
## Any Antibiotic Use

dlong1 <- dlong %>%
  filter(!is.na(NEUTAB))

temp1 <- ipwtm(exposure = anti_start, family = "survival", numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type, denominator = ~ NEUTAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
               id = ANON_ID, tstart = tstart, timevar = tstop, type = "first", data = dlong1)

ipwtm_anc_ever <- tbl_regression(temp1$den.mod, exponentiate = T)
saveRDS(ipwtm_anc_ever, './output/ipwtm_anc_ever.rds')

OMa_ever1 <- coxph(Surv(tstart, tstop, death) ~ anti_start + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp1$ipw.weights, robust=T)
BCMa_ever1 <- coxph(Surv(tstart, tstop, bc_death) ~ anti_start + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type + cluster(ANON_ID), data = dlong1, weights = temp1$ipw.weights, robust=T)

OMu_ever1 <- update(OMa_ever1, weights = rep(1, length(OMa_ever1$weights)))
BCMu_ever1 <- update(BCMa_ever1, weights = rep(1, length(BCMa_ever1$weights)))

tb_OM_ever1_uwt <- tbl_merge(list(
  tbl_regression(OMu_ever1, exponentiate = T) %>% 
    modify_header(estimate = '**Unweighted**'),
  tbl_regression(OMa_ever1, exponentiate = T) %>% 
    modify_header(estimate = '**MSM**')))

tb_BCM_ever1_uwt <- tbl_merge(list(
  tbl_regression(BCMu_ever1, exponentiate = T) %>% 
    modify_header(estimate = '**Unweighted**'),
  tbl_regression(BCMa_ever1, exponentiate = T) %>% 
    modify_header(estimate = '**MSM**')))

adj_unwt_ever_anc <- tbl_merge(
  list(tb_OM_ever1_uwt, tb_BCM_ever1_uwt),
  tab_spanner = c('**Overall Survival**', '**Breast Cancer-Specific Survival**')) %>% 
  modify_column_hide(starts_with('p.value')) %>% 
  modify_footnote(everything() ~ NA,abbreviation = TRUE) %>% 
  modify_table_body(~ .x %>% mutate(
    label = ifelse(label == 'N', 'No', ifelse(label == 'Y', 'Yes', label)),
    label = ifelse(label == '9', 'Unknown', label),
    label = ifelse(label == 'anti_start', 'Any antimicrobial exposure', label),
    label = ifelse(label == 'stage', 'Cancer stage', label),
    label = ifelse(label == 'ses', 'Socioeconomic status quintile', label)))

saveRDS(adj_unwt_ever_anc, './output/adj_unwt_ever_anc.rds')



# LANDMARK ANALYSIS -------------------------------------------------------
dlong1 <- dlong %>%
  filter(!is.na(LYMAB)) %>% 
  group_by(ANON_ID) %>%
  mutate(max_tstop = max(tstop)) %>%
  as.data.frame()

out <- data.frame()
for (i in 1:6) {
  temp1 <- ipwtm(exposure = anti_start, family = "survival",
                 numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                 denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                 id = ANON_ID, tstart = tstart, timevar = tstop, type = "first",
                 data = dlong1 %>% filter(max_tstop > (i-1)*12))
  
  om <- coxph(Surv(tstart, tstop, death) ~ anti_start +
                dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type +
                cluster(ANON_ID),
              data = dlong1 %>% filter(max_tstop > (i-1)*12),
              weights = temp1$ipw.weights, robust=T)
  
  out <- bind_rows(out, data.frame(s = 'OS', anti = 'Any', yr = i-1, t(coef(summary(om))[1,])))
  
  bcm <- coxph(Surv(tstart, tstop, bc_death) ~ anti_start +
                 dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type +
                 cluster(ANON_ID),
               data = dlong1 %>% filter(max_tstop > (i-1)*12),
               weights = temp1$ipw.weights, robust=T)
  
  out <- bind_rows(out, data.frame(s = 'BCS', anti = 'Any', yr = i-1, t(coef(summary(bcm))[1,])))
}

for (i in 1:6) {
  temp.ntot <- poisson_ipwtm(exposure = n_tot_anti, family = "poisson", 
                        numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type, 
                        denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                        id = ANON_ID, timevar = tstop, type = "all", corstr="exchangeable", trunc=0.05,
                        data = dlong1 %>% filter(max_tstop > (i-1)*12))
  
  
  om <- coxph(Surv(tstart, tstop, death) ~ n_tot_anti + 
                dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type +
                cluster(ANON_ID),
              data = dlong1 %>% filter(max_tstop > (i-1)*12),
              weights = temp.ntot$weights.trunc, robust=T)
  
  out <- bind_rows(out, data.frame(s = 'OS', anti = 'Total', yr = i-1, t(coef(summary(om))[1,])))
  
  bcm <- coxph(Surv(tstart, tstop, bc_death) ~ n_tot_anti +
                 dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type +
                 cluster(ANON_ID),
               data = dlong1 %>% filter(max_tstop > (i-1)*12),
               weights = temp.ntot$weights.trunc, robust=T)
  
  out <- bind_rows(out, data.frame(s = 'BCS', anti = 'Total', yr = i-1, t(coef(summary(bcm))[1,])))
}

for (i in 1:6) {
  temp.nuni <- poisson_ipwtm(exposure = n_uniq_anti, family = "poisson", 
                        numerator = ~ dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type, 
                        denominator = ~ LYMAB + dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type,
                        id = ANON_ID, timevar = tstop, type = "all", corstr="exchangeable", trunc=0.01,
                        data = dlong1 %>% filter(max_tstop > (i-1)*12))
  
  
  om <- coxph(Surv(tstart, tstop, death) ~ n_uniq_anti +
                dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type +
                cluster(ANON_ID),
              data = dlong1 %>% filter(max_tstop > (i-1)*12),
              weights = temp.nuni$weights.trunc, robust=T)
  
  out <- bind_rows(out, data.frame(s = 'OS', anti = 'Unique', yr = i-1, t(coef(summary(om))[1,])))
  
  bcm <- coxph(Surv(tstart, tstop, bc_death) ~ n_uniq_anti +
                 dx_age + race_eth + ses + stage + grade + CT_1yr + RT + GF_ever + surg_type +
                 cluster(ANON_ID),
               data = dlong1 %>% filter(max_tstop > (i-1)*12),
               weights = temp.nuni$ipw.weights, robust=T)
  
  out <- bind_rows(out, data.frame(s = 'BCS', anti = 'Unique', yr = i-1, t(coef(summary(bcm))[1,])))
}

saveRDS(out, file = './output/landmark.rds')


