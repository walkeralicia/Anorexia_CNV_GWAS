

#============= This script meta-analyses the UKB and ANGI pleiotropic CNV burden results =============================

# read in R libraries =================
library(logistf)
library(readxl)
library(stringr)
library(dplyr)
options(scipen=999)

## Read in UKB and ANGI results ================
wkdir="/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs"
ukb <- read.table(paste(wkdir, "UKBB_drCNVs_burden.txt", sep="/"), header=T)
angi <- read.table("/QRISdata/Q4399/Anorexia/ANGI/results/ANGI_drCNVs_burden.txt", header=T)

## split data by method and reformat ==========
methods <- c("Logistic", "Unadjusted_Logistic", "Firth", "Unadjusted_Firth")

results <- lapply(methods, function(method) {
  ukb_method <- ukb %>% filter(Method == method)
  angi_method <- angi %>% filter(Method == method)
  comb <- merge(angi_method, ukb_method, by = "Index", all = TRUE)
  cols_to_replace_na <- c("Case_N.y", "Cont_N.y", "Case_N.x", "Cont_N.x", "contaunzusa2yes", "contseyes", "caseaunzusa2yes", "caseseyes", "Cont_Freq_per10K.x", "Cont_Freq_per10K.y")
  comb[cols_to_replace_na] <- lapply(comb[cols_to_replace_na], function(x) ifelse(is.na(x), 0, x))
  return(as.data.frame(comb))
})

# .x is ANGI and .y is the UKB.
#caseaunzusa2yes = Case + aunzusa2 + CNV present
#contaunzusa2yes = Control + aunzusa2 + CNV present
#caseseyes = Case + SWE + CNV present
#contseyes = Control + SWE + CNV present 



#================ Meta-analyse the Firth-adjusted regression results with Stouffer's method ====================================

# stouffer meta-analyses ===============================

stouffer <- results[[3]]
stouffer$Zvalue.x <- (-1)*sign(log(stouffer$OR.x))*qnorm(stouffer$Pvalue.x/2)
stouffer$Zvalue.y <- (-1)*sign(log(stouffer$OR.y))*qnorm(stouffer$Pvalue.y/2)
stouffer$Stouffer_Zvalue <- (stouffer$Zvalue.x + stouffer$Zvalue.y)/sqrt(2)
stouffer$Stouffer_Pvalue <- 2*pnorm(-(abs(stouffer$Stouffer_Zvalue)))
stouffer$Stouffer_onesidedPvalue <- ifelse(stouffer$Stouffer_Zvalue>0, stouffer$Stouffer_Pvalue/2, 1-(stouffer$Stouffer_Pvalue/2))
stouffer$beta_meta <- (stouffer$beta.x/(stouffer$se.x*stouffer$se.x) + stouffer$beta.y/(stouffer$se.y*stouffer$se.y)) / (1/(stouffer$se.x*stouffer$se.x)+1/(stouffer$se.y*stouffer$se.y))
stouffer$se_meta <- 1/sqrt(1/(stouffer$se.x*stouffer$se.x) +1/(stouffer$se.y*stouffer$se.y))
stouffer$OR_meta <- signif(exp(stouffer$beta_meta),3)
stouffer$dp <- nchar(stouffer$OR_meta)-2
stouffer$LCI_meta <- round(exp(stouffer$beta_meta-1.96*stouffer$se_meta), stouffer$dp)
stouffer$HCI_meta <- round(exp(stouffer$beta_meta+1.96*stouffer$se_meta), stouffer$dp)
stouffer$CI_meta <- paste0("(", stouffer$LCI_meta, "-", stouffer$HCI_meta, ")")
stouffer <-stouffer[order(stouffer$Stouffer_onesidedPvalue),]


# annotate the information about the pleiotropic CNVs to the result file ========================================================================

collins <- read_excel(paste(wkdir, "1-s2.0-S0092867422007887-mmc4.xlsx", sep="/"), sheet=1)
collins <- as.data.frame(collins)
collins$Location <- paste0(collins$Chrom, ":", collins$Start, ":", collins$End)
collins$Genotype <- ifelse(collins$'CNV Type' == "DEL", "Deletion", "Duplication")


# Replace NA values with 0 ======================
comb <- merge(collins, stouffer, by.x = "Segment ID", by.y="Index", all.x=TRUE)
cols_to_replace_na <- c("Case_N.y", "Cont_N.y", "Case_N.x", "Cont_N.x", "contaunzusa2yes", "contseyes", "caseaunzusa2yes", "caseseyes", "Cont_Freq_per10K.x", "Cont_Freq_per10K.y")
comb[cols_to_replace_na] <- lapply(comb[cols_to_replace_na], function(x) ifelse(is.na(x), 0, x))

# Select and rename columns ========================
res <- comb %>%
  select(Cytoband, Location, Genotype, 
         ANGI_Case_N := Case_N.x, ANGI_Cont_N := Cont_N.x, ANGI_Cont_SE := Cont_SE.x, 
         ANGI_Cont_Freq_10K := Cont_Freq_per10K.x, ANGI_Cont_CI := Cont_CI.x, ANGI_OR := OR.x, 
         ANGI_CI := CI.x, ANGI_Pvalue := Pvalue.x, ANGI_Zvalue := Zvalue.x, 
         UKB_Case_N := Case_N.y, UKB_Cont_N := Cont_N.y, UKB_Cont_SE := Cont_SE.y, 
         UKB_Cont_Freq_10K := Cont_Freq_per10K.y, UKB_Cont_CI := Cont_CI.y, 
         UKB_OR := OR.y, UKB_CI := CI.y, UKB_Pvalue := Pvalue.y, 
         UKB_Zvalue := Zvalue.y, Stouffer_Zvalue, Stouffer_Pvalue, 
         Stouffer_onesidedPvalue, Meta_OR := OR_meta, Meta_CI := CI_meta)


# Fisher's exact test to compare CNV control frequency between ANGI and the UKB ===========
RES <- res 
ukb_controls = 385930
angi_controls = 5044
RES$ANGI_Cont_notN <- angi_controls-as.numeric(RES$ANGI_Cont_N)
RES$UKB_Cont_notN <- ukb_controls-as.numeric(RES$UKB_Cont_N)
n=length(RES$Cytoband)
for (i in 1:n){
  FM=matrix(c(as.numeric(RES$ANGI_Cont_notN[i]),as.numeric(RES$ANGI_Cont_N[i]),
              as.numeric(RES$UKB_Cont_notN[i]),as.numeric(RES$UKB_Cont_N[i])),nrow=2,ncol=2)
  FT=fisher.test(FM)
  RES$Cont_freq_Pvalue[i]=signif(FT$p.value,2)
}


# Reformating and saving final results ================================================

final <- RES %>%
  arrange(Stouffer_onesidedPvalue) %>%
  select(-ANGI_Cont_notN, -UKB_Cont_notN) %>%
  mutate(across(c(Stouffer_onesidedPvalue, ANGI_Zvalue, UKB_Zvalue, ANGI_Pvalue, UKB_Pvalue, Stouffer_Zvalue, Stouffer_Pvalue, Cont_freq_Pvalue), ~signif(., 3)))


write.csv(final, paste(wkdir, "collins_all.csv", sep="/"),  row.names=F)

paper <- final %>% filter(Stouffer_onesidedPvalue<0.05)
write.csv(paper, paste(wkdir, "collins_sig.csv",  sep="/"), row.names=F)




############# EXTRA CODE/ANALYSES ##############################################################################

#========= Firth Logistic Regression with dummy data for ANGI and the UKB (with ANGI study region covariate) ==========

# Counts in each group from ANGI
angi_cont_aunzusa2 = 1373
angi_cont_se = 3654
angi_case_aunzusa2 = 3725
angi_case_se = 3660
# Counts in each group from the UKB
ukb_case = 1246
ukb_cont = 385930

comb <- results[[4]]
res <- list()
for (i in 1:dim(comb)[1]){
  print(i)
  CNV <- c(rep(1,comb$Case_N.y[i]),rep(0,(ukb_case-(comb$Case_N.y[i]))),
           rep(1,comb$Cont_N.y[i]),rep(0,(ukb_cont-(comb$Cont_N.y[i]))), 
           rep(1,comb$caseaunzusa2yes[i]),rep(0,angi_case_aunzusa2-comb$caseaunzusa2yes[i]), 
           rep(1,comb$contaunzusa2yes[i]), rep(0, angi_cont_aunzusa2-comb$contaunzusa2yes[i]),
           rep(1,comb$caseseyes[i]),rep(0,angi_case_se-comb$caseseyes[i]), 
           rep(1,comb$contseyes[i]), rep(0, angi_cont_se-comb$contseyes[i]))
  AN <- c(rep(1, ukb_case), 
          rep(0, ukb_cont),
          rep(1,comb$caseaunzusa2yes[i]),rep(1,angi_case_aunzusa2-comb$caseaunzusa2yes[i]), 
          rep(0, comb$contaunzusa2yes[i]), rep(0, angi_cont_aunzusa2-comb$contaunzusa2yes[i]),
          rep(1,comb$caseseyes[i]),rep(1,angi_case_se-comb$caseseyes[i]), 
          rep(0, comb$contseyes[i]), rep(0, angi_cont_se-comb$contseyes[i]))
  STUDY <- c(rep(2, 382718),
             rep(1,comb$caseaunzusa2yes[i]),rep(1,angi_case_aunzusa2-comb$caseaunzusa2yes[i]), 
             rep(1, comb$contaunzusa2yes[i]), rep(1,angi_cont_aunzusa2-comb$contaunzusa2yes[i]),
             rep(0,comb$caseseyes[i]),rep(0,angi_case_se-comb$caseseyes[i]), 
             rep(0, comb$contseyes[i]), rep(0, angi_cont_se-comb$contseyes[i]))
  if (length(table(CNV))==1){
    vec <- c(NA, comb[i,1], comb[i,2], comb$Case_N.x[i]+comb$Case_N.y[i], comb$Cont_N.x[i]+comb$Cont_N.y[i], NA, NA, NA)
    vec <- unname(vec)
    res[[i]] <- vec
  } else {
    dum <- data.frame(AN, STUDY, CNV)
    #  Firth Logistic Regression
    mfirth = logistf(AN ~ STUDY + CNV, data = dum)
    firth = flac(mfirth, data=dum)
    OR <- signif(exp(coef(firth)[3]),3)
    dp2 <- nchar(OR)-2
    CI <- paste0("(", round(exp(confint(firth)[3,1]), dp2), "-", round(exp(confint(firth)[3,2]), dp2), ")")
    pval <- signif(firth$prob[3],dp2)
    vec <- c("Firth", comb[i,1], comb[i,2], comb$Case_N.x[i]+comb$Case_N.y[i], comb$Cont_N.x[i]+comb$Cont_N.y[i], OR, CI, pval)
    vec <- unname(vec)
    res[[i]] <- vec
  }
}

res_all <- do.call(rbind, res)
res_all <- as.data.frame(res_all)
colnames(res_all) <- c("Method", "Index", "Syndrome", "Case_N", "Cont_N", "OR", "CI", "Pvalue")
write.table(res_all, paste(wkdir, "meta_analyses/Firth_Dummy_RegionCovariate.txt", sep="/"), sep ="\t", quote=T, row.names=F, col.names=T)


#======= Comparing Firth logistic regression to normal logistic regression (without covariates) using dummy data ========================================

comb <- results[[4]]
res2 <- data.frame()
for (i in 1:dim(comb)[1]){
  print(i)
  CNV <- c(rep(1,comb$Case_N[i]),rep(0,(ukb_case-(comb$Case_N[i]))),
           rep(1,comb$Cont_N[i]),rep(0,(ukb_cont-(comb$Cont_N[i]))), 
           rep(1,comb$caseaunzusa2yes[i]),rep(0,angi_case_aunzusa2-comb$caseaunzusa2yes[i]), 
           rep(1,comb$contaunzusa2yes[i]), rep(0, angi_cont_aunzusa2-comb$contaunzusa2yes[i]),
           rep(1,comb$caseseyes[i]),rep(0,angi_case_se-comb$caseseyes[i]), 
           rep(1,comb$contseyes[i]), rep(0, angi_cont_se-comb$contseyes[i]))
  AN <- c(rep(1, ukb_case), 
          rep(0, ukb_cont),
          rep(1,comb$caseaunzusa2yes[i]),rep(1,angi_case_aunzusa2-comb$caseaunzusa2yes[i]), 
          rep(0, comb$contaunzusa2yes[i]), rep(0, angi_cont_aunzusa2-comb$contaunzusa2yes[i]),
          rep(1,comb$caseseyes[i]),rep(1,angi_case_se-comb$caseseyes[i]), 
          rep(0, comb$contseyes[i]), rep(0, angi_cont_se-comb$contseyes[i]))
  if (length(table(CNV))==1) next 
  dum <- data.frame(AN, CNV)
  firth = logistf(AN ~ CNV, data = dum)
  m <- glm(AN~CNV, data = dum, family = "binomial")
  counts <- table(dum$AN,dum$CNV)
  caseyes=counts[2,2]
  caseno=counts[2,1]
  contyes=counts[1,2]
  contno=counts[1,1]
  ###NORMAL
  pyes=caseyes/(caseyes+contyes)
  pno=caseno/(caseno+contno)
  beta0=-log((1/pno)-1)
  beta1=-log((1/pyes)-1)-beta0
  OR <- exp(beta1)
  pval <- signif(summary(m)$coefficients[2,4],3)
  ### FIRTH
  firth_beta0=coef(firth)[1]
  firth_beta1=coef(firth)[2]
  firth_OR <- exp(firth_beta1)
  firth_pno=1/(exp(-firth_beta0)+1)
  firth_pyes=1/(exp(firth_beta1+firth_beta0)+1) 
  firth_pval <- signif(firth$prob[2],3)
  vec <- c(angi[i,1],caseyes, caseno, contyes, contno, pyes, pno, beta0, beta1, OR, pval, firth_beta0, firth_beta1, firth_OR, firth_pno, firth_pyes, firth_pval)
  res2 <- rbind(res2, vec)
}

colnames(res2) <- c("Index", "caseyes", "caseno", "contyes", "contno", "pyes", "pno", "beta0", "beta1", "OR","pval", "firth_beta0","firth_beta1", "firth_OR", "firth_pno", "firth_pyes", "firth_pval")
write.table(res2, paste(wkdir, "meta_analyses/Logistic_vs_Firth_Dummy.txt", sep="/"), sep ="\t", quote=F, row.names=F, col.names=T)




