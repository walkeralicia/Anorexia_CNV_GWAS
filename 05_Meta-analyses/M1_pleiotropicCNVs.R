### need to change case and control numbers in this script

module load R/3.6
export R_LIBS=/scratch/90days/uqawal15/R_libraries/rlib_3.6
R
library(logistf)
library(readxl)
library(stringr)
library(dplyr)
options(scipen=999)

## Read in UKB and ANGI results 
ukb <- read.table("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/drCNVs_unweighted.txt", header=T)
angi <- read.table("/QRISdata/Q4399/Anorexia/ANGI/results/angi_drCNVs_unweighted.txt", header=T)

## split data by method and reformat
results <- list()
methods <- c("Logistic", "Unadjusted_Logistic", "Firth", "Unadjusted_Firth")
for (i in 1:length(methods)){
  ukb_method <- ukb %>% filter(Method==methods[i])
  angi_method <- angi %>% filter(Method==methods[i])
  comb <- merge(angi_method, ukb_method, by = c("Index"), all=TRUE) #ANGI is x and UKB is Y
  comb$Case_N.y <- ifelse(is.na(comb$Case_N.y), 0, comb$Case_N.y)
  comb$Cont_N.y <- ifelse(is.na(comb$Cont_N.y), 0, comb$Cont_N.y)
  comb$Case_N.x <- ifelse(is.na(comb$Case_N.x), 0, comb$Case_N.x)
  comb$Cont_N.x <- ifelse(is.na(comb$Cont_N.x), 0, comb$Cont_N.x)
  comb$contaunzusa2yes <- ifelse(is.na(comb$contaunzusa2yes), 0, comb$contaunzusa2yes)
  comb$contseyes <- ifelse(is.na(comb$contseyes), 0, comb$contseyes)
  comb$caseaunzusa2yes <- ifelse(is.na(comb$caseaunzusa2yes), 0, comb$caseaunzusa2yes)
  comb$caseseyes <- ifelse(is.na(comb$caseseyes), 0, comb$caseseyes)
  comb$Cont_Freq_per10K.x <- ifelse(is.na(comb$Cont_Freq_per10K.x), 0, comb$Cont_Freq_per10K.x)
  comb$Cont_Freq_per10K.y <- ifelse(is.na(comb$Cont_Freq_per10K.y), 0, comb$Cont_Freq_per10K.y)
  results[[i]] <- as.data.frame(comb)
}


# Counts in each group from ANGI
#Status   group        n
#1     aunzusa2  1373
#1     se        3654
#2    aunzusa2  3725
#2     se        3660

#caseaunzusa2yes = Case + aunzusa2 + CNV
#contaunzusa2yes = Control + aunzusa2 + CNV
#caseseyes = Case + SWE + CNV
#contseyes = Control + SWE + CNV

### Firth Logistic Regression with dummy data and study region as covariate
comb <- results[[4]]
res <- list()
for (i in 1:dim(comb)[1]){
  print(i)
  CNV <- c(rep(1,comb$Case_N.y[i]),rep(0,(1246-(comb$Case_N.y[i]))),
           rep(1,comb$Cont_N.y[i]),rep(0,(381472-(comb$Cont_N.y[i]))), 
           rep(1,comb$caseaunzusa2yes[i]),rep(0,3725-comb$caseaunzusa2yes[i]), 
           rep(1,comb$contaunzusa2yes[i]), rep(0, 1373-comb$contaunzusa2yes[i]),
           rep(1,comb$caseseyes[i]),rep(0,3660-comb$caseseyes[i]), 
           rep(1,comb$contseyes[i]), rep(0, 3654-comb$contseyes[i]))
  AN <- c(rep(1, 1246), 
          rep(0, 381472),
          rep(1,comb$caseaunzusa2yes[i]),rep(1,3725-comb$caseaunzusa2yes[i]), 
          rep(0, comb$contaunzusa2yes[i]), rep(0, 1373-comb$contaunzusa2yes[i]),
          rep(1,comb$caseseyes[i]),rep(1,3660-comb$caseseyes[i]), 
          rep(0, comb$contseyes[i]), rep(0, 3654-comb$contseyes[i]))
  STUDY <- c(rep(2, 382718),
             rep(1,comb$caseaunzusa2yes[i]),rep(1,3725-comb$caseaunzusa2yes[i]), 
             rep(1, comb$contaunzusa2yes[i]), rep(1,1373-comb$contaunzusa2yes[i]),
             rep(0,comb$caseseyes[i]),rep(0,3660-comb$caseseyes[i]), 
             rep(0, comb$contseyes[i]), rep(0, 3654-comb$contseyes[i]))
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
write.table(res_all, "/afm01/UQ/Q4399/Anorexia/UKB/burden_analysis/known_loci/Both_RegionCovariate.txt", sep ="\t", quote=T, row.names=F, col.names=T)


## Comparing FIRTH logistic regression to normal logistic regression (without covariates)

comb <- results[[4]]
res2 <- data.frame()
for (i in 1:dim(comb)[1]){
  print(i)
  CNV <- c(rep(1,comb$Case_N[i]),rep(0,(1246-(comb$Case_N[i]))),
           rep(1,comb$Cont_N[i]),rep(0,(381472-(comb$Cont_N[i]))), 
           rep(1,comb$caseaunzusa2yes[i]),rep(0,3725-comb$caseaunzusa2yes[i]), 
           rep(1,comb$contaunzusa2yes[i]), rep(0, 1373-comb$contaunzusa2yes[i]),
           rep(1,comb$caseseyes[i]),rep(0,3660-comb$caseseyes[i]), 
           rep(1,comb$contseyes[i]), rep(0, 3654-comb$contseyes[i]))
  AN <- c(rep(1, 1246), 
          rep(0, 381472),
          rep(1,comb$caseaunzusa2yes[i]),rep(1,3725-comb$caseaunzusa2yes[i]), 
          rep(0, comb$contaunzusa2yes[i]), rep(0, 1373-comb$contaunzusa2yes[i]),
          rep(1,comb$caseseyes[i]),rep(1,3660-comb$caseseyes[i]), 
          rep(0, comb$contseyes[i]), rep(0, 3654-comb$contseyes[i]))
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
  pyes=caseyes/(caseyes+contyes) # prob of case if you have cnv
  pno=caseno/(caseno+contno) # prob of case if you have dont have cnv
  beta0=-log((1/pno)-1)
  beta1=-log((1/pyes)-1)-beta0
  OR <- exp(beta1)
  pval <- signif(summary(m)$coefficients[2,4],3)
  ### FIRTH
  firth_beta0=coef(firth)[1] # doesnt change to normal logistic regression. 
  firth_beta1=coef(firth)[2] # increases, thus pyes increases and cnv counts for those with AN increases and those without AN decreases. 
  firth_OR <- exp(firth_beta1)
  firth_pno=1/(exp(-firth_beta0)+1) # pno doesnt change. Thus, no cnv counts doesnt change. 
  firth_pyes=1/(exp(firth_beta1+firth_beta0)+1) 
  firth_pval <- signif(firth$prob[2],3)
  vec <- c(angi[i,1],caseyes, caseno, contyes, contno, pyes, pno, beta0, beta1, OR, pval, firth_beta0, firth_beta1, firth_OR, firth_pno, firth_pyes, firth_pval)
  res2 <- rbind(res2, vec)
}

colnames(res2) <- c("Index", "caseyes", "caseno", "contyes", "contno", "pyes", "pno", "beta0", "beta1", "OR","pval", "firth_beta0","firth_beta1", "firth_OR", "firth_pno", "firth_pyes", "firth_pval")
write.table(res2, "/scratch/90days/uqawal15/anx_ukb_cnvs/burden_analyses/AN/known_loci/Logistic_vs_Firth.txt", sep ="\t", quote=F, row.names=F, col.names=T)



##############Pooling effect size and standard errors of the simple analysis via weighted effective sample size##########################

pooled <- results[[2]]
colnames(pooled) <- c("Index", "Syndrome", "Locus", "Type", "Critical_Region", "Coe_significant", "Sz-CNV", "Cognitive_significant",
                      "Method.x", "ANGI_Case_N", "ANGI_Cont_N", "ANGI_Cont_SE", "ANGI_Cont_Freq_10K", "ANGI_Cont_CI", "ANGI_OR", 
                      "ANGI_CI", "ANGI_Pvalue", colnames(all)[18:25], "ANGI_Neff", "Method.y", "UKB_Case_N", "UKB_Cont_N", "UKB_Cont_SE", "UKB_Cont_Freq_10K", 
                      "UKB_Cont_CI", "UKB_OR", "UKB_CI", "UKB_Pvalue", "UKB_Neff")
pooled$Pooled_OR <- signif((pooled$ANGI_OR*pooled$ANGI_Neff+pooled$UKB_OR*pooled$UKB_Neff)/(pooled$ANGI_Neff+pooled$UKB_Neff),3)
pooled$dp <- nchar(pooled$Pooled_OR)-2
pooled$Pooled_SE <- sqrt(1/(pooled$ANGI_Neff+pooled$UKB_Neff))
pooled$Pooled_LCI <- round(pooled$Pooled_OR-1.96*pooled$Pooled_SE, pooled$dp)
pooled$Pooled_UCI <- round(pooled$Pooled_OR+1.96*pooled$Pooled_SE, pooled$dp)
pooled$Pooled_CI <- paste0("(", pooled$Pooled_LCI, "-", pooled$Pooled_UCI, ")")
pooled$Pooled_beta <- log(pooled$Pooled_OR)
pooled$Pooled_TestStatistic <- pooled$Pooled_beta/log(pooled$Pooled_SE)
pooled$Pooled_Pvalue <- 2*pnorm(-(abs(pooled$Pooled_TestStatistic)))

final_pooled <- pooled %>% select(-Index, -Type, -Method.x, -Method.y, -contaunzusa2yes, -contseyes,
                                  -caseaunzusa2yes, -caseseyes, -contaunzusa2no, -contseno,-caseaunzusa2no, -caseseno, -dp)
write.csv(final_pooled, "/scratch/90days/uqawal15/anx_ukb_cnvs/burden_analyses/AN/known_loci/all_pooled.csv",  row.names=F)


###################Meta-analyse the Firth-adjusted regression results with Stouffer's method #############################
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

collins <- read_excel("/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/1-s2.0-S0092867422007887-mmc4.xlsx", sheet=1)
collins <- as.data.frame(collins)
collins$Location <- paste0(collins$Chrom, ":", collins$Start, ":", collins$End)
collins$Genotype <- ifelse(collins$'CNV Type' == "DEL", "Deletion", "Duplication")


comb <- merge(collins, stouffer, by.x = "Segment ID", by.y="Index", all.x=TRUE)
comb$Case_N.y <- ifelse(is.na(comb$Case_N.y), 0, comb$Case_N.y)
comb$Cont_N.y <- ifelse(is.na(comb$Cont_N.y), 0, comb$Cont_N.y)
comb$Case_N.x <- ifelse(is.na(comb$Case_N.x), 0, comb$Case_N.x)
comb$Cont_N.x <- ifelse(is.na(comb$Cont_N.x), 0, comb$Cont_N.x)
comb$contaunzusa2yes <- ifelse(is.na(comb$contaunzusa2yes), 0, comb$contaunzusa2yes)
comb$contseyes <- ifelse(is.na(comb$contseyes), 0, comb$contseyes)
comb$caseaunzusa2yes <- ifelse(is.na(comb$caseaunzusa2yes), 0, comb$caseaunzusa2yes)
comb$caseseyes <- ifelse(is.na(comb$caseseyes), 0, comb$caseseyes)
comb$Cont_Freq_per10K.x <- ifelse(is.na(comb$Cont_Freq_per10K.x), 0, comb$Cont_Freq_per10K.x)
comb$Cont_Freq_per10K.y <- ifelse(is.na(comb$Cont_Freq_per10K.y), 0, comb$Cont_Freq_per10K.y)



res2 <- comb %>% 
  select(Cytoband, Location, Genotype, Case_N.x, Cont_N.x, 
         Cont_SE.x, Cont_Freq_per10K.x, Cont_CI.x, OR.x, CI.x, Zvalue.x, Pvalue.x,
         Case_N.y, Cont_N.y, Cont_SE.y, Cont_Freq_per10K.y, Cont_CI.y, OR.y, CI.y, 
         Pvalue.y, Zvalue.y, Stouffer_Zvalue, Stouffer_Pvalue, Stouffer_onesidedPvalue, OR_meta, CI_meta)


colnames(res2) <- c("Cytoband", "Location", "Genotype", "ANGI_Case_N", "ANGI_Cont_N", "ANGI_Cont_SE", "ANGI_Cont_Freq_10K", "ANGI_Cont_CI", "ANGI_OR", 
                      "ANGI_CI", "ANGI_Pvalue","ANGI_Zvalue", "UKB_Case_N", "UKB_Cont_N", "UKB_Cont_SE", "UKB_Cont_Freq_10K", 
                      "UKB_Cont_CI", "UKB_OR", "UKB_CI", "UKB_Pvalue","UKB_Zvalue", 
                    "Stouffer_Zvalue","Stouffer_Pvalue", "Stouffer_onesidedPvalue", "Meta_OR", "Meta_CI")


RES <- res2 
RES$ANGI_Cont_notN <- 5044-as.numeric(RES$ANGI_Cont_N)
RES$UKB_Cont_notN <- 385930-as.numeric(RES$UKB_Cont_N)
n=length(RES$Cytoband)
for (i in 1:n){
  FM=matrix(c(as.numeric(RES$ANGI_Cont_notN[i]),as.numeric(RES$ANGI_Cont_N[i]),
              as.numeric(RES$UKB_Cont_notN[i]),as.numeric(RES$UKB_Cont_N[i])),nrow=2,ncol=2)
  FT=fisher.test(FM)
  RES$Cont_freq_Pvalue[i]=signif(FT$p.value,2)
}


final <- RES[order(RES$Stouffer_onesidedPvalue),]
final <- final %>% select(-ANGI_Cont_notN, -UKB_Cont_notN)
final$Stouffer_onesidedPvalue <- signif(final$Stouffer_onesidedPvalue, 3)
final$ANGI_Zvalue <- signif(final$ANGI_Zvalue, 3)
final$UKB_Zvalue <- signif(final$UKB_Zvalue, 3)
final$ANGI_Pvalue <- signif(final$ANGI_Pvalue, 3)
final$UKB_Pvalue <- signif(final$UKB_Pvalue, 3)
final$Stouffer_Zvalue <- signif(final$Stouffer_Zvalue, 3)
final$Stouffer_Pvalue <- signif(final$Stouffer_Pvalue, 3)
final$Cont_freq_Pvalue <- signif(final$Cont_freq_Pvalue,3)

write.csv(final, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/collins_all.csv",  row.names=F)

paper <- final %>% filter(Stouffer_onesidedPvalue<0.05)
write.csv(paper, "/QRISdata/Q4399/Anorexia/UKB/burden_analysis/drCNVs/collins_sig.csv",  row.names=F)








