
#============= This script meta-analyses the UKB and ANGI developmental CNV burden results =============================

#== Read in required R libraries =================
library(logistf)
library(readxl)
library(stringr)
library(dplyr)
options(scipen=999)

#== Set paths ===
wkdir<-"/Anorexia/UKB/burden_analysis/devCNVs"
data_path<-"data" ## path to CNV_List.xlsx
ANGI_data<-"ANGI_data" ## path to folder with ANGI meta-analyses results


#== Read in UKB and ANGI results ================
ukb <- read.table(paste(wkdir, "UKBB_devCNVs_burden.txt", sep="/"), header=T)
angi <- read.table(paste(ANGI_data, "ANGI_devCNVs_burden.txt", sep="/"), header=T)

#== Split data by method and reformat ==========
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

# Stouffer's method ==========================

stouffer <- results[[3]]
stouffer$Zvalue.x <- (-1)*sign(log(stouffer$OR.x))*qnorm(stouffer$Pvalue.x/2)
stouffer$Zvalue.y <- (-1)*sign(log(stouffer$OR.y))*qnorm(stouffer$Pvalue.y/2)
stouffer$Stouffer_Zvalue <- (stouffer$Zvalue.x + stouffer$Zvalue.y)/sqrt(2)
stouffer$Stouffer_Pvalue <- 2*pnorm(-(abs(stouffer$Stouffer_Zvalue)))
stouffer$beta_meta <- (stouffer$beta.x/(stouffer$se.x*stouffer$se.x) + stouffer$beta.y/(stouffer$se.y*stouffer$se.y)) / (1/(stouffer$se.x*stouffer$se.x)+1/(stouffer$se.y*stouffer$se.y))
stouffer$se_meta <- 1/sqrt(1/(stouffer$se.x*stouffer$se.x) +1/(stouffer$se.y*stouffer$se.y))
stouffer$OR_meta <- signif(exp(stouffer$beta_meta),3)
stouffer$dp <- nchar(stouffer$OR_meta)-2
stouffer$LCI_meta <- round(exp(stouffer$beta_meta-1.96*stouffer$se_meta), stouffer$dp)
stouffer$HCI_meta <- round(exp(stouffer$beta_meta+1.96*stouffer$se_meta), stouffer$dp)
stouffer$CI_meta <- paste0("(", stouffer$LCI_meta, "-", stouffer$HCI_meta, ")")

# Annotate the information about the pleiotropic CNVs to the result file ========================================================================

cnvs <- read_excel(paste(data_path, "CNV_List.xlsx", sep="/"), sheet=1)
cnvs <- as.data.frame(cnvs)
cnvs$ID <- 1:nrow(cnvs)

# Replace NA values with 0 ======================
comb <- merge(cnvs, stouffer, by.x = "ID", by.y="Index", all.x=TRUE)
cols_to_replace_na <- c("Case_N.y", "Cont_N.y", "Case_N.x", "Cont_N.x", "contaunzusa2yes", "contseyes", "caseaunzusa2yes", "caseseyes", "Cont_Freq_per10K.x", "Cont_Freq_per10K.y")
comb[cols_to_replace_na] <- lapply(comb[cols_to_replace_na], function(x) ifelse(is.na(x), 0, x))


res <- comb %>% 
  select(ID, CNV, Location, Genotype, Size, pHaplo, pTtriplo, Rees, Marshall,
         ANGI_Case_N := Case_N.x, ANGI_Cont_N := Cont_N.x, 
         ANGI_Cont_SE := Cont_SE.x, ANGI_Cont_Freq_10K := Cont_Freq_per10K.x, 
         ANGI_Cont_CI := Cont_CI.x, ANGI_OR := OR.x, ANGI_CI := CI.x, 
         ANGI_Zvalue := Zvalue.x, ANGI_Pvalue := Pvalue.x,
         UKB_Case_N := Case_N.y, UKB_Cont_N := Cont_N.y, 
         UKB_Cont_SE := Cont_SE.y, UKB_Cont_Freq_10K := Cont_Freq_per10K.y, 
         UKB_Cont_CI := Cont_CI.y, UKB_OR := OR.y, UKB_CI := CI.y, 
         UKB_Zvalue := Zvalue.y, UKB_Pvalue := Pvalue.y,
         Stouffer_Zvalue, Stouffer_Pvalue, Meta_OR := OR_meta, Meta_CI := CI_meta)



# Fisher's exact test to compare CNV control frequency between ANGI and the UKB ===========
RES <- res 
ukb_controls = 385930 # UKB Control Sample Size
angi_controls = 5044 ## ANGI Control Sample Size
RES$ANGI_Cont_notN <- angi_controls-as.numeric(RES$ANGI_Cont_N)
RES$UKB_Cont_notN <- ukb_controls-as.numeric(RES$UKB_Cont_N)
n=length(RES$ID)
for (i in 1:n){
  FM=matrix(c(as.numeric(RES$ANGI_Cont_notN[i]),as.numeric(RES$ANGI_Cont_N[i]),
              as.numeric(RES$UKB_Cont_notN[i]),as.numeric(RES$UKB_Cont_N[i])),nrow=2,ncol=2)
  FT=fisher.test(FM)
  RES$Cont_freq_Pvalue[i]=signif(FT$p.value,2)
}

# Reformating and saving final results ================================================

final <- RES %>%
  arrange(Stouffer_Pvalue) %>%
  select(-ANGI_Cont_notN, -UKB_Cont_notN) %>%
  mutate(across(c(Stouffer_Pvalue, ANGI_Zvalue, UKB_Zvalue, ANGI_Pvalue, UKB_Pvalue, Stouffer_Zvalue, Stouffer_Pvalue, Cont_freq_Pvalue), ~signif(., 3)))
final$Rees <- ifelse(is.na(final$Rees), "no", final$Rees)
final$Marshall <- ifelse(is.na(final$Marshall), "no", final$Marshall)


write.csv(final, paste(wkdir, "meta_analyses/devCNVs_all.csv", sep="/"),  row.names=F)
paper <- final %>% filter(Stouffer_Pvalue<0.05)
write.csv(paper, paste(wkdir, "meta_analyses/devCNVs_sig.csv", sep="/"),  row.names=F)

